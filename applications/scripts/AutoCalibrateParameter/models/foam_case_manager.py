"""
OpenFOAM Case Manager

Manages the running and post-processing of OpenFOAM cases.
Migrated from the original meltpool_model.py.
"""

from __future__ import annotations

import math
import os
import re
import shlex
import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd


# ========== Utility functions ==========


def _run(cmd: str, workdir: Path) -> None:
    """Execute a shell command"""
    print(f"[cmd] {cmd}")
    subprocess.run(cmd, cwd=str(workdir), shell=True, check=True)


def _run_script(script: str, workdir: Path, runner: Optional[str]) -> None:
    """Execute a script with optional custom runner"""
    script = script.strip() + ("\n" if not script.endswith("\n") else "")
    label = runner or "bash -lc"
    print(f"[cmd:{label}] <<'EOF'\n{script}EOF")

    if runner:
        subprocess.run(
            runner, cwd=str(workdir), shell=True, check=True, text=True, input=script
        )
    else:
        subprocess.run(["bash", "-lc", script], cwd=str(workdir), check=True)


def _update_dict_value(path: Path, key: str, value: float) -> None:
    """Update a value in an OpenFOAM dictionary file"""
    lines = path.read_text().splitlines()
    pattern = re.compile(rf"^\s*{re.escape(key)}\b")

    for i, line in enumerate(lines):
        if pattern.match(line):
            indent = line[: len(line) - len(line.lstrip())]
            lines[i] = f"{indent}{key:<23}{value};"
            path.write_text("\n".join(lines) + "\n")
            return

    lines.append(f"{key:<23}{value};")
    path.write_text("\n".join(lines) + "\n")


def _update_T_internal(path: Path, temperature: float) -> None:
    """Update the initial value of the temperature field"""
    lines = path.read_text().splitlines()
    pattern = re.compile(r"^\s*internalField\s+uniform\s+([0-9eE+\-\.]+)")

    for i, line in enumerate(lines):
        if pattern.search(line):
            indent = line[: len(line) - len(line.lstrip())]
            lines[i] = f"{indent}internalField   uniform {temperature};"
            path.write_text("\n".join(lines) + "\n")
            return
    raise RuntimeError(f"internalField line not found in {path}")


def _rewrite_power_file(path: Path, power_w: float) -> None:
    """Modify the laser power file"""
    lines = path.read_text().splitlines()
    pattern = re.compile(r"\(\s*([0-9eE+\-\.]+)\s+([0-9eE+\-\.]+)\s*\)")

    new_lines = []
    for line in lines:
        m = pattern.search(line)
        if m:
            t_str, val_str = m.group(1), m.group(2)
            new_val = power_w if float(val_str) != 0.0 else 0.0
            new_lines.append(f"    ({t_str}           {new_val})")
        else:
            new_lines.append(line)
    path.write_text("\n".join(new_lines) + "\n")


# ========== OpenFOAM Case Manager ==========


class OpenFOAMCaseManager:
    """
    Manages the running and post-processing of OpenFOAM cases

    Attributes
    ----------
    case_dir : Path
        Case directory
    postproc_script : Path
        Post-processing script path
    foam_bashrc : str, optional
        OpenFOAM bashrc path
    foam_runner : str, optional
        OpenFOAM runner command
    n_proc : int
        Number of MPI parallel processes
    hpc_mode : bool
        Whether HPC mode is enabled
    """

    def __init__(
        self,
        case_dir: Path,
        postproc_script: Path,
        foam_bashrc: Optional[str] = None,
        postproc_python: str = "python",
        postproc_runner: Optional[str] = None,
        pvpython: Optional[str] = None,
        n_proc: int = 12,
        foam_runner: Optional[str] = None,
        mpirun_flags: Tuple[str, ...] = ("--oversubscribe",),
        hpc_mode: bool = False,
        laser_diameter: float = 70e-6,
        cell_size: float = 1e-6,
        x_domain: Tuple[float, float] = (0.0, 0.0003),
        y_begin_track: float = 150e-6,
        y_end_track: float = 300e-6,
        plot_geometry: bool = True,
        archive_mode: str = "latest",
    ):
        """
        Initialize the case manager

        Parameters
        ----------
        case_dir : Path
            OpenFOAM case directory
        postproc_script : Path
            Post-processing script path
        foam_bashrc : str, optional
            OpenFOAM bashrc path
        postproc_python : str
            Post-processing Python interpreter
        postproc_runner : str, optional
            Post-processing runner
        pvpython : str, optional
            Path to pvpython executable
        n_proc : int
            Number of MPI parallel processes
        foam_runner : str, optional
            OpenFOAM runner command (e.g. "of2506")
        mpirun_flags : tuple
            Extra mpirun arguments
        hpc_mode : bool
            Whether HPC mode is enabled
        """
        self.case_dir = Path(case_dir).resolve()
        self.postproc_script = Path(postproc_script).expanduser().resolve()
        self.foam_bashrc = foam_bashrc
        self.postproc_python = postproc_python
        self.postproc_runner = postproc_runner
        self.pvpython = pvpython
        self.n_proc = n_proc
        self.foam_runner = foam_runner
        self.mpirun_flags = tuple(mpirun_flags)
        self.hpc_mode = hpc_mode
        self.laser_diameter = laser_diameter
        self.cell_size = cell_size
        self.x_domain = tuple(x_domain)
        self.y_begin_track = y_begin_track
        self.y_end_track = y_end_track
        self.plot_geometry = plot_geometry

        normalized_archive_mode = (
            "latest" if archive_mode is None else str(archive_mode).strip().lower()
        )
        if normalized_archive_mode not in {"latest", "all"}:
            raise ValueError("archive_mode must be 'latest' or 'all'")
        self.archive_mode = normalized_archive_mode

        # Key file paths
        self.transport_props = self.case_dir / "constant" / "transportProperties"
        self.temp_field = self.case_dir / "initial" / "T"
        self.time_vs_power = self.case_dir / "constant" / "timeVsLaserPower"
        self.laser_props = self.case_dir / "constant" / "LaserProperties"
        self.main_foam = self.case_dir / "main.foam"

    def _build_reconstruct_command(self) -> str:
        archive_mode = str(getattr(self, "archive_mode", "latest")).strip().lower()
        if archive_mode == "all":
            return "reconstructPar -noZero >log.reconstructPar 2>&1"
        return "reconstructPar -latestTime >log.reconstructPar 2>&1"

    def update_recoil_coeff(self, value: float) -> None:
        """
        Update the recoil pressure coefficient in transportProperties

        Parameters
        ----------
        value : float
            Recoil pressure coefficient
        """
        _update_dict_value(self.transport_props, "recoilCoeff", value)

    def update_parameters(
        self,
        sigma: float,
        marangoni: float,
        substrate_temp: float,
        absorptivity: float = 1.0,
        recoil_coeff: float = 1.0,
        radius_flavour: float = 2.0,
        laser_radius: float = 50e-6,
    ) -> None:
        """
        Update material parameters

        Parameters
        ----------
        sigma : float
            Surface tension coefficient
        marangoni : float
            Marangoni constant
        substrate_temp : float
            Substrate temperature (K)
        absorptivity : float, optional
            Laser absorptivity coefficient (default: 1.0)
        recoil_coeff : float, optional
            Recoil pressure coefficient (default: 1.0)
        radius_flavour : float, optional
            Radius_Flavour in LaserProperties (default: 2.0)
        laser_radius : float, optional
            laserRadius in LaserProperties, in meters (default: 50e-6)
        """
        _update_dict_value(self.transport_props, "sigma", sigma)
        _update_dict_value(self.transport_props, "Marangoni_Constant", marangoni)
        _update_T_internal(self.temp_field, substrate_temp)
        _update_dict_value(self.laser_props, "absorptivity", absorptivity)
        _update_dict_value(self.transport_props, "recoilCoeff", recoil_coeff)
        _update_dict_value(self.laser_props, "Radius_Flavour", radius_flavour)
        _update_dict_value(self.laser_props, "laserRadius", laser_radius)

    def set_power(self, power_w: float, absorptivity: float = 1.0) -> None:
        """
        Set the laser power (accounting for absorptivity)

        Effective power = nominal power x absorptivity

        Parameters
        ----------
        power_w : float
            Nominal laser power (W)
        absorptivity : float
            Absorptivity coefficient
        """
        effective_power = power_w * absorptivity
        _rewrite_power_file(self.time_vs_power, effective_power)

    def run_simulation(self) -> None:
        """Run the OpenFOAM simulation"""
        if self.hpc_mode:
            self._run_simulation_hpc()
        else:
            self._run_simulation_pc()

    def _run_simulation_hpc(self) -> None:
        """Run simulation in HPC mode"""
        mpirun_cmd = [
            "mpirun",
            *self.mpirun_flags,
            "-np",
            str(self.n_proc),
            "laserbeamFoam",
            "-parallel",
            ">log.laserbeamFoam",
            "2>&1",
        ]

        lines = ["set -eo pipefail"]
        if self.foam_bashrc:
            lines.append(f"source {shlex.quote(self.foam_bashrc)}")
        lines.extend(
            [
                "bash ./Allclean || true",
                "rm -rf 0 processor*",
                "cp -r initial 0",
                "touch main.foam",
                "blockMesh >log.blockMesh 2>&1",
                "setSolidFraction >log.setSolidFraction 2>&1",
                "decomposePar >log.decomposePar 2>&1",
                " ".join(mpirun_cmd),
                self._build_reconstruct_command(),
            ]
        )
        _run_script("\n".join(lines), self.case_dir, self.foam_runner)

    def _run_simulation_pc(self) -> None:
        """Run simulation in PC (personal computer) mode"""
        DEFAULT_BASHRC = "/usr/lib/openfoam/openfoam2506/etc/bashrc"
        bashrc = shlex.quote(self.foam_bashrc or DEFAULT_BASHRC)
        case_dir = shlex.quote(str(self.case_dir))
        mpirun_flags = (
            " ".join(self.mpirun_flags) if self.mpirun_flags else "--oversubscribe"
        )

        reconstruct_cmd = self._build_reconstruct_command()

        cmd = (
            f"bash -lc 'set -eo pipefail; "
            f"export WM_PROJECT_SITE=${{WM_PROJECT_SITE-}}; source {bashrc} && "
            f"cd {case_dir} && bash ./Allclean || true; rm -rf 0 processor*; "
            f"cp -r initial 0; touch main.foam; "
            f"blockMesh >log.blockMesh 2>&1; "
            f"setSolidFraction >log.setSolidFraction 2>&1; "
            f"decomposePar >log.decomposePar 2>&1; "
            f"mpirun -np {self.n_proc} {mpirun_flags} laserbeamFoam -parallel >log.laserbeamFoam 2>&1; "
            f"{reconstruct_cmd}'"
        )
        _run(cmd, self.case_dir)

    def run_postprocess(self) -> None:
        """Run the post-processing script"""
        self._write_input_data()
        if self.hpc_mode:
            self._run_postprocess_hpc()
        else:
            self._run_postprocess_pc()

    def _write_input_data(self) -> None:
        """Generate input_data.py in the case directory, overriding defaults in the post-processing script"""
        path = self.case_dir / "input_data.py"
        of_loc = self.foam_bashrc or ""
        path.write_text(
            f"LASER_DIAMETER = {self.laser_diameter!r}\n"
            f"OF_LOCATION = {of_loc!r}\n"
            f"CELL_SIZE = {self.cell_size!r}\n"
            f"X_MIN_AND_MAX_DOMAIN = [{self.x_domain[0]!r}, {self.x_domain[1]!r}]\n"
            f"Y_COORD_BEGIN_TRACK = {self.y_begin_track!r}\n"
            f"Y_COORD_END_TRACK = {self.y_end_track!r}\n"
            f"PLOT_GEOMETRY_VS_Y_LOCATION = {self.plot_geometry!r}\n"
        )

    def _run_postprocess_hpc(self) -> None:
        """Run post-processing in HPC mode"""
        if not self.main_foam.exists():
            self.main_foam.touch()

        env = os.environ.copy()
        if self.postproc_runner:
            env["POSTPROC_RUNNER"] = self.postproc_runner
        if self.pvpython:
            env["PVPYTHON"] = self.pvpython

        cmd = shlex.split(self.postproc_python) + [str(self.postproc_script)]
        print(f"[cmd] {' '.join(shlex.quote(c) for c in cmd)}")
        subprocess.run(cmd, cwd=str(self.case_dir), check=True, env=env)

    def _run_postprocess_pc(self) -> None:
        """Run post-processing in PC mode"""
        _run(
            f"{self.postproc_python} {shlex.quote(str(self.postproc_script))}",
            self.case_dir,
        )

    def read_metrics(self) -> Dict[str, float]:
        """
        Read post-processing results

        Returns
        -------
        dict
            Dictionary containing width_mean_m, depth_mean_m, area_mean_m2
        """
        metrics_path = self.case_dir / "cross_sections_statistics.csv"
        if not metrics_path.exists():
            raise FileNotFoundError(f"Missing {metrics_path}; did post-processing succeed?")

        df = pd.read_csv(metrics_path)
        return {
            "width_mean_m": float(df["width"].mean()),
            "depth_mean_m": float(df["depth"].mean()),
            "area_mean_m2": float(df["area"].mean()) if "area" in df.columns else math.nan,
        }

    def _get_numeric_time_dirs(self) -> list[Path]:
        """Get and sort the numeric time-step directories under the case directory."""
        time_dirs: list[Path] = []
        for path in self.case_dir.iterdir():
            if not path.is_dir():
                continue
            try:
                float(path.name)
                time_dirs.append(path)
            except ValueError:
                continue
        return sorted(time_dirs, key=lambda p: float(p.name))

    def archive_latest_time(self, dest_parent: Path, new_name: Optional[str] = None) -> None:
        """
        Archive the last time-step directory

        Parameters
        ----------
        dest_parent : Path
            Destination parent directory
        new_name : str, optional
            New directory name; if None, keep original name
        """
        import shutil

        # Find the latest numeric directory
        time_dirs = self._get_numeric_time_dirs()

        if not time_dirs:
            print("  [Archive] No time directories found")
            return

        # Sort by numeric value
        latest_dir = sorted(time_dirs, key=lambda p: float(p.name))[-1]

        dest_name = new_name if new_name else latest_dir.name
        dest_path = dest_parent / dest_name

        if dest_path.exists():
            shutil.rmtree(dest_path)

        try:
            shutil.copytree(latest_dir, dest_path)
            print(f"  [Archive] Saved {latest_dir.name} -> {dest_path}")
        except Exception as e:
            print(f"  [Archive] Save failed: {e}")

    def archive_all_results(self, dest_parent: Path) -> None:
        """
        Archive all results (all time steps + key output files).

        Parameters
        ----------
        dest_parent : Path
            Destination directory
        """
        import shutil

        dest_parent = Path(dest_parent)
        dest_parent.mkdir(parents=True, exist_ok=True)

        try:
            # 1) All numeric time steps
            time_dirs = self._get_numeric_time_dirs()
            if not time_dirs:
                print("  [Archive] No time directories found")
            for src in time_dirs:
                dst = dest_parent / src.name
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)

            # 2) Key directories
            for dirname in ("constant", "system"):
                src = self.case_dir / dirname
                if not src.is_dir():
                    continue
                dst = dest_parent / dirname
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)

            # 3) Logs and result CSVs
            for pattern in ("log.*", "*.csv"):
                for src in self.case_dir.glob(pattern):
                    if src.is_file():
                        shutil.copy2(src, dest_parent)

            # 4) Additional files for visualization/reproduction
            for filename in ("main.foam", "input_data.py"):
                src = self.case_dir / filename
                if src.is_file():
                    shutil.copy2(src, dest_parent)

            print(f"  [Archive] All results saved -> {dest_parent}")
        except Exception as e:
            print(f"  [Archive] Save failed: {e}")

    def archive_results(self, dest_parent: Path, mode: str = "latest") -> None:
        """
        Save results according to the archive mode.

        Parameters
        ----------
        dest_parent : Path
            Destination directory
        mode : str
            Archive mode: "latest" or "all"
        """
        normalized = str(mode).strip().lower()
        if normalized == "all":
            self.archive_all_results(dest_parent)
            return
        if normalized == "latest":
            self.archive_latest_time(dest_parent)
            return
        raise ValueError(f"Unknown archive mode: {mode}")

    @classmethod
    def from_config(cls, config) -> "OpenFOAMCaseManager":
        """
        Create a case manager from a configuration object

        Parameters
        ----------
        config : BaseConfig
            Configuration object

        Returns
        -------
        OpenFOAMCaseManager
        """
        return cls(
            case_dir=config.case_dir,
            postproc_script=config.postproc_script,
            foam_bashrc=config.foam_bashrc,
            postproc_python=config.postproc_python,
            postproc_runner=getattr(config, "postproc_runner", None),
            pvpython=config.pvpython,
            n_proc=config.n_proc,
            foam_runner=config.foam_runner,
            mpirun_flags=config.mpirun_flags,
            hpc_mode=config.hpc_mode,
            laser_diameter=config.laser_diameter,
            cell_size=config.cell_size,
            x_domain=config.x_domain,
            y_begin_track=config.y_begin_track,
            y_end_track=config.y_end_track,
            plot_geometry=config.plot_geometry,
            archive_mode=getattr(config, "archive_mode", "latest"),
        )
