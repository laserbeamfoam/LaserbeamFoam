#!/usr/bin/env python3
"""Generate phase1 beam-shape simulation cases for LaserbeamFoam."""

import os
import shutil
import math

WORKBOOK = os.path.dirname(os.path.abspath(__file__)) + "/workbook"
OUTPUT_ROOT = os.path.dirname(os.path.abspath(__file__)) + "/phase1"
PHASE1_SET = f"{WORKBOOK}/Glumann/data/adresses/phase1/phase1_set.txt"
ADR_DIR = f"{WORKBOOK}/Glumann/data/adresses/koord_doe"
TEMPLATE_DIR = f"{WORKBOOK}/IFSW/opa_cbc_movingFrame_ss316L"
SOLVER_DIR = f"{WORKBOOK}/IFSW/opa_cbc_movingFrame_6060"

POWER_W = 3300  # 3.3 kW in W
VELOCITY_MS = 0.2  # 0.2 m/s
VELOCITY_STR = "0200ms"  # naming convention for 0.2 m/s

FREQUENCIES = [
    ("1kHz", 1e3),
    ("10kHz", 10e3),
    ("100kHz", 100e3),
]

# Map from phase1_set row index (0-based) to address file number
# Row 0 -> 001, Row 1 -> 004, ..., Row i -> 1 + i*3
def addr_file_number(row_idx):
    return 1 + row_idx * 3


def parse_phase1_set(path):
    shapes = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            # column 3 (index 2) has the shape name with .png
            shape_png = parts[2].strip()
            shape_name = shape_png.replace(".png", "")
            shapes.append(shape_name)
    return shapes


def read_address_file(path):
    coords = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # format: "x,y"
            x_str, y_str = line.split(",")
            x = int(x_str.strip())
            y = int(y_str.strip())
            coords.append((x, y))
    return coords


def make_beamShapeDict(shape_name, coords, freq_hz):
    n_addr = len(coords)
    duration = 1.0 / freq_hz  # shapeDuration = 1/shapeFrequency
    addr_lines = "\n".join(f"    ( {x}  {y}  )" for x, y in coords)
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  6
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      shapeNamesDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

shapeNames
(
    "{shape_name}"
);

shapeFrequency    ({freq_hz});       // Hz
shapeDuration     ({duration:.15e});  // s

{shape_name}
(
{addr_lines}
);
"""


def make_timeVsLaserPosition(velocity_ms, end_time=1.0):
    # scan in -z direction at given velocity
    dist = velocity_ms * end_time
    return f"""(
    (0.0    (0 -2e-3 0.0))
    ({end_time}    (0 -2e-3 {-dist}))
)
"""


def make_timeVsLaserPower(power_w, end_time=1.0):
    return f"""(
    (0.0      {power_w})
    ({end_time}      {power_w})
)
"""


def write_file(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(content)


def copy_file(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)


def create_initial_files(case_dir, velocity_ms):
    """Create initial/ directory with boundary conditions from ss316L template,
    adjusting velocity to the required value."""
    src_initial = f"{TEMPLATE_DIR}/initial"
    dst_initial = f"{case_dir}/initial"
    if not os.path.exists(dst_initial):
        shutil.copytree(src_initial, dst_initial)

    # Update U with correct velocity (negative z-direction)
    u_path = f"{dst_initial}/U"
    vel_str = f"{velocity_ms:.15g}"
    with open(u_path) as f:
        content = f.read()
    # Replace velocity values
    content = content.replace("0.166667", vel_str).replace("-0.166667", f"-{vel_str}")
    with open(u_path, "w") as f:
        f.write(content)


def create_case(shape_name, freq_label, freq_hz):
    """Create a complete case directory for one shape/frequency combination."""
    # Determine address file
    shapes = parse_phase1_set(PHASE1_SET)
    try:
        row_idx = shapes.index(shape_name)
    except ValueError:
        print(f"  WARNING: shape {shape_name} not found in phase1_set.txt, skipping")
        return

    file_num = addr_file_number(row_idx)
    addr_file = f"{ADR_DIR}/{file_num:03d}_{shape_name}.txt"
    if not os.path.exists(addr_file):
        print(f"  WARNING: address file not found: {addr_file}, skipping")
        return

    coords = read_address_file(addr_file)
    print(f"  Shape {shape_name}: {len(coords)} addresses from {addr_file}")

    # Case directory name
    case_name = f"{shape_name}_{freq_label}_7500W_{VELOCITY_STR}"
    case_dir = f"{OUTPUT_ROOT}/{case_name}"

    # Create directory structure
    os.makedirs(case_dir, exist_ok=True)

    # --- system/ (from 6060 solver config) ---
    sys_src = f"{SOLVER_DIR}/system"
    sys_dst = f"{case_dir}/system"
    if not os.path.exists(sys_dst):
        shutil.copytree(sys_src, sys_dst)

    # --- constant/ (mix of 6060 solver config + custom files) ---
    const_src = f"{SOLVER_DIR}/constant"
    const_dst = f"{case_dir}/constant"
    if not os.path.exists(const_dst):
        shutil.copytree(const_src, const_dst)

    # Override transportProperties with ss316L material (from template)
    tp_src = f"{TEMPLATE_DIR}/constant/transportProperties"
    tp_dst = f"{case_dir}/constant/transportProperties"
    copy_file(tp_src, tp_dst)

    # Write custom beamShapeDict
    beam_dict = make_beamShapeDict(shape_name, coords, freq_hz)
    write_file(f"{case_dir}/constant/beamShapeDict", beam_dict)

    # Write custom timeVsLaserPosition0
    pos_content = make_timeVsLaserPosition(VELOCITY_MS)
    write_file(f"{case_dir}/constant/timeVsLaserPosition0", pos_content)

    # Write custom timeVsLaserPower0
    power_content = make_timeVsLaserPower(POWER_W)
    write_file(f"{case_dir}/constant/timeVsLaserPower0", power_content)

    # --- initial/ (from ss316L template with adjusted velocity) ---
    create_initial_files(case_dir, VELOCITY_MS)

    # --- Allrun and Allclean from template ---
    copy_file(f"{SOLVER_DIR}/Allrun", f"{case_dir}/Allrun")
    copy_file(f"{TEMPLATE_DIR}/Allclean", f"{case_dir}/Allclean")
    os.chmod(f"{case_dir}/Allrun", 0o755)
    os.chmod(f"{case_dir}/Allclean", 0o755)

    print(f"  Created: {case_dir}")


def main():
    shapes = parse_phase1_set(PHASE1_SET)
    print(f"Found {len(shapes)} shapes in phase1_set.txt:")
    for s in shapes:
        print(f"  - {s}")

    print(f"\nGenerating {len(shapes) * len(FREQUENCIES)} cases...\n")

    for shape_name in shapes:
        for freq_label, freq_hz in FREQUENCIES:
            print(f"\nShape: {shape_name}, Frequency: {freq_label} ({freq_hz} Hz)")
            create_case(shape_name, freq_label, freq_hz)

    print(f"\nDone! Cases created under {OUTPUT_ROOT}/")
    print(f"\nSummary: {len(shapes)} shapes x {len(FREQUENCIES)} frequencies = {len(shapes) * len(FREQUENCIES)} cases")


if __name__ == "__main__":
    main()
