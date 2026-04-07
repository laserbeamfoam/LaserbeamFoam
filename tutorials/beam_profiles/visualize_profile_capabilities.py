#!/usr/bin/env python3
"""Generate an HTML gallery that visualizes beam profile capabilities.

The gallery mirrors the profile and modifier mechanics implemented in
`src/laserHeatSource/HeatSourceProfiles.C`:

- sampling strategies: radial-polar and Cartesian
- profiles: superGaussian, annular, bessel, csv, and imported beam-shape data
- modifiers: planeBias, drumNoise, cartesianModeNoise, and gridNoise
- weighted profile combinations
- bare noise fields and the same noise fields applied to beam profiles

The script is self-contained and uses only the Python standard library.
It writes an HTML index plus a set of SVG panels into an output directory.
"""

from __future__ import annotations

import argparse
import html
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "_generated_profile_capabilities"
SMALL = 1.0e-12
TAU = 2.0 * math.pi


@dataclass
class Panel:
    title: str
    caption: str
    filename: str


def clamp(value: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, value))


def lerp(a: float, b: float, t: float) -> float:
    return a + t * (b - a)


def blend_factor(x: float) -> float:
    clamped = clamp(x, 0.0, 1.0)
    return clamped * clamped * (3.0 - 2.0 * clamped)


def rotate(x: float, z: float, degrees: float) -> tuple[float, float]:
    angle = math.radians(degrees)
    c = math.cos(angle)
    s = math.sin(angle)
    return c * x + s * z, -s * x + c * z


def estimate_bessel_zero(azimuthal_order: int, radial_order: int) -> float:
    known = (
        (2.4048255577, 5.5200781103, 8.6537279129, 11.7915344391),
        (3.8317059702, 7.0155866698, 10.1734681351, 13.3236919363),
        (5.1356223018, 8.4172441404, 11.6198411721, 14.7959517824),
        (6.3801618952, 9.7610231299, 13.0152007217, 16.2234661603),
    )
    if 0 <= azimuthal_order < 4 and 1 <= radial_order <= 4:
        return known[azimuthal_order][radial_order - 1]
    return math.pi * (radial_order + 0.5 * azimuthal_order - 0.25)


def bessel_j(order: int, value: float, terms: int = 80) -> float:
    """Integer-order Bessel J_n using the power-series definition."""
    half = 0.5 * value
    total = 0.0
    for m in range(terms):
        sign = -1.0 if m % 2 else 1.0
        numer = half ** (2 * m + order)
        denom = math.factorial(m) * math.factorial(m + order)
        term = sign * numer / max(denom, 1.0)
        total += term
        if abs(term) < 1.0e-14:
            break
    return total


def max_abs(values: list[list[float]]) -> float:
    maximum = 0.0
    for row in values:
        for value in row:
            maximum = max(maximum, abs(value))
    return maximum


def grid_max(values: list[list[float]]) -> float:
    maximum = 0.0
    for row in values:
        for value in row:
            maximum = max(maximum, value)
    return maximum


def normalize_integral(
    values: list[list[float]],
    dx: float,
    dz: float,
) -> list[list[float]]:
    integral = 0.0
    for row in values:
        for value in row:
            integral += max(value, 0.0) * dx * dz

    if integral <= SMALL:
        return [[0.0 for _ in row] for row in values]

    scale = 1.0 / integral
    return [[max(value, 0.0) * scale for value in row] for row in values]


def combine_weighted_profiles(
    profiles: list[list[list[float]]],
    weights: list[float],
    dx: float,
    dz: float,
) -> list[list[float]]:
    normalized = [normalize_integral(profile, dx, dz) for profile in profiles]
    positive = [max(weight, 0.0) for weight in weights]
    weight_sum = sum(positive)

    if weight_sum <= SMALL:
        return [[0.0 for _ in row] for row in normalized[0]]

    nz = len(normalized[0])
    nx = len(normalized[0][0])
    result = [[0.0 for _ in range(nx)] for _ in range(nz)]

    for profile, weight in zip(normalized, positive):
        factor = weight / weight_sum
        for iz in range(nz):
            for ix in range(nx):
                result[iz][ix] += factor * profile[iz][ix]
    return result


def sample_grid(
    func: Callable[[float, float], float],
    half_width_x: float,
    half_width_z: float,
    nx: int = 81,
    nz: int = 81,
) -> tuple[list[list[float]], float, float]:
    dx = 2.0 * half_width_x / max(nx - 1, 1)
    dz = 2.0 * half_width_z / max(nz - 1, 1)
    grid: list[list[float]] = []
    for iz in range(nz):
        z = -half_width_z + iz * dz
        row: list[float] = []
        for ix in range(nx):
            x = -half_width_x + ix * dx
            row.append(func(x, z))
        grid.append(row)
    return grid, dx, dz


def read_scalar_matrix(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    with path.open("r", encoding="ascii", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            rows.append([float(value) for value in stripped.replace(",", " ").split()])
    if not rows:
        raise ValueError(f"{path} does not contain any numeric rows")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError(f"{path} contains ragged rows")
    return rows


def read_flat_scalar_grid(path: Path, nx: int, nz: int) -> list[list[float]]:
    values: list[float] = []
    with path.open("r", encoding="ascii", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            values.extend(float(value) for value in stripped.replace(",", " ").split())

    expected = nx * nz
    if len(values) != expected:
        raise ValueError(f"{path} contains {len(values)} values, expected {expected}")

    return [values[row * nx : (row + 1) * nx] for row in range(nz)]


def bilinear_sample(values: list[list[float]], fx: float, fz: float) -> float:
    nz = len(values)
    nx = len(values[0])
    fx = clamp(fx, 0.0, nx - 1.0)
    fz = clamp(fz, 0.0, nz - 1.0)

    ix0 = min(int(fx), nx - 1)
    iz0 = min(int(fz), nz - 1)
    ix1 = min(ix0 + 1, nx - 1)
    iz1 = min(iz0 + 1, nz - 1)
    tx = fx - ix0
    tz = fz - iz0

    v00 = values[iz0][ix0]
    v10 = values[iz0][ix1]
    v01 = values[iz1][ix0]
    v11 = values[iz1][ix1]
    interp_x0 = v00 + tx * (v10 - v00)
    interp_x1 = v01 + tx * (v11 - v01)
    return interp_x0 + tz * (interp_x1 - interp_x0)


def infer_square_shape(flat_count: int) -> int:
    side = int(round(math.sqrt(flat_count)))
    if side * side != flat_count:
        raise ValueError(f"cannot infer square shape from {flat_count} values")
    return side


def read_beamshape_matrix(path: Path) -> list[list[float]]:
    tokens: list[float] = []
    with path.open("r", encoding="ascii", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            tokens.extend(float(value) for value in stripped.split())
    side = infer_square_shape(len(tokens))
    return [tokens[row * side : (row + 1) * side] for row in range(side)]


def super_gaussian_profile(
    radius_x: float,
    radius_z: float,
    order: float,
    cutoff: float = 1.0e-6,
    rotation_deg: float = 0.0,
    center_shift: tuple[float, float] = (0.0, 0.0),
) -> Callable[[float, float], float]:
    cutoff = clamp(cutoff, 1.0e-12, 0.999999)
    major_radius = max(radius_x, SMALL)
    minor_radius = max(radius_z, SMALL)
    angle = math.radians(rotation_deg)
    c = math.cos(angle)
    s = math.sin(angle)
    inv_major2 = 1.0 / (major_radius * major_radius)
    inv_minor2 = 1.0 / (minor_radius * minor_radius)
    q_xx = c * c * inv_major2 + s * s * inv_minor2
    q_xz = c * s * (inv_major2 - inv_minor2)
    q_zz = s * s * inv_major2 + c * c * inv_minor2

    def evaluate(x: float, z: float) -> float:
        x_local = x - center_shift[0]
        z_local = z - center_shift[1]
        quadratic = q_xx * x_local * x_local + 2.0 * q_xz * x_local * z_local + q_zz * z_local * z_local
        radius = math.sqrt(max(quadratic, 0.0))
        return max(math.exp(-radius**order), 0.0)

    return evaluate


def annular_profile(
    inner_radius: float,
    outer_radius: float,
    falloff: float,
    rotation_deg: float = 0.0,
    peak_at_inner: bool = False,
) -> Callable[[float, float], float]:
    falloff = max(falloff, SMALL)

    def evaluate(x: float, z: float) -> float:
        x_local, z_local = rotate(x, z, rotation_deg)
        r = math.sqrt(x_local * x_local + z_local * z_local)
        inner_transition = 0.5 * (1.0 + math.tanh((r - inner_radius) / falloff))
        outer_transition = 0.5 * (1.0 - math.tanh((r - outer_radius) / falloff))
        intensity = inner_transition * outer_transition
        if peak_at_inner:
            intensity *= inner_radius / max(r, SMALL)
        return max(intensity, 0.0)

    return evaluate


def bessel_profile(
    order: int,
    first_zero_radius: float,
    support_radius: float,
) -> Callable[[float, float], float]:
    estimated_first_zero = estimate_bessel_zero(order, 1)
    radius = first_zero_radius / max(estimated_first_zero, SMALL)

    def evaluate(x: float, z: float) -> float:
        r = math.sqrt(x * x + z * z)
        if r > support_radius:
            return 0.0
        value = bessel_j(order, r / radius)
        return max(value * value, 0.0)

    return evaluate


def csv_profile(
    values: list[list[float]],
    dx: float,
    dz: float,
    half_x: float,
    half_z: float,
) -> Callable[[float, float], float]:
    nx = len(values[0])
    nz = len(values)

    def evaluate(x: float, z: float) -> float:
        if abs(x) > half_x + dx or abs(z) > half_z + dz:
            return 0.0
        fx = clamp((x + half_x) / max(dx, SMALL), 0.0, nx - 1.0)
        fz = clamp((z + half_z) / max(dz, SMALL), 0.0, nz - 1.0)
        return max(bilinear_sample(values, fx, fz), 0.0)

    return evaluate


def beamshape_profile(
    values: list[list[float]],
    half_length: float,
) -> Callable[[float, float], float]:
    nx = len(values[0])
    nz = len(values)

    def evaluate(x: float, z: float) -> float:
        if abs(x) > half_length or abs(z) > half_length:
            return 0.0
        fx = (x + half_length) / max(2.0 * half_length, SMALL) * (nx - 1)
        fz = (z + half_length) / max(2.0 * half_length, SMALL) * (nz - 1)
        return max(bilinear_sample(values, fx, fz), 0.0)

    return evaluate


def plane_bias_modifier(
    amplitude: float,
    offset: float,
    reference_scale: float,
    direction: tuple[float, float],
    origin: tuple[float, float] = (0.0, 0.0),
    min_factor: float = 0.0,
) -> Callable[[float, float], float]:
    dx, dz = direction
    length = math.sqrt(dx * dx + dz * dz)
    if length <= SMALL:
        dx, dz = 1.0, 0.0
    else:
        dx /= length
        dz /= length

    def evaluate(x: float, z: float) -> float:
        coord = ((x - origin[0]) * dx + (z - origin[1]) * dz) / max(reference_scale, SMALL)
        return max(1.0 + offset + amplitude * coord, min_factor)

    return evaluate


def drum_noise_modifier(
    modes: list[dict[str, float | int]],
    radius: float,
    min_factor: float = 0.0,
    time_value: float = 0.0,
) -> Callable[[float, float], float]:
    prepared: list[dict[str, float]] = []
    for mode in modes:
        azimuthal_order = int(mode["azimuthalOrder"])
        radial_order = int(mode["radialOrder"])
        alpha = estimate_bessel_zero(azimuthal_order, radial_order)
        n_samples = 512
        max_abs_value = SMALL
        mean_value = 0.0
        for sample in range(n_samples + 1):
            xi = sample / n_samples
            raw = bessel_j(azimuthal_order, alpha * xi)
            max_abs_value = max(max_abs_value, abs(raw))
            if azimuthal_order == 0:
                weight = 0.5 if sample in (0, n_samples) else 1.0
                mean_value += weight * raw * xi

        mean_value = (
            2.0 * mean_value / n_samples / max_abs_value if azimuthal_order == 0 else 0.0
        )
        base_amplitude = float(mode["amplitude"])
        phase = math.radians(float(mode.get("phase", 0.0)))
        frequency = float(mode.get("frequency", 0.0))
        temporal_phase = TAU * frequency * time_value
        prepared.append(
            {
                "azimuthal_order": azimuthal_order,
                "alpha": alpha,
                "normalisation": max_abs_value,
                "mean_value": mean_value,
                "current_phase": phase + temporal_phase,
                "current_amplitude": base_amplitude * math.cos(temporal_phase),
            }
        )

    def evaluate(x: float, z: float) -> float:
        r = math.sqrt(x * x + z * z)
        if r > radius:
            return 1.0
        theta = math.atan2(z, x) if r > SMALL else 0.0
        modulation = 0.0
        for mode in prepared:
            radial_value = bessel_j(mode["azimuthal_order"], mode["alpha"] * r / max(radius, SMALL))
            radial_value /= mode["normalisation"]
            angular_value = (
                1.0
                if mode["azimuthal_order"] == 0
                else math.cos(mode["azimuthal_order"] * theta + mode["current_phase"])
            )
            value = radial_value * angular_value
            if mode["azimuthal_order"] == 0:
                value -= mode["mean_value"]
            modulation += mode["current_amplitude"] * value
        return max(1.0 + modulation, min_factor)

    return evaluate


def cartesian_mode_noise_modifier(
    modes: list[dict[str, float | int]],
    half_width: float,
    half_height: float,
    min_factor: float = 0.0,
    time_value: float = 0.0,
) -> Callable[[float, float], float]:
    prepared: list[dict[str, float]] = []
    for mode in modes:
        phase = math.radians(float(mode.get("phase", 0.0)))
        frequency = float(mode.get("frequency", 0.0))
        temporal_phase = TAU * frequency * time_value
        prepared.append(
            {
                "x_order": int(mode["xOrder"]),
                "z_order": int(mode["zOrder"]),
                "current_phase": phase + temporal_phase,
                "current_amplitude": float(mode["amplitude"]) * math.cos(temporal_phase),
            }
        )

    def evaluate(x: float, z: float) -> float:
        if abs(x) > half_width + SMALL or abs(z) > half_height + SMALL:
            return 1.0
        xi = clamp((x + half_width) / max(2.0 * half_width, SMALL), 0.0, 1.0)
        zi = clamp((z + half_height) / max(2.0 * half_height, SMALL), 0.0, 1.0)
        modulation = 0.0
        for mode in prepared:
            x_basis = (
                1.0
                if mode["x_order"] == 0
                else math.cos(TAU * mode["x_order"] * xi + mode["current_phase"])
            )
            z_basis = (
                1.0
                if mode["z_order"] == 0
                else math.cos(TAU * mode["z_order"] * zi + mode["current_phase"])
            )
            modulation += mode["current_amplitude"] * x_basis * z_basis
        return max(1.0 + modulation, min_factor)

    return evaluate


def generate_grid_mask(
    nx: int,
    nz: int,
    half_width: float,
    half_height: float,
    amplitude: float,
    correlation_length: float,
    smooth_passes: int,
    seed: int,
) -> list[list[float]]:
    rng = random.Random(seed)
    mask = [[rng.gauss(0.0, 1.0) for _ in range(nx)] for _ in range(nz)]

    mean = sum(sum(row) for row in mask) / (nx * nz)
    for iz in range(nz):
        for ix in range(nx):
            mask[iz][ix] -= mean

    passes = max(smooth_passes, 0)
    if correlation_length > SMALL:
        dx = 2.0 * half_width / max(nx - 1, 1)
        dz = 2.0 * half_height / max(nz - 1, 1)
        h = max(dx, dz)
        passes = max(passes, int(correlation_length / max(h, SMALL) + 0.5))

    for _ in range(passes):
        blurred = [[0.0 for _ in range(nx)] for _ in range(nz)]
        for iz in range(nz):
            for ix in range(nx):
                total = 0.0
                count = 0
                for j in range(max(iz - 1, 0), min(iz + 1, nz - 1) + 1):
                    for i in range(max(ix - 1, 0), min(ix + 1, nx - 1) + 1):
                        total += mask[j][i]
                        count += 1
                blurred[iz][ix] = total / count
        mask = blurred

    mean = sum(sum(row) for row in mask) / (nx * nz)
    max_value = SMALL
    for iz in range(nz):
        for ix in range(nx):
            mask[iz][ix] -= mean
            max_value = max(max_value, abs(mask[iz][ix]))

    for iz in range(nz):
        for ix in range(nx):
            mask[iz][ix] = amplitude * mask[iz][ix] / max_value

    return mask


def grid_noise_modifier(
    mask: list[list[float]],
    half_width: float,
    half_height: float,
    min_factor: float = 0.0,
) -> Callable[[float, float], float]:
    nz = len(mask)
    nx = len(mask[0])

    def evaluate(x: float, z: float) -> float:
        if abs(x) > half_width + SMALL or abs(z) > half_height + SMALL:
            return 1.0
        fx = (x + half_width) / max(2.0 * half_width, SMALL) * (nx - 1)
        fz = (z + half_height) / max(2.0 * half_height, SMALL) * (nz - 1)
        sampled = bilinear_sample(mask, fx, fz)
        return max(1.0 + sampled, min_factor)

    return evaluate


def apply_modifier(
    base: Callable[[float, float], float],
    modifier: Callable[[float, float], float],
) -> Callable[[float, float], float]:
    def evaluate(x: float, z: float) -> float:
        return max(base(x, z), 0.0) * max(modifier(x, z), 0.0)

    return evaluate


def make_signed_field(
    modifier: Callable[[float, float], float],
) -> Callable[[float, float], float]:
    def evaluate(x: float, z: float) -> float:
        return modifier(x, z) - 1.0

    return evaluate


def intensity_color(value: float, vmax: float) -> tuple[int, int, int]:
    # Approximation of matplotlib's "plasma" colormap.
    stops = [
        (0.0, (13, 8, 135)),
        (0.15, (84, 3, 160)),
        (0.30, (139, 10, 165)),
        (0.45, (186, 54, 122)),
        (0.60, (221, 94, 65)),
        (0.75, (244, 136, 35)),
        (0.90, (253, 188, 43)),
        (1.0, (240, 249, 33)),
    ]
    t = 0.0 if vmax <= SMALL else clamp(value / vmax, 0.0, 1.0)
    for index in range(len(stops) - 1):
        x0, c0 = stops[index]
        x1, c1 = stops[index + 1]
        if t <= x1:
            local = 0.0 if x1 <= x0 else (t - x0) / (x1 - x0)
            return (
                int(round(lerp(c0[0], c1[0], local))),
                int(round(lerp(c0[1], c1[1], local))),
                int(round(lerp(c0[2], c1[2], local))),
            )
    return stops[-1][1]


def signed_color(value: float, vmax: float) -> tuple[int, int, int]:
    if vmax <= SMALL:
        return intensity_color(0.0, 1.0)
    # Map [-vmax, vmax] onto [0, 1] using the same plasma ramp.
    return intensity_color(value + vmax, 2.0 * vmax)


def rgb_string(rgb: tuple[int, int, int]) -> str:
    return f"rgb({rgb[0]},{rgb[1]},{rgb[2]})"


def make_canvas(width: int, height: int, background: tuple[int, int, int]) -> bytearray:
    pixels = bytearray(width * height * 3)
    for y in range(height):
        for x in range(width):
            index = (y * width + x) * 3
            pixels[index : index + 3] = bytes(background)
    return pixels


def set_pixel(
    pixels: bytearray,
    width: int,
    height: int,
    x: int,
    y: int,
    color: tuple[int, int, int],
) -> None:
    if x < 0 or x >= width or y < 0 or y >= height:
        return
    index = (y * width + x) * 3
    pixels[index : index + 3] = bytes(color)


def fill_rect(
    pixels: bytearray,
    width: int,
    height: int,
    x0: int,
    y0: int,
    rect_width: int,
    rect_height: int,
    color: tuple[int, int, int],
) -> None:
    x1 = min(x0 + rect_width, width)
    y1 = min(y0 + rect_height, height)
    for y in range(max(y0, 0), y1):
        for x in range(max(x0, 0), x1):
            set_pixel(pixels, width, height, x, y, color)


def stroke_rect(
    pixels: bytearray,
    width: int,
    height: int,
    x0: int,
    y0: int,
    rect_width: int,
    rect_height: int,
    color: tuple[int, int, int],
) -> None:
    for x in range(x0, x0 + rect_width):
        set_pixel(pixels, width, height, x, y0, color)
        set_pixel(pixels, width, height, x, y0 + rect_height - 1, color)
    for y in range(y0, y0 + rect_height):
        set_pixel(pixels, width, height, x0, y, color)
        set_pixel(pixels, width, height, x0 + rect_width - 1, y, color)


def draw_circle(
    pixels: bytearray,
    width: int,
    height: int,
    cx: float,
    cy: float,
    radius: float,
    color: tuple[int, int, int],
) -> None:
    r = max(int(math.ceil(radius)), 1)
    r2 = radius * radius
    for y in range(int(cy) - r, int(cy) + r + 1):
        for x in range(int(cx) - r, int(cx) + r + 1):
            dx = (x + 0.5) - cx
            dy = (y + 0.5) - cy
            if dx * dx + dy * dy <= r2:
                set_pixel(pixels, width, height, x, y, color)


def draw_ellipse_outline(
    pixels: bytearray,
    width: int,
    height: int,
    cx: float,
    cy: float,
    rx: float,
    ry: float,
    color: tuple[int, int, int],
) -> None:
    steps = 720
    for i in range(steps):
        angle = TAU * i / steps
        x = int(round(cx + rx * math.cos(angle)))
        y = int(round(cy + ry * math.sin(angle)))
        set_pixel(pixels, width, height, x, y, color)


def write_bmp(path: Path, width: int, height: int, pixels: bytearray) -> None:
    row_stride = width * 3
    padding = (4 - (row_stride % 4)) % 4
    image_size = (row_stride + padding) * height
    file_size = 14 + 40 + image_size

    with path.open("wb") as handle:
        handle.write(b"BM")
        handle.write(file_size.to_bytes(4, "little"))
        handle.write((0).to_bytes(2, "little"))
        handle.write((0).to_bytes(2, "little"))
        handle.write((54).to_bytes(4, "little"))

        handle.write((40).to_bytes(4, "little"))
        handle.write(width.to_bytes(4, "little", signed=True))
        handle.write(height.to_bytes(4, "little", signed=True))
        handle.write((1).to_bytes(2, "little"))
        handle.write((24).to_bytes(2, "little"))
        handle.write((0).to_bytes(4, "little"))
        handle.write(image_size.to_bytes(4, "little"))
        handle.write((2835).to_bytes(4, "little", signed=True))
        handle.write((2835).to_bytes(4, "little", signed=True))
        handle.write((0).to_bytes(4, "little"))
        handle.write((0).to_bytes(4, "little"))

        pad = b"\x00" * padding
        for row in range(height - 1, -1, -1):
            start = row * row_stride
            row_data = pixels[start : start + row_stride]
            bgr = bytearray()
            for x in range(width):
                idx = x * 3
                r, g, b = row_data[idx : idx + 3]
                bgr.extend((b, g, r))
            handle.write(bgr)
            handle.write(pad)


def scatter_image(
    points: list[tuple[float, float, float]],
    half_width_x: float,
    half_width_z: float,
    path: Path,
    circle_radius: float | None = None,
) -> None:
    width = 320
    height = 320
    margin = 14
    plot_width = width - 2 * margin
    plot_height = height - 2 * margin
    pixels = make_canvas(width, height, (249, 247, 242))
    accent = intensity_color(0.88, 1.0)
    outline = intensity_color(0.55, 1.0)

    def px(x: float) -> float:
        return margin + (x + half_width_x) / max(2.0 * half_width_x, SMALL) * plot_width

    def pz(z: float) -> float:
        return margin + (half_width_z - z) / max(2.0 * half_width_z, SMALL) * plot_height

    fill_rect(pixels, width, height, margin, margin, plot_width, plot_height, (255, 255, 255))
    stroke_rect(pixels, width, height, margin, margin, plot_width, plot_height, (150, 150, 150))

    if circle_radius is not None:
        rx = circle_radius / max(half_width_x, SMALL) * (plot_width / 2.0)
        rz = circle_radius / max(half_width_z, SMALL) * (plot_height / 2.0)
        draw_ellipse_outline(
            pixels,
            width,
            height,
            px(0.0),
            pz(0.0),
            rx,
            rz,
            outline,
        )

    for x, z, weight in points:
        radius = 1.2 + 2.8 * weight
        draw_circle(pixels, width, height, px(x), pz(z), radius, accent)
    write_bmp(path, width, height, pixels)


def heatmap_image(
    values: list[list[float]],
    half_width_x: float,
    half_width_z: float,
    path: Path,
    signed: bool = False,
    circle_radius: float | None = None,
) -> None:
    width = 320
    height = 320
    margin = 14
    plot_width = width - 2 * margin
    plot_height = height - 2 * margin
    nz = len(values)
    nx = len(values[0])
    cell_width = plot_width / nx
    cell_height = plot_height / nz
    vmax = max_abs(values) if signed else grid_max(values)
    pixels = make_canvas(width, height, (249, 247, 242))
    fill_rect(pixels, width, height, margin, margin, plot_width, plot_height, (255, 255, 255))
    stroke_rect(pixels, width, height, margin, margin, plot_width, plot_height, (150, 150, 150))

    for iz, row in enumerate(values):
        y = margin + iz * cell_height
        for ix, value in enumerate(row):
            x = margin + ix * cell_width
            color = signed_color(value, vmax) if signed else intensity_color(value, vmax)
            fill_rect(
                pixels,
                width,
                height,
                int(round(x)),
                int(round(y)),
                max(int(math.ceil(cell_width + 0.3)), 1),
                max(int(math.ceil(cell_height + 0.3)), 1),
                color,
            )

    if circle_radius is not None:
        rx = circle_radius / max(half_width_x, SMALL) * (plot_width / 2.0)
        rz = circle_radius / max(half_width_z, SMALL) * (plot_height / 2.0)
        draw_ellipse_outline(
            pixels,
            width,
            height,
            margin + plot_width / 2.0,
            margin + plot_height / 2.0,
            rx,
            rz,
            (235, 235, 235),
        )
    write_bmp(path, width, height, pixels)


def heatmap_sequence_image(
    fields: list[list[list[float]]],
    path: Path,
    signed: bool = False,
) -> None:
    frame_count = len(fields)
    if frame_count < 1:
        raise ValueError("heatmap_sequence_image requires at least one field")

    frame_size = 180
    gap = 10
    margin = 14
    width = 2 * margin + frame_count * frame_size + (frame_count - 1) * gap
    height = 2 * margin + frame_size
    pixels = make_canvas(width, height, (249, 247, 242))

    if signed:
        vmax = max(max_abs(field) for field in fields)
    else:
        vmax = max(grid_max(field) for field in fields)

    for frame_index, values in enumerate(fields):
        x0 = margin + frame_index * (frame_size + gap)
        y0 = margin
        fill_rect(pixels, width, height, x0, y0, frame_size, frame_size, (255, 255, 255))
        stroke_rect(pixels, width, height, x0, y0, frame_size, frame_size, (150, 150, 150))

        nz = len(values)
        nx = len(values[0])
        cell_width = frame_size / nx
        cell_height = frame_size / nz

        for iz, row in enumerate(values):
            y = y0 + iz * cell_height
            for ix, value in enumerate(row):
                x = x0 + ix * cell_width
                color = signed_color(value, vmax) if signed else intensity_color(value, vmax)
                fill_rect(
                    pixels,
                    width,
                    height,
                    int(round(x)),
                    int(round(y)),
                    max(int(math.ceil(cell_width + 0.3)), 1),
                    max(int(math.ceil(cell_height + 0.3)), 1),
                    color,
                )

    write_bmp(path, width, height, pixels)


def radial_polar_samples(
    radius: float,
    n_radial: int,
    n_angular: int,
) -> list[tuple[float, float, float]]:
    boundaries = [radius * i / n_radial for i in range(n_radial + 1)]
    annulus_areas = [math.pi * (boundaries[i + 1] ** 2 - boundaries[i] ** 2) for i in range(n_radial)]
    radial_points = [0.0 for _ in range(n_radial)]
    for i in range(1, n_radial):
        inner = boundaries[i]
        outer = boundaries[i + 1]
        radial_points[i] = (2.0 / 3.0) * (outer**3 - inner**3) / max(outer * outer - inner * inner, SMALL)

    total_samples = 1 + (n_radial - 1) * n_angular if n_radial > 0 else 0
    points: list[tuple[float, float, float]] = []
    max_area = max(annulus_areas[0], *(area / n_angular for area in annulus_areas[1:]))
    for global_idx in range(total_samples):
        if global_idx == 0:
            area = annulus_areas[0]
            points.append((0.0, 0.0, area / max_area))
            continue
        shifted = global_idx - 1
        ring_count = n_radial - 1
        i_theta = shifted // ring_count
        i_r = 1 + shifted % ring_count
        theta = TAU * i_theta / n_angular
        x = radial_points[i_r] * math.cos(theta)
        z = radial_points[i_r] * math.sin(theta)
        area = annulus_areas[i_r] / n_angular
        points.append((x, z, area / max_area))
    return points


def cartesian_samples(
    half_width_x: float,
    half_width_z: float,
    n_x: int,
    n_z: int,
) -> list[tuple[float, float, float]]:
    dx = 2.0 * half_width_x / n_x
    dz = 2.0 * half_width_z / n_z
    points: list[tuple[float, float, float]] = []
    for iz in range(n_z):
        for ix in range(n_x):
            x = -half_width_x + (ix + 0.5) * dx
            z = -half_width_z + (iz + 0.5) * dz
            points.append((x, z, 1.0))
    return points


def gated_samples(
    samples: list[tuple[float, float, float]],
    area: float | Callable[[tuple[float, float, float]], float],
    evaluator: Callable[[float, float], float],
    min_fraction: float,
) -> list[tuple[float, float, float]]:
    ray_powers: list[float] = []
    for sample in samples:
        local_area = area(sample) if callable(area) else area
        ray_powers.append(max(evaluator(sample[0], sample[1]), 0.0) * local_area)

    peak = max(ray_powers, default=0.0)
    if peak <= SMALL:
        return []
    threshold = min_fraction * peak
    retained: list[tuple[float, float, float]] = []
    for sample, power in zip(samples, ray_powers):
        if power + SMALL < threshold:
            continue
        retained.append((sample[0], sample[1], power / peak))
    return retained


def create_gallery(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    panels_by_section: list[tuple[str, str, list[Panel]]] = []

    csv_path = ROOT / "laserbeamFoam" / "csv_example" / "constant" / "intensity_moon.csv"
    csv_values = read_flat_scalar_grid(csv_path, 161, 161)
    csv_half = 1.6e-3
    csv_dx = 2.0e-5
    csv_dz = 2.0e-5
    csv_profile_fn = csv_profile(csv_values, csv_dx, csv_dz, csv_half, csv_half)

    beamshape_path = ROOT / "laserbeamFoam" / "opa_cbc_example" / "beamShape.inp"
    beamshape_values = read_beamshape_matrix(beamshape_path) if beamshape_path.exists() else None
    beamshape_half = 7.099084025e-4
    beamshape_profile_fn = (
        beamshape_profile(beamshape_values, beamshape_half) if beamshape_values is not None else None
    )

    # Section 1: sampling strategies.
    section_panels: list[Panel] = []
    radial_points = radial_polar_samples(radius=1.0, n_radial=9, n_angular=24)
    radial_path = output_dir / "sampling_radial_polar_layout.bmp"
    scatter_image(
        radial_points,
        1.1,
        1.1,
        radial_path,
        circle_radius=1.0,
    )
    section_panels.append(
        Panel(
            "Radial-polar sample layout",
            "Matches the ring-sector seeding used by `sampling.type radialPolar`.",
            radial_path.name,
        )
    )

    radial_profile = super_gaussian_profile(radius_x=0.42, radius_z=0.42, order=6.0)
    boundaries = [1.0 * i / 9 for i in range(10)]
    annulus_areas = [math.pi * (boundaries[i + 1] ** 2 - boundaries[i] ** 2) for i in range(9)]

    def radial_area(sample: tuple[float, float, float]) -> float:
        r = math.sqrt(sample[0] * sample[0] + sample[1] * sample[1])
        if r <= SMALL:
            return annulus_areas[0]
        for index in range(1, len(boundaries) - 1):
            if boundaries[index] <= r <= boundaries[index + 1] + SMALL:
                return annulus_areas[index] / 24.0
        return annulus_areas[-1] / 24.0

    radial_gated = gated_samples(radial_points, radial_area, radial_profile, min_fraction=0.08)
    radial_gated_path = output_dir / "sampling_radial_polar_gated_rays.bmp"
    scatter_image(
        radial_gated,
        1.1,
        1.1,
        radial_gated_path,
        circle_radius=1.0,
    )
    section_panels.append(
        Panel(
            "Radial-polar gated rays",
            "Power-based gating removes the dim outer sectors while keeping the on-axis sample.",
            radial_gated_path.name,
        )
    )

    cart_points = cartesian_samples(1.0, 0.7, 21, 15)
    cart_path = output_dir / "sampling_cartesian_lattice.bmp"
    scatter_image(
        cart_points,
        1.1,
        0.8,
        cart_path,
    )
    section_panels.append(
        Panel(
            "Cartesian sample lattice",
            "The sampling block supplies `halfWidth{X,Z}` and `nCartesian{X,Z}` directly.",
            cart_path.name,
        )
    )

    csv_points = cartesian_samples(csv_half, csv_half, 41, 41)
    csv_gated = gated_samples(csv_points, (2.0 * csv_half / 41) ** 2, csv_profile_fn, min_fraction=0.10)
    csv_gated_path = output_dir / "sampling_cartesian_csv_gated_rays.bmp"
    scatter_image(
        csv_gated,
        csv_half * 1.05,
        csv_half * 1.05,
        csv_gated_path,
    )
    section_panels.append(
        Panel(
            "Cartesian gated rays from CSV beam",
            "Cartesian sampling is the practical fit for imported spatial maps such as CSV intensity fields.",
            csv_gated_path.name,
        )
    )

    panels_by_section.append(
        (
            "Sampling strategies",
            "The seeding logic is now profile-local. Each profile carries its own `sampling` block.",
            section_panels,
        )
    )

    # Section 2: profiles and modifiers.
    section_panels = []
    profile_specs: list[tuple[str, str, Callable[[float, float], float], float, float, float | None]] = [
        (
            "Circular Gaussian",
            "The `gaussian` alias is the super-Gaussian implementation with order = 2.",
            super_gaussian_profile(
                radius_x=0.34,
                radius_z=0.34,
                order=2.0,
            ),
            0.95,
            0.95,
            None,
        ),
        (
            "Rotated elliptic super-Gaussian",
            "Elliptic support with rotation and center shift.",
            super_gaussian_profile(
                radius_x=0.42,
                radius_z=0.28,
                order=6.0,
                rotation_deg=28.0,
                center_shift=(0.08, -0.05),
            ),
            0.95,
            0.95,
            None,
        ),
        (
            "Soft annular ring",
            "Soft inner and outer transitions from the annular profile.",
            annular_profile(0.28, 0.56, 0.045, rotation_deg=15.0, peak_at_inner=False),
            0.9,
            0.9,
            None,
        ),
        (
            "First-order Bessel beam",
            "Ring structure from `bessel` with explicit first-zero radius.",
            bessel_profile(order=1, first_zero_radius=0.42, support_radius=0.86),
            0.92,
            0.92,
            None,
        ),
        (
            "CSV crescent import",
            "Imported intensity map from `csv_example/intensity_moon.csv`.",
            csv_profile_fn,
            csv_half,
            csv_half,
            None,
        ),
    ]

    if beamshape_profile_fn is not None:
        profile_specs.append(
            (
                "Beam-shape import map",
                "Imported primitive beam-shape data from `opa_cbc_example/beamShape.inp`.",
                beamshape_profile_fn,
                beamshape_half,
                beamshape_half,
                None,
            )
        )

    plane_bias = plane_bias_modifier(
        amplitude=0.65,
        offset=0.0,
        reference_scale=0.65,
        direction=(1.0, 0.35),
        origin=(-0.12, 0.0),
        min_factor=0.15,
    )
    cart_modes = cartesian_mode_noise_modifier(
        modes=[
            {"xOrder": 2, "zOrder": 0, "amplitude": 0.20, "phase": 10.0, "frequency": 0.8},
            {"xOrder": 1, "zOrder": 3, "amplitude": 0.16, "phase": 75.0, "frequency": 1.6},
        ],
        half_width=0.85,
        half_height=0.85,
        min_factor=0.08,
        time_value=0.23,
    )
    drum_modes = drum_noise_modifier(
        modes=[
            {"azimuthalOrder": 0, "radialOrder": 1, "amplitude": 0.18, "frequency": 0.6, "phase": 0.0},
            {"azimuthalOrder": 3, "radialOrder": 1, "amplitude": 0.14, "frequency": 1.1, "phase": 35.0},
        ],
        radius=0.58,
        min_factor=0.05,
        time_value=0.31,
    )
    grid_mask = generate_grid_mask(
        nx=33,
        nz=33,
        half_width=0.85,
        half_height=0.85,
        amplitude=0.28,
        correlation_length=0.12,
        smooth_passes=0,
        seed=19,
    )
    grid_modifier = grid_noise_modifier(grid_mask, half_width=0.85, half_height=0.85, min_factor=0.05)

    modified_specs = [
        (
            "Plane-biased super-Gaussian",
            "Modifiers multiply the base profile. Plane bias applies a directed linear tilt.",
            apply_modifier(
                super_gaussian_profile(radius_x=0.42, radius_z=0.42, order=4.0),
                plane_bias,
            ),
            0.95,
            0.95,
            None,
        ),
        (
            "Drum-noise annular ring",
            "Drum noise redistributes ring intensity in radial and azimuthal modes.",
            apply_modifier(
                annular_profile(0.30, 0.60, 0.04),
                drum_modes,
            ),
            0.92,
            0.92,
            None,
        ),
        (
            "Cartesian-mode super-Gaussian",
            "Cartesian modal noise is useful when the support is naturally rectangular.",
            apply_modifier(
                super_gaussian_profile(radius_x=0.46, radius_z=0.34, order=5.0, rotation_deg=-18.0),
                cart_modes,
            ),
            0.95,
            0.95,
            None,
        ),
        (
            "Grid-noise CSV beam",
            "Grid noise applies a smoothly correlated random mask on top of an imported map.",
            apply_modifier(csv_profile_fn, grid_modifier),
            csv_half,
            csv_half,
            None,
        ),
    ]

    for index, spec in enumerate(profile_specs + modified_specs, start=1):
        title, caption, func, half_x, half_z, circle_radius = spec
        values, _, _ = sample_grid(func, half_x, half_z)
        filename = f"profile_panel_{index:02d}.bmp"
        heatmap_image(values, half_x, half_z, output_dir / filename, circle_radius=circle_radius)
        section_panels.append(Panel(title, caption, filename))

    panels_by_section.append(
        (
            "Profiles and modifiers",
            "These panels map directly to the base profile types and multiplicative modifier stack in `HeatSourceProfiles.C`.",
            section_panels,
        )
    )

    # Section 3: weighted combinations.
    section_panels = []
    core_profile = super_gaussian_profile(radius_x=0.42, radius_z=0.42, order=6.0)
    ring_profile = annular_profile(0.45, 0.72, 0.04)
    base_half = 0.95
    core_values, dx, dz = sample_grid(core_profile, base_half, base_half)
    ring_values, _, _ = sample_grid(ring_profile, base_half, base_half)

    sweeps = [
        ("Core-ring mix 100/0", [1.0, 0.0]),
        ("Core-ring mix 75/25", [0.75, 0.25]),
        ("Core-ring mix 50/50", [0.50, 0.50]),
        ("Core-ring mix 25/75", [0.25, 0.75]),
        ("Core-ring mix 0/100", [0.0, 1.0]),
    ]
    for index, (title, weights) in enumerate(sweeps, start=1):
        combined = combine_weighted_profiles([core_values, ring_values], weights, dx, dz)
        filename = f"combination_core_ring_{index:02d}.bmp"
        caption = "Weights act on per-profile power allocation after each component is normalized."
        heatmap_image(combined, base_half, base_half, output_dir / filename)
        section_panels.append(Panel(title, caption, filename))

    bessel_values, bessel_dx, bessel_dz = sample_grid(
        bessel_profile(order=1, first_zero_radius=0.48, support_radius=0.95),
        base_half,
        base_half,
    )
    envelope_values, _, _ = sample_grid(
        super_gaussian_profile(radius_x=0.36, radius_z=0.36, order=4.0),
        base_half,
        base_half,
    )
    hybrid = combine_weighted_profiles([bessel_values, envelope_values], [0.70, 0.30], bessel_dx, bessel_dz)
    hybrid_name = "combination_bessel_supergaussian_hybrid.bmp"
    hybrid_caption = "A mixed beam using the same additive weighting mechanics as the multi-profile path."
    heatmap_image(hybrid, base_half, base_half, output_dir / hybrid_name)
    section_panels.append(Panel("Bessel-superGaussian hybrid 70/30", hybrid_caption, hybrid_name))

    panels_by_section.append(
        (
            "Weighted combinations",
            "Combinations are additive. The physically relevant part is that the weights distribute power between individually normalized components.",
            section_panels,
        )
    )

    # Section 4: noise functions.
    section_panels = []
    noise_half = 0.9
    bare_noise_specs = [
        (
            "Bare drum-mode field",
            "Signed modulation field: `drumNoise.evaluate(...) - 1`.",
            make_signed_field(
                drum_noise_modifier(
                    modes=[
                        {"azimuthalOrder": 0, "radialOrder": 1, "amplitude": 0.18, "frequency": 0.9, "phase": 0.0},
                        {"azimuthalOrder": 2, "radialOrder": 1, "amplitude": 0.15, "frequency": 1.4, "phase": 30.0},
                    ],
                    radius=0.62,
                    min_factor=0.05,
                    time_value=0.27,
                )
            ),
            noise_half,
            noise_half,
        ),
        (
            "Bare Cartesian modal field",
            "Signed modulation field from orthogonal cosine modes.",
            make_signed_field(
                cartesian_mode_noise_modifier(
                    modes=[
                        {"xOrder": 1, "zOrder": 2, "amplitude": 0.18, "phase": 25.0, "frequency": 0.7},
                        {"xOrder": 3, "zOrder": 1, "amplitude": 0.14, "phase": 110.0, "frequency": 1.3},
                    ],
                    half_width=noise_half,
                    half_height=noise_half,
                    min_factor=0.05,
                    time_value=0.18,
                )
            ),
            noise_half,
            noise_half,
        ),
        (
            "Bare correlated grid field",
            "Signed correlated random mask sampled bilinearly from a smoothed grid.",
            make_signed_field(
                grid_noise_modifier(
                    generate_grid_mask(
                        nx=33,
                        nz=33,
                        half_width=noise_half,
                        half_height=noise_half,
                        amplitude=0.32,
                        correlation_length=0.16,
                        smooth_passes=0,
                        seed=7,
                    ),
                    half_width=noise_half,
                    half_height=noise_half,
                    min_factor=0.05,
                )
            ),
            noise_half,
            noise_half,
        ),
    ]

    for index, (title, caption, func, half_x, half_z) in enumerate(bare_noise_specs, start=1):
        values, _, _ = sample_grid(func, half_x, half_z)
        filename = f"noise_bare_{index:02d}.bmp"
        heatmap_image(values, half_x, half_z, output_dir / filename, signed=True)
        section_panels.append(Panel(title, caption, filename))

    applied_specs = [
        (
            "Drum-mode super-Gaussian",
            apply_modifier(
                super_gaussian_profile(radius_x=0.36, radius_z=0.36, order=4.0),
                drum_noise_modifier(
                    modes=[
                        {"azimuthalOrder": 0, "radialOrder": 1, "amplitude": 0.16, "frequency": 0.9, "phase": 0.0},
                        {"azimuthalOrder": 3, "radialOrder": 1, "amplitude": 0.14, "frequency": 1.8, "phase": 55.0},
                    ],
                    radius=0.48,
                    min_factor=0.05,
                    time_value=0.22,
                ),
            ),
        ),
        (
            "Cartesian-mode super-Gaussian beam",
            apply_modifier(
                super_gaussian_profile(radius_x=0.46, radius_z=0.30, order=6.0),
                cartesian_mode_noise_modifier(
                    modes=[
                        {"xOrder": 2, "zOrder": 0, "amplitude": 0.18, "phase": 20.0, "frequency": 0.5},
                        {"xOrder": 1, "zOrder": 2, "amplitude": 0.13, "phase": 90.0, "frequency": 1.0},
                    ],
                    half_width=0.85,
                    half_height=0.85,
                    min_factor=0.05,
                    time_value=0.25,
                ),
            ),
        ),
        (
            "Correlated grid noise on CSV beam",
            apply_modifier(
                csv_profile_fn,
                grid_noise_modifier(
                    generate_grid_mask(
                        nx=33,
                        nz=33,
                        half_width=csv_half,
                        half_height=csv_half,
                        amplitude=0.30,
                        correlation_length=2.4e-4,
                        smooth_passes=0,
                        seed=11,
                    ),
                    half_width=csv_half,
                    half_height=csv_half,
                    min_factor=0.05,
                ),
            ),
        ),
    ]

    for index, (title, func) in enumerate(applied_specs, start=1):
        half_x = csv_half if "CSV" in title else noise_half
        half_z = csv_half if "CSV" in title else noise_half
        values, _, _ = sample_grid(func, half_x, half_z)
        filename = f"noise_applied_{index:02d}.bmp"
        caption = "Same modifier field as above, now multiplied by a base profile."
        heatmap_image(values, half_x, half_z, output_dir / filename)
        section_panels.append(Panel(title, caption, filename))

    time_sequence_specs = [
        (
            "Continuous drum-mode evolution",
            "Left to right: t = 0.00, 0.20, 0.40, 0.60. Same plasma scale across all frames.",
            [
                sample_grid(
                    make_signed_field(
                        drum_noise_modifier(
                            modes=[
                                {"azimuthalOrder": 0, "radialOrder": 1, "amplitude": 0.18, "frequency": 0.9, "phase": 0.0},
                                {"azimuthalOrder": 2, "radialOrder": 1, "amplitude": 0.15, "frequency": 1.4, "phase": 30.0},
                            ],
                            radius=0.62,
                            min_factor=0.05,
                            time_value=time_value,
                        )
                    ),
                    noise_half,
                    noise_half,
                )[0]
                for time_value in (0.00, 0.20, 0.40, 0.60)
            ],
            True,
        ),
        (
            "Continuous Cartesian-mode evolution",
            "Left to right: t = 0.00, 0.25, 0.50, 0.75. The phase motion comes from the modal frequencies.",
            [
                sample_grid(
                    make_signed_field(
                        cartesian_mode_noise_modifier(
                            modes=[
                                {"xOrder": 1, "zOrder": 2, "amplitude": 0.18, "phase": 25.0, "frequency": 0.7},
                                {"xOrder": 3, "zOrder": 1, "amplitude": 0.14, "phase": 110.0, "frequency": 1.3},
                            ],
                            half_width=noise_half,
                            half_height=noise_half,
                            min_factor=0.05,
                            time_value=time_value,
                        )
                    ),
                    noise_half,
                    noise_half,
                )[0]
                for time_value in (0.00, 0.25, 0.50, 0.75)
            ],
            True,
        ),
        (
            "Stepwise drum noise on super-Gaussian",
            "Left to right: four discrete updates representing successive stepwise random refreshes.",
            [
                sample_grid(
                    apply_modifier(
                        super_gaussian_profile(radius_x=0.36, radius_z=0.36, order=4.0),
                        drum_noise_modifier(
                            modes=[
                                {"azimuthalOrder": 0, "radialOrder": 1, "amplitude": amplitudes[0], "phase": phases[0]},
                                {"azimuthalOrder": 2, "radialOrder": 1, "amplitude": amplitudes[1], "phase": phases[1]},
                                {"azimuthalOrder": 4, "radialOrder": 1, "amplitude": amplitudes[2], "phase": phases[2]},
                            ],
                            radius=0.48,
                            min_factor=0.05,
                            time_value=0.0,
                        ),
                    ),
                    noise_half,
                    noise_half,
                )[0]
                for phases, amplitudes in (
                    ((0.0, 20.0, 45.0), (0.12, 0.10, 0.08)),
                    ((135.0, 65.0, 170.0), (0.14, 0.09, 0.07)),
                    ((250.0, 185.0, 25.0), (0.10, 0.12, 0.09)),
                    ((40.0, 300.0, 220.0), (0.13, 0.08, 0.10)),
                )
            ],
            False,
        ),
    ]

    for index, (title, caption, fields, signed) in enumerate(time_sequence_specs, start=1):
        filename = f"noise_time_sequence_{index:02d}.bmp"
        heatmap_sequence_image(fields, output_dir / filename, signed=signed)
        section_panels.append(Panel(title, caption, filename))

    panels_by_section.append(
        (
            "Noise fields",
            "The bare fields are shown as signed modulation. The applied fields show the same modulation after multiplication with a base profile. Time-sequence panels read left to right.",
            section_panels,
        )
    )

    html_lines = [
        "<!doctype html>",
        "<html lang='en'>",
        "<head>",
        "<meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>Beam Profile Capability Gallery</title>",
        "<style>",
        ":root { color-scheme: light; --bg: #f4efe6; --card: #fffdfa; --ink: #221c17; --muted: #64584d; --line: #d7cbb9; --accent: #935a2d; }",
        "body { margin: 0; font-family: Georgia, 'Times New Roman', serif; background: radial-gradient(circle at top, #fbf7ef 0%, var(--bg) 58%, #e9e0d4 100%); color: var(--ink); }",
        "main { max-width: 1360px; margin: 0 auto; padding: 28px 24px 60px; }",
        "h1 { margin: 0 0 10px; font-size: 34px; }",
        "p.lead { max-width: 980px; color: var(--muted); line-height: 1.5; margin: 0 0 28px; }",
        "section { margin-top: 34px; padding-top: 18px; border-top: 1px solid var(--line); }",
        "h2 { margin: 0 0 6px; font-size: 22px; }",
        "p.section-copy { margin: 0 0 16px; color: var(--muted); line-height: 1.45; }",
        ".grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 16px; }",
        ".card { background: linear-gradient(180deg, #fffdfa 0%, #f6f0e5 100%); border: 1px solid var(--line); border-radius: 14px; padding: 14px; box-shadow: 0 10px 30px rgba(66, 45, 17, 0.06); }",
        ".card h3 { margin: 0 0 8px; font-size: 17px; }",
        ".card p { margin: 10px 0 0; line-height: 1.45; color: var(--muted); font-size: 14px; }",
        "img { display: block; width: 100%; height: auto; border-radius: 8px; background: #fff; }",
        "code { font-family: 'SFMono-Regular', Consolas, monospace; font-size: 0.95em; color: var(--accent); }",
        "</style>",
        "</head>",
        "<body>",
        "<main>",
        "<h1>Beam Profile Capability Gallery</h1>",
        "<p class='lead'>This gallery mirrors the beam profile machinery implemented in <code>src/laserHeatSource/HeatSourceProfiles.C</code>. It covers profile-local sampling, the base profile families, modifier stacking, additive weighted combinations, and both bare and applied noise fields. Imported <code>csv</code> and <code>beamShape</code> examples are loaded from the tutorial assets already in this repository.</p>",
    ]

    for title, copy, panels in panels_by_section:
        html_lines.append("<section>")
        html_lines.append(f"<h2>{html.escape(title)}</h2>")
        html_lines.append(f"<p class='section-copy'>{html.escape(copy)}</p>")
        html_lines.append("<div class='grid'>")
        for panel in panels:
            html_lines.append("<article class='card'>")
            html_lines.append(f"<h3>{html.escape(panel.title)}</h3>")
            html_lines.append(
                f"<img src='{html.escape(panel.filename)}' alt='{html.escape(panel.title)}'>"
            )
            html_lines.append(f"<p>{html.escape(panel.caption)}</p>")
            html_lines.append("</article>")
        html_lines.append("</div>")
        html_lines.append("</section>")

    html_lines.extend(["</main>", "</body>", "</html>"])
    (output_dir / "index.html").write_text("\n".join(html_lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate an HTML gallery that visualizes beam profile capabilities."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Directory for the generated HTML and SVG files (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    create_gallery(args.output_dir.resolve())
    print(f"Wrote beam profile capability gallery to {args.output_dir.resolve() / 'index.html'}")


if __name__ == "__main__":
    main()
