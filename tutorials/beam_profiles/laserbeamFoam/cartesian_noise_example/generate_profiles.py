#!/usr/bin/env python3

import math
from pathlib import Path

try:
    from PIL import Image
except ImportError as exc:
    raise SystemExit(
        "generate_profiles.py requires Pillow. Run it with "
        "`conda run -n ds python generate_profiles.py`."
    ) from exc


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "constant"

NX = 161
NZ = 161
SPAN_X = 3.2e-3
SPAN_Z = 3.2e-3
DX = SPAN_X / (NX - 1)
DZ = SPAN_Z / (NZ - 1)
HALF_X = 0.5 * DX * (NX - 1)
HALF_Z = 0.5 * DZ * (NZ - 1)
CONTENT_FRACTION = 0.88


def smooth_step(d, width):
    if width <= 0.0:
        return 1.0 if d <= 0.0 else 0.0
    return 0.5 * (1.0 - math.tanh(d / width))


def rect(x, z, cx, cz, hx, hz, blur):
    dx = abs(x - cx) - hx
    dz = abs(z - cz) - hz
    outside = math.hypot(max(dx, 0.0), max(dz, 0.0))
    inside = min(max(dx, dz), 0.0)
    return smooth_step(outside + inside, blur)


def circle(x, z, cx, cz, r, blur):
    return smooth_step(math.hypot(x - cx, z - cz) - r, blur)


def ring(x, z, cx, cz, r_inner, r_outer, blur):
    rr = math.hypot(x - cx, z - cz)
    return smooth_step(rr - r_inner, blur) * smooth_step(r_outer - rr, blur)


def segment(x, z, x1, z1, x2, z2, radius, blur):
    vx = x2 - x1
    vz = z2 - z1
    denom = vx * vx + vz * vz
    if denom <= 0.0:
        return circle(x, z, x1, z1, radius, blur)
    t = ((x - x1) * vx + (z - z1) * vz) / denom
    t = min(1.0, max(0.0, t))
    px = x1 + t * vx
    pz = z1 + t * vz
    return smooth_step(math.hypot(x - px, z - pz) - radius, blur)


def polygon(x, z, points, blur):
    inside = False
    min_dist = 1.0e9
    n = len(points)
    for i in range(n):
        x1, z1 = points[i]
        x2, z2 = points[(i + 1) % n]
        if ((z1 > z) != (z2 > z)):
            x_cross = (x2 - x1) * (z - z1) / (z2 - z1 + 1.0e-30) + x1
            if x < x_cross:
                inside = not inside

        vx = x2 - x1
        vz = z2 - z1
        denom = vx * vx + vz * vz
        t = 0.0 if denom <= 0.0 else ((x - x1) * vx + (z - z1) * vz) / denom
        t = min(1.0, max(0.0, t))
        px = x1 + t * vx
        pz = z1 + t * vz
        min_dist = min(min_dist, math.hypot(x - px, z - pz))

    return smooth_step(-min_dist if inside else min_dist, blur)


def add_max(accum, value, weight=1.0):
    return max(accum, weight * value)


def moon_profile(x, z):
    blur = 3.5e-5
    main = circle(x, z, 0.0, 0.0, 8.8e-4, blur)
    subtract = circle(x, z, 2.8e-4, 0.0, 8.0e-4, blur)
    glow = circle(x, z, -1.8e-4, 0.0, 4.2e-4, 1.2e-4)
    return max(main * (1.0 - 0.95 * subtract), 0.0) + 0.25 * glow


def star_profile(x, z):
    blur = 3.0e-5
    val = 0.0
    for angle_deg in (90, 162, 234, 306, 18):
        a = math.radians(angle_deg)
        x2 = 9.0e-4 * math.cos(a)
        z2 = 9.0e-4 * math.sin(a)
        val = add_max(val, segment(x, z, 0.0, 0.0, x2, z2, 1.1e-4, blur))
    val = add_max(val, circle(x, z, 0.0, 0.0, 2.4e-4, blur), 0.95)
    return val


def clover_profile(x, z):
    blur = 3.0e-5
    val = 0.0
    for angle_deg in (45, 135, 225, 315):
        a = math.radians(angle_deg)
        cx = 4.9e-4 * math.cos(a)
        cz = 4.9e-4 * math.sin(a)
        val = add_max(val, circle(x, z, cx, cz, 3.7e-4, blur))
    val = add_max(val, circle(x, z, 0.0, 0.0, 2.6e-4, blur), 0.8)
    return val


def pinwheel_profile(x, z):
    blur = 3.0e-5
    val = ring(x, z, 0.0, 0.0, 1.2e-4, 2.6e-4, blur)
    for angle_deg in (35, 125, 215, 305):
        a = math.radians(angle_deg)
        outer = (
            7.8e-4 * math.cos(a),
            7.8e-4 * math.sin(a)
        )
        inner1 = (
            2.8e-4 * math.cos(a + 0.62),
            2.8e-4 * math.sin(a + 0.62)
        )
        inner2 = (
            2.8e-4 * math.cos(a - 0.18),
            2.8e-4 * math.sin(a - 0.18)
        )
        poly = [(0.0, 0.0), inner1, outer, inner2]
        val = add_max(val, polygon(x, z, poly, blur))
    return val


PROFILES = {
    "intensity_moon.csv": moon_profile,
    "intensity_star.csv": star_profile,
    "intensity_pinwheel.csv": pinwheel_profile,
    "intensity_clover.csv": clover_profile,
}

IMAGE_PROFILES = {
    "intensity_ifsw.csv": OUT / "image.png",
}


def normalise(values):
    vmax = max(values)
    if vmax <= 0.0:
        return values
    return [v / vmax for v in values]


def blue_ink_strength(rgba):
    r, g, b, a = rgba
    if a <= 0:
        return 0.0

    white_distance = math.sqrt(
        (255 - r) ** 2 + (255 - g) ** 2 + (255 - b) ** 2
    ) / (math.sqrt(3.0) * 255.0)

    blue_excess = max(0.0, b - 0.5 * (r + g)) / 255.0
    return (a / 255.0) * white_distance * max(0.35, blue_excess)


def crop_content_bounds(image, threshold=0.03):
    width, height = image.size
    xs = []
    ys = []

    for y in range(height):
        for x in range(width):
            if blue_ink_strength(image.getpixel((x, y))) > threshold:
                xs.append(x)
                ys.append(y)

    if not xs:
        raise ValueError("No non-white content found in image profile.")

    return (min(xs), min(ys), max(xs) + 1, max(ys) + 1)


def write_image_profile(path, image_path):
    image = Image.open(image_path).convert("RGBA")
    bounds = crop_content_bounds(image)
    image = image.crop(bounds)

    resampling = getattr(Image, "Resampling", Image)
    max_width = max(1, int(round(CONTENT_FRACTION * NX)))
    max_height = max(1, int(round(CONTENT_FRACTION * NZ)))
    image.thumbnail((max_width, max_height), resampling.LANCZOS)

    canvas = Image.new("RGBA", (NX, NZ), (255, 255, 255, 0))
    offset_x = (NX - image.size[0]) // 2
    offset_y = (NZ - image.size[1]) // 2
    canvas.paste(image, (offset_x, offset_y), image)

    values = []
    for iz in range(NZ):
        for ix in range(NX):
            values.append(blue_ink_strength(canvas.getpixel((ix, iz))))

    values = normalise(values)

    with path.open("w", encoding="ascii") as f:
        for v in values:
            f.write(f"{v:.8f}\n")


def write_profile(path, func):
    values = []
    for iz in range(NZ):
        z = -HALF_Z + iz * DZ
        for ix in range(NX):
            x = -HALF_X + ix * DX
            values.append(max(func(x, z), 0.0))

    values = normalise(values)

    with path.open("w", encoding="ascii") as f:
        for v in values:
            f.write(f"{v:.8f}\n")


def main():
    for name, func in PROFILES.items():
        write_profile(OUT / name, func)

    for name, image_path in IMAGE_PROFILES.items():
        write_image_profile(OUT / name, image_path)


if __name__ == "__main__":
    main()
