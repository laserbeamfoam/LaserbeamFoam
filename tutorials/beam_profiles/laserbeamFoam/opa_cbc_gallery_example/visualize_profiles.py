#!/usr/bin/env python3
"""Visualize CSV beam intensity profiles with quick ASCII previews."""

from pathlib import Path

def read_csv_profile(filename, nx, nz):
    """Read CSV profile into 2D list"""
    with open(filename, 'r') as f:
        values = [float(line.strip()) for line in f if line.strip()]

    expected = nx * nz
    if len(values) < expected:
        raise ValueError(
            f"{filename} contains {len(values)} values, expected {expected}"
        )
    values = values[:expected]
    
    # Reshape to 2D (row-major: z is outer loop)
    profile = []
    for zi in range(nz):
        row = values[zi*nx:(zi+1)*nx]
        profile.append(row)
    return profile

def ascii_render(profile, width=60):
    """Render 2D profile as ASCII art"""
    nz = len(profile)
    nx = len(profile[0])
    
    # Find max for normalization
    max_val = max(max(row) for row in profile)
    if max_val == 0:
        return "Empty profile"
    
    # ASCII characters from dark to bright
    chars = " .':-=+*#%@"
    
    # Downsample if needed
    skip_x = max(1, nx // width)
    skip_z = max(1, nz // width)
    
    result = []
    for zi in range(0, nz, skip_z):
        line = ""
        for xi in range(0, nx, skip_x):
            val = profile[zi][xi] / max_val
            char_idx = min(int(val * len(chars)), len(chars)-1)
            line += chars[char_idx]
        result.append(line)
    
    return '\n'.join(result)

def main():
    base_dir = Path(__file__).resolve().parent / "constant"
    
    shapes = [
        ('intensity_moon.csv', 161, 161, 'MOON (Smooth crescent)'),
        ('intensity_star.csv', 161, 161, 'STARBURST (Radial sharp spokes)'),
        ('intensity_ifsw.csv', 161, 161, 'IFSW (Logo-style text with reflection)'),
        ('intensity_pinwheel.csv', 161, 161, 'PINWHEEL (Swirling asymmetric lobes)'),
        ('intensity_clover.csv', 161, 161, 'CLOVER (Three-lobed beam)')
    ]
    
    print("=" * 70)
    print("CSV BEAM PROFILE VISUALIZATION")
    print("=" * 70)
    print()
    
    for filename, nx, nz, description in shapes:
        filepath = base_dir / filename
        if not filepath.exists():
            print(f"Missing {filename}")
            continue
        
        print(f"{description}")
        print(f"File: {filename} ({nx}x{nz} grid)")
        print("-" * 70)
        
        profile = read_csv_profile(filepath, nx, nz)
        ascii_art = ascii_render(profile, width=60)
        print(ascii_art)
        print()
        
        # Statistics
        flat = [val for row in profile for val in row]
        nonzero = sum(1 for v in flat if v > 0.01)
        print(f"Stats: max={max(flat):.3f}, nonzero points={nonzero}/{len(flat)}")
        print("=" * 70)
        print()

if __name__ == '__main__':
    main()
