# CUDA Support for `laserHeatSource`

The `laserHeatSource` library includes optional GPU acceleration for ray
 tracing. CUDA support is **disabled by default** and is enabled automatically
 when the environment variable `CUDA_DIR` is set before building.

## 1. Setting `CUDA_DIR`

### If CUDA was installed via Ubuntu (`nvidia-cuda-toolkit`)

```bash
export CUDA_DIR=/usr
```

### If CUDA was installed via NVIDIA’s official installer

```bash
export CUDA_DIR=/usr/local/cuda
```

### How to determine which value to use

Run:

```bash
which nvcc
```

- If it prints `/usr/bin/nvcc` → use `CUDA_DIR=/usr`
- If it prints `/usr/local/cuda/bin/nvcc` → use `CUDA_DIR=/usr/local/cuda`

## 2. Building the library

```bash
./Allwmake
```

If `CUDA_DIR` is set correctly, the build system will:

- add CUDA include paths
- link against the CUDA runtime (`libcudart`)
- define `FOAM_USE_CUDA` so GPU code is compiled

If `CUDA_DIR` is **not** set, the library builds in CPU-only mode.

## 3. Runtime behaviour

If compiled with CUDA support, `laserHeatSource` will:

- automatically detect an available NVIDIA GPU
- fall back to CPU ray tracing if no GPU is present

CUDA support is optional and does not affect normal CPU-only builds.
