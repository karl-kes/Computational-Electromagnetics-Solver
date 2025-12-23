# 3D Finite-Difference Time-Domain Electromagnetic Wave Simulator (WIP)
 
A high-performance Finite-Difference Time-Domain (FDTD) solver for Maxwell's equations in 3D, implemented in C++ with OpenMP parallelization and interactive Plotly visualization.

![C++](https://img.shields.io/badge/C++-17-blue?logo=c%2B%2B)
![OpenMP](https://img.shields.io/badge/OpenMP-Parallel-green)
![Python](https://img.shields.io/badge/Python-3.x-yellow?logo=python)
![License](https://img.shields.io/badge/License-MIT-brightgreen)

## Overview

This project implements the **Yee algorithm** for solving Maxwell's curl equations on a staggered grid. The FDTD method discretizes both space and time, updating electric and magnetic fields in a leapfrog scheme that naturally satisfies Faraday's and Ampère's laws.

### Maxwell's Equations (Vacuum)

The solver implements the two curl equations:

$$\frac{\partial \mathbf{B}}{\partial t} = -\nabla \times \mathbf{E}$$

$$\frac{\partial \mathbf{E}}{\partial t} = c^2 \nabla \times \mathbf{B}$$

### CFL Stability Condition

The time step is automatically calculated to satisfy the Courant-Friedrichs-Lewy condition:

$$\Delta t = \frac{0.99}{c \sqrt{\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2} + \frac{1}{\Delta z^2}}}$$

## Quick Start

### Prerequisites

- C++17 compiler (GCC recommended)
- OpenMP support
- Python 3.x with NumPy and Plotly

### Build & Run

```bash
# Compile with OpenMP (recommended)
g++ -std=c++17 main.cpp -o main.exe -fopenmp
```

### Without OpenMP

```bash
g++ -std=c++17 main.cpp -o main.exe
```

### Grid Constructor Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Nx, Ny, Nz` | Grid dimensions | 5 |
| `dx, dy, dz` | Spatial step sizes | 1.0 |
| `c` | Speed of light | 1.0 |

## Visualization

The included Python script generates an interactive 3D animation of the Ez field evolution:

```python
python render.py
```

## Project Structure

```
├── main.cpp          # Grid class and simulation driver
├── visualize.py      # Plotly animation generator
├── output/           # CSV slice exports
└── README.md
```

## Theory

### Yee Algorithm

The Yee scheme staggers E and B fields both spatially and temporally:

1. **B-field update** (half time step): Uses E-field curl at current positions
2. **E-field update** (full time step): Uses B-field curl at offset positions

This leapfrog integration preserves the divergence-free nature of the magnetic field and maintains second-order accuracy in both space and time.

### Boundary Conditions

Currently implements **absorbing boundaries** at grid edges (fields decay at boundaries). The update loops exclude boundary cells to prevent out-of-bounds access.

## Performance

*Results vary by hardware. The `collapse(2)` and `simd` pragmas optimize the nested loops for parallel execution.*

## Future Improvements

- [ ] Perfectly Matched Layer (PML) absorbing boundaries
- [ ] Material properties (permittivity, permeability, conductivity)
- [ ] Total-Field/Scattered-Field (TF/SF) source injection
- [ ] Near-to-far-field transformation
- [ ] GPU acceleration (CUDA/OpenCL)

## References

1. Yee, K. S. (1966). "Numerical solution of initial boundary value problems involving Maxwell's equations in isotropic media." *IEEE Transactions on Antennas and Propagation*.
2. Taflove, A., & Hagness, S. C. (2005). *Computational Electrodynamics: The Finite-Difference Time-Domain Method*. Artech House.