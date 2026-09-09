# 3D heat equation

The second Kestrel reference application solves the 3D heat equation with a
forward Euler, seven-point finite-difference stencil and insulated (zero-flux
Neumann) boundaries. It retains the original CPU/OpenMP and CUDA implementations,
Gaussian and cosine initial conditions, and binary VTK output.

Imported from [karl-kes/3d-heat-solver](https://github.com/karl-kes/3d-heat-solver)
at commit `f10a9c3ecbfdf691507caf3e4b5410dccbbdb540`; see [LICENSE](LICENSE).
The original numerical kernels and regression tests are preserved. The example
uses Kestrel's build/backend selection, but still owns its original grid and
storage classes. Migration to `kestrel::grid`, fields, and XPU iteration is a
subsequent step once their layout and indexing contracts are implemented.

## Build and run

From the Kestrel repository root, for the CPU backend:

```bash
cmake -S . -B build-heat-test \
    -DCMAKE_CXX_COMPILER=g++-15 \
    -DCMAKE_BUILD_TYPE=Test \
    -DKESTREL_ENABLE_CUDA=OFF
cmake --build build-heat-test --parallel
ctest --test-dir build-heat-test -L heat --output-on-failure

cd build-heat-test/examples/heat
OMP_NUM_THREADS=4 ./kestrel-heat --nx 32 --ny 24 --nz 16 --steps 100
./kestrel-heat --help
```

Use a separate build directory with `-DKESTREL_ENABLE_CUDA=ON` for CUDA. The
parent project selects CUDA when its compiler is available, otherwise CPU.
Without a visible GPU, specify the intended target explicitly, for example
`-DCMAKE_CUDA_ARCHITECTURES=75` for a compute-capability 7.5 device, to compile
without native architecture detection. Execution still requires a supported GPU.

Heat defaults to `float`, independently of the developing framework precision
API. Set `-DKESTREL_HEAT_PRECISION=double` to build the example and its tests in
double precision. `KESTREL_BUILD_EXAMPLES=OFF` excludes both reference examples.

## Configuration and output

| Argument | Meaning | Default |
|----------|---------|---------|
| `--nx`, `--ny`, `--nz` | Point counts, including boundary/ghost layers; each at least 3 | 64 |
| `--steps` | Integration steps | 1000 |
| `--alpha` | Positive thermal diffusivity | 1.0 |
| `--dx`, `--dy`, `--dz` | Positive spatial spacing | 1.0 |
| `--ic` | `gaussian` or `cosine` | `gaussian` |
| `--output-interval` | VTK interval; zero disables output | 0 |

The solver derives its timestep from diffusivity and spacing using its explicit
stability limit. VTK files are written to `out/` relative to the working directory
and can be opened in ParaView. Run from the build directory to keep output there:

```bash
OMP_NUM_THREADS=4 ./kestrel-heat --ic cosine --steps 20 --output-interval 10
```

## Validation

Debug and Test builds register ten `heat-*` CTest cases: Gaussian diffusion and
conservation, the Neumann cosine solution, configuration, aligned allocation,
initialization, convergence order, a convergence probe, and three invalid CLI
inputs. CUDA numerical tests report a skip when no CUDA device is available.
CPU tests use four OpenMP threads to avoid oversubscribing small workloads.

The upstream paper and benchmarking scripts remain in the original repository;
their recorded performance numbers are not measurements of this integration.
