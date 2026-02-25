# PSC3 / GOTPM --- Cosmological N-body Simulation Code

**PSC3** (Parallel Simulation Code 3), also known as **GOTPM** (Grid-Of-Oct-Tree Particle Mesh), is a cosmological N-body simulation code that computes gravitational dynamics of dark matter particles in an expanding universe using the **TreePM** method.

## Features

- **TreePM method**: PM (Particle-Mesh) for long-range forces via FFT-based Poisson solver + Barnes-Hut tree for short-range corrections
- **CPL dark energy** equation of state: w(z) = w0 + wa * z/(1+z)
- **2nd-order Lagrangian Perturbation Theory (2LPT)** initial conditions
- **Hybrid MPI+OpenMP** parallelization
- **Variable PM** (VarPM) domain decomposition for load balancing
- **Index-based position encoding** (XYZDBL+INDEX) for sub-cell displacements preserving double-precision accuracy
- **TSC** (Triangular-Shaped Cloud) mass assignment
- **FFTW3**-based distributed FFT for the Poisson solver
- **Poisson correction factor** for non-LCDM cosmologies
- Custom **stack-based memory pool** allocator

## Version History

| Version | Date | Notes |
|---------|------|-------|
| 1.4 | Mar 2024 | Newer CAMB for various DarkEnergyModel |
| 1.3 | Jul 2009 | User-friendly I/O, double-precision tree, no CPU/grid limitation |
| 1.2 | 2009 | Double-precision accuracy via XYZDBL index encoding |
| 1.1 | Apr 2009 | CAMB power spectrum integration |

## Prerequisites

- **Intel C/Fortran Compilers** (`mpiicx`, `mpiifx`) with MPI wrappers
- **OpenMP** support (`-qopenmp`)
- **FFTW3** library (single-precision, with MPI and OpenMP thread support)
- **SLURM** job scheduler (for cluster execution)
- A power spectrum file (CAMB output or ASCII format)

## Building

1. Edit `Rules.make` to set compiler paths and options:
   - `FFTWDIR` --- path to FFTW3 installation
   - `NMEG` --- memory pool size in megabytes per MPI rank (e.g., `8000L`)
   - `CC`, `FC`, `F90C` --- compiler wrappers

2. Build:
   ```bash
   make clean
   make this
   ```

This produces the executable `namu.exe`.

## Running a Simulation

### Directory Setup

```bash
mkdir run_directory
cd run_directory
cp /path/to/PSC3_CPL_mv_new/namu.exe .
cp params.dat .
cp camb.z=47_0.31_-1.2_-2.4.ascii .   # power spectrum file
```

### Parameter File (`params.dat`)

The simulation is configured via a plain-text parameter file using `define KEY = VALUE` format:

```
define INITIAL CONDITION   = 2          # 1=Zeldovich, 2=2LPT
define OmegaMatter0        = 0.31
define OmegaLambda0        = 0.69
define w0 of DE (CPL)      = -1.2
define wa of DE (CPL)      = -2.4
define Hubble              = 0.72
define nPS                 = 0.96
define Boxsize(Mpc/h)      = 1024
define Nx                  = 1024
define Ny                  = 1024
define Nz                  = 1024
define Astep               = 0.025
define Nstep               = 1881
define Iseed               = -56
define Powerflag           = 2
define Ascii Powerfile     = camb.z=47_0.31_-1.2_-2.4.ascii
```

### Job Submission (SLURM)

```bash
#!/bin/bash
#SBATCH --job-name=MV311224
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH --time=1000000:00
#SBATCH --output outfile.o%j.grammar
#SBATCH --error errfile.e%j.grammar
#SBATCH -p normal

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun ./namu.exe params.dat
```

> **Note**: `ntasks` must divide `Nz` evenly. The product `ntasks x cpus-per-task` should match available CPU cores.

### Restart

```bash
srun ./namu.exe params.00101    # restart from step 101
```

## Compile-Time Options

Set in `Rules.make` via `CDFLAGS` and `COMFLAGS`:

| Flag | Description |
|------|-------------|
| `-DNMEG=nL` | Memory pool size per MPI rank in MB |
| `-DWGROUPSIZE=n` | I/O group size |
| `-DINCLUDE_TREE_FORCE` | Enable tree force corrections |
| `-DINDEX` | Index-based position encoding (mandatory) |
| `-DXYZDBL` | Double-precision position recovery |
| `-DVarPM` | Variable PM domain decomposition (mandatory) |
| `-DPMSEEDFORCE` | PM seed force for initialization |
| `-DSAVESLICE` | Save x-z slice density files |
| `-D_LARGE_FILES` | Large file support (>4GB) |
| `-DMPI_SLOW` | Conservative MPI communication |
| `-DDEMODEL=n` | Dark energy Poisson correction (0=scale-independent, 1=scale-dependent) |

## Source Files

### Core

| File | Description |
|------|-------------|
| `main.c` | Program entry point; MPI/FFTW init; parameter reading |
| `namu.zeldovich.c` | Zeldovich IC; `pm2treetype()`/`tree2pmtype()` |
| `namu.2lpt.c` | 2LPT IC and main simulation time loop |
| `pmmain.c` | PM force: TSC, FFT, Poisson solver, FDA4 kick |
| `treemain.c` | Tree force correction driver |
| `kjhtpm.mod1.gpu.omp.c` | Tree correction: cell linked lists, tree walk (OpenMP) |
| `Treewalk.c` | Barnes-Hut tree traversal |
| `Memory2.c` | Stack-based memory pool allocator |

### PM Components

| File | Description |
|------|-------------|
| `tsc_omp2.c` | TSC mass assignment (OpenMP) |
| `fda4.c` | 4-point finite-difference acceleration |
| `p_solver.c` | Poisson solver in Fourier space |
| `fft.c` | FFTW3 wrapper functions |
| `VarPM.c` | Variable PM domain decomposition |

### Parallelization

| File | Description |
|------|-------------|
| `pmigrate.kisti.omp.c` | Particle migration between MPI ranks |
| `mydomaindecomposition.omp.c` | Dynamic domain rebalancing |
| `slicemigrate.c` | Slab-based mesh data migration |

### Cosmology

| File | Description |
|------|-------------|
| `FLRW.F` | Growth factors, H(z), Omega(z), CPL DE model |
| `pcorr.c` | Poisson correction factor interpolation |

### I/O

| File | Description |
|------|-------------|
| `kjhrw.c` | Particle data read/write |
| `dumpmain.c` | Synchronized particle data output |
| `savexzslice.c` | x-z slice density output |
| `header.c` | Parameter file I/O |

## Flow Control Files

- **`Suddenstop.flag`**: Controls graceful shutdown. Number followed by `+` (continue) or `-` (stop and save).
- **`WriteSync&WholeDen.flag`**: Redshift list for synchronized data output, power spectrum, and density dumps.

## Documentation

The full technical reference manual is available in LaTeX format under `docs/manual/`:
```bash
cd docs/manual
pdflatex main.tex
```

## Notes

- The tree opening angle `theta` is adaptively reduced at early steps (high redshift) for accuracy.
- FFTW3 single-precision libraries are required (`-lfftw3f_threads -lfftw3f_mpi -lfftw3f -lfftw3f_omp`).
- Use `-nofor-main` when linking mixed C/Fortran objects with Intel compilers.
