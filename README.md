# Linear Isotropic Elasticity Numerical Simulation Code with CFE method

This project provides numerical simulation code for analyzing the mechanical behavior of linear isotropic elastic bodies using SPH methods. The code is based on the CFE (Central Force Elasticity) method developed in Utsumi et al. (2025, in preparation).

## Features

- 3D linear isotropic elasticity analysis
- User friendly APIs
- Highly readable input and output format
- Parallel calculation with MPI and OpenMP
- Various sample simulation codes
  - Uniaxial tensile test
  - Linear waves propagation
  - Torsional test
  - Plate oscillation
  - Particle system relaxation
- Visualization tools for results

## Requirements

- C++20 or higher
- CMake 3.10 or higher
- Python 3.6 or higher
- MPI (Message Passing Interface)
- OpenMP (Open Multi-Processing)
- FDPS (Framework for Developing Particle Simulator) - <https://github.com/FDPS/FDPS.git>

## Installation

1. Clone and install FDPS:

```bash
git clone https://github.com/FDPS/FDPS.git
```

2. Clone this repository:

```bash
git clone [repository-url]
cd [repository-name]
```

1. Build the project (basically not required):

```bash
mkdir build
cmake -S . -B build
cmake --build build --target [your-target]
```

## Quick Start

1. Move to project directory (e.g., project/tensile)
2. Configure parameters in the settings file (tensile.py):
   - Material properties
   - Simulation parameters
   - Output settings
3. Run the simulation:

```bash
python3 run.py {n_procs}
```

where `{n_procs}` is the number of MPI processes to use
This scripts automatically build the required sources

4. Check the results in `runs/{timestamp}/`:

- Deformation data
- Stress-strain curves
- Visualization files

## Documentation

- [API Reference](docs/API.md): Detailed API documentation
- [Examples](docs/examples/): Example implementations and tutorials
  - [Tensile Test](docs/examples/tensile_test.md)
  - More examples coming soon...

## Sample Projects

1. **Tensile Test** (`project/tensile/`)
   - Uniaxial tensile test simulation
   - Calculates Young's modulus and Poisson's ratio
   - Configurable boundary conditions and material properties

2. **Linear Analysis** (`project/linear/`)
   - Linear wace propagation
   - Sound speed calculation of longitude and transverse
   - Sound speed ratio validation with respect to elastic dynamics

3. **Other Analysis Types**
   - Torsional lod test (`project/torsion/`)
   - Plate oscillation analysis (`project/plate/`)
   - Particle system relaxation (`project/relax/`)

## Output

The simulation generates the following outputs:

- Deformation data in binary format
- Stress-strain data in CSV format
- Visualization files for post-processing
- Log files with simulation parameters and performance metrics

## License

MIT License

This project uses FDPS (Framework for Developing Particle Simulator) which is also licensed under the MIT License.
See the LICENSE file for details.

## Contributing

Pull requests and issue reports are welcome. Please follow these steps:

1. Fork the repository
2. Create your feature branch
3. Commit your changes
4. Push to the branch
5. Create a Pull Request
