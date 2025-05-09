# Elastic Body Mechanics Numerical Simulation Code

This project provides numerical simulation code for analyzing the mechanical behavior of elastic bodies using particle-based methods. The code is based on the CFE (Central Force Elasticity) method developed in Utsumi et al. (2025, in preparation).

## Features

- 3D elastic body deformation analysis
- Tensile test simulation
- Poisson's ratio calculation
- Young's modulus calculation
- Parallel computation support
- Visualization tools for results

## Requirements

- C++17 or higher
- CMake 3.10 or higher
- Python 3.6 or higher
- MPI (Message Passing Interface)
- FDPS (Framework for Developing Particle Simulator) - <https://github.com/FDPS/FDPS.git>

## Installation

1. Clone and install FDPS:

```bash
git clone https://github.com/FDPS/FDPS.git
cd FDPS
./install.sh
```

2. Clone this repository:

```bash
git clone [repository-url]
cd [repository-name]
```

3. Build the project:

```bash
mkdir build
cd build
cmake ..
make
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

## Available Projects

1. **Tensile Test** (`project/tensile/`)
   - Uniaxial tensile test simulation
   - Calculates Young's modulus and Poisson's ratio
   - Configurable boundary conditions and material properties

2. **Linear Analysis** (`project/linear/`)
   - Linear elastic analysis
   - Stress-strain relationship analysis
   - Material property validation

3. **Other Analysis Types**
   - Biaxial test (`project/biaxial/`)
   - Torsion test (`project/torsion/`)
   - Plate analysis (`project/plate/`)
   - Eigenvalue analysis (`project/eigen_fn/`)

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

## Citation

If you use this code in your research, please cite:

```
@article{utsumi2025cfe,
  title={Coupled Finite Element Method for Elastic Body Simulation},
  author={Utsumi, Shusuke and others},
  journal={in preparation},
  year={2025}
}
```
