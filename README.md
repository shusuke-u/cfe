# Elastic Body Mechanics Numerical Simulation Code

This project provides numerical simulation code for analyzing the mechanical behavior of elastic bodies using particle-based methods.

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

## Usage

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
4. Check the results in `runs/{timestamp}/`:

- Deformation data
- Stress-strain curves
- Visualization files

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
[Add citation information here]
