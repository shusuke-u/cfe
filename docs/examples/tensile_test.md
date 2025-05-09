# Tensile Test Example

This example demonstrates how to perform a uniaxial tensile test simulation.

## Configuration

Create a configuration file `tensile.py`:

```python
# Material properties
direction = "x"  # Tensile axis
tensile_length_ratio = 0.1  # 10% stretch
spring_to_sph_ratio = 1.0  # Spring to SPH ratio

# Simulation parameters
msa_converge_bound = 1e-6  # Mean square acceleration convergence bound
min_step_bound = 1000  # Minimum number of steps
step_out = 100  # Output interval
clear_vel_interval = 1000  # Clear velocities every 1000 steps
clear_vel_lower = 100  # Start clearing velocities after 100 steps

# Convergence criteria
poisson_converge_bound = 1e-3  # Poisson's ratio convergence bound
```

## Implementation

```cpp
#include <cfe.h>

int main(int argc, char* argv[]) {
    // Initialize
    cfe::initialize(argc, argv);
    sh::syntax conf("tensile.py");
    
    // Setup
    auto field = std::make_shared<cfe::field>(conf);
    auto material = field->make_material(conf);
    field->clear_vel();
    field->finish_setup();
    
    // Set boundary conditions
    auto tensile_axis = sh::dir(conf["direction"]);
    auto& lower_side = material->subset("lower." + conf["direction"]);
    auto& upper_side = material->subset("upper." + conf["direction"]);
    
    lower_side.fix(tensile_axis);
    upper_side.fix(tensile_axis);
    
    // Scale material
    auto poisson_estimated = calculate_poisson_ratio(conf);
    for(auto i : cfe::vector::indices()) {
        field->scale(
            i == tensile_axis ? 1 + conf["tensile_length_ratio"]
                            : 1 - conf["tensile_length_ratio"] * poisson_estimated,
            i
        );
    }
    
    // Main loop
    for(auto step : std::views::iota(0UL)) {
        if(step % conf["step_out"] == 0) {
            // Calculate and output results
            output_results(field, material, step);
            
            // Check convergence
            if(is_converged(field, material, conf)) break;
        }
        
        // Clear velocities if needed
        if(should_clear_velocities(step, conf)) {
            material->clear_vel();
        }
        
        field->next_step();
    }
    
    cfe::finalize();
}
```

## Results

The simulation generates the following outputs:

1. **Physical Quantities**
   - Young's modulus
   - Poisson's ratio
   - Stress-strain relationship
   - Total energy
   - Kinetic energy
   - Angular momentum

2. **Output Files**
   - `result.yml`: YAML format results
   - `result.tsv`: Tab-separated values
   - `dump/{step}.tsv`: Particle data at each output step

## Analysis

To analyze the results:

1. Check the convergence of physical quantities
2. Plot stress-strain curves
3. Calculate material properties
4. Verify energy conservation

## Common Issues

1. **Non-convergence**
   - Increase `min_step_bound`
   - Adjust `msa_converge_bound`
   - Check boundary conditions

2. **Numerical Instability**
   - Decrease `tensile_length_ratio`
   - Adjust `clear_vel_interval`
   - Check material properties

3. **Performance Issues**
   - Adjust `step_out`
   - Optimize number of MPI processes
   - Check system resources
