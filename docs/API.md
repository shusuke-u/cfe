# API Reference

## Core Classes

### Field Class

The `cfe::field` class is the main class for managing the simulation field.

```cpp
class field {
public:
    // Constructor
    field(const sh::syntax& conf);
    
    // Material creation
    std::shared_ptr<material> make_material(const sh::syntax& conf);
    
    // Velocity management
    void clear_vel();
    
    // Simulation control
    void finish_setup();
    void compute_force();
    void next_step();
    
    // Physical quantities
    double mean_square_acc() const;
    double total_eng() const;
    double kinetic_eng() const;
    double angular_mom() const;
    double time() const;
    
    // Scaling
    void scale(double factor, size_t axis);
};
```

### Material Class

The `material` class represents the physical material in the simulation.

```cpp
class material {
public:
    // Subset management
    subset& subset(const std::string& name);
    
    // Velocity management
    void clear_vel();
    
    // Domain information
    domain domain_by_edge_ptcls() const;
    
    // Output
    void write(const std::string& filename, double time) const;
};
```

### Subset Class

The `subset` class represents a subset of particles in the material.

```cpp
class subset {
public:
    // Boundary conditions
    void fix(size_t axis);
    
    // Force calculation
    vector sum(double (particles::full_ptcl::*member)() const,
              std::function<double(const particles::full_ptcl&)> weight) const;
};
```

## Configuration

The simulation is configured using a Python file that defines various parameters:

```python
# Example configuration (tensile.py)
direction = "x"  # Tensile axis
tensile_length_ratio = 0.1  # 10% stretch
spring_to_sph_ratio = 1.0  # Spring to SPH ratio
msa_converge_bound = 1e-6  # Mean square acceleration convergence bound
min_step_bound = 1000  # Minimum number of steps
step_out = 100  # Output interval
```

## Common Parameters

- `direction`: The axis along which the force is applied ("x", "y", or "z")
- `tensile_length_ratio`: The ratio of stretch in the tensile direction
- `spring_to_sph_ratio`: The ratio between spring and SPH forces
- `msa_converge_bound`: Convergence criterion for mean square acceleration
- `min_step_bound`: Minimum number of simulation steps
- `step_out`: Interval for outputting results

## Error Handling

The simulation uses the following error handling mechanisms:

- `cfe::abort()`: Terminates the simulation with an error message
- `cfe::cout`: Thread-safe output stream for logging
- Exception handling for configuration errors

## Performance Considerations

- Use appropriate number of MPI processes based on your system
- Adjust `step_out` based on your needs (more frequent output = slower simulation)
- Consider using `clear_vel_interval` for long simulations to prevent numerical drift
