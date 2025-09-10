# OOP Coagulation Simulation Framework

This document describes the object-oriented refactoring of the original MATLAB coagulation simulation code.

## Overview

The original procedural code has been refactored into a modular, object-oriented architecture that provides better organization, maintainability, and reusability while preserving the exact mathematical formulations.

## Class Structure

### Core Configuration Classes

#### `SimulationConfig`
- **Purpose**: Manages all simulation parameters and configuration
- **Key Properties**: Physical constants, grid parameters, time settings, kernel options
- **Methods**: `derive()`, `validate()`
- **Usage**: Primary interface for setting simulation parameters

#### `DerivedGrid`
- **Purpose**: Precomputed constants and grids derived from configuration
- **Key Properties**: Volume bounds, fractal parameters, settling constants
- **Methods**: `getFractalRadii()`, `getConservedRadii()`, `getImageDiameters()`
- **Usage**: Provides efficient access to derived quantities

### Computational Classes

#### `KernelLibrary`
- **Purpose**: Implements all coagulation kernel functions
- **Methods**: `brownian()`, `curvilinearDS()`, `curvilinearShear()`, `fractalDS()`, `fractalShear()`, `rectilinearDS()`, `rectilinearShear()`
- **Usage**: Stateless service for kernel evaluations

#### `SettlingVelocityService`
- **Purpose**: Calculates particle settling velocities
- **Methods**: `velocity()`, `velocityWithConfig()`, `velocityForSections()`
- **Usage**: Utility service for settling calculations

#### `LinearProcessBuilder`
- **Purpose**: Builds linear operators for growth, sinking, and disaggregation
- **Methods**: `growthMatrix()`, `sinkingMatrix()`, `disaggregationMatrices()`, `linearMatrix()`
- **Usage**: Creates matrices for linear processes

#### `InitialSpectrumBuilder`
- **Purpose**: Generates initial particle size spectra
- **Methods**: `initialSpectrum()`, `customSpectrum()`, `exponentialSpectrum()`, `powerLawSpectrum()`
- **Usage**: Flexible initial condition generation

### Coagulation Engine Classes

#### `BetaAssembler`
- **Purpose**: Computes sectionally integrated coagulation kernel matrices
- **Methods**: `computeFor()`, `combineAndScale()`, `scaleBetas()`
- **Usage**: Core coagulation matrix computation

#### `BetaMatrices`
- **Purpose**: Container for coagulation kernel matrices
- **Properties**: `b1`, `b2`, `b3`, `b4`, `b5`, `b25`
- **Methods**: `validate()`, `displaySummary()`
- **Usage**: Holds and validates beta matrices

#### `CoagulationRHS`
- **Purpose**: ODE right-hand side evaluation and Jacobian computation
- **Methods**: `evaluate()`, `jacobian()`, `rateTerms()`, `validate()`
- **Usage**: Core numerical integration interface

### Integration and Analysis Classes

#### `ODESolver`
- **Purpose**: Adapter for MATLAB ODE solvers
- **Methods**: `solve()`, `setOptions()`, `setSolver()`
- **Usage**: Flexible time integration

#### `MassBalanceAnalyzer`
- **Purpose**: Analyzes mass balance and diagnostics
- **Methods**: `sectional()`, `total()`, `sectionalWithRates()`, `displayBalanceSummary()`
- **Usage**: Post-simulation analysis

#### `OutputGenerator`
- **Purpose**: Generates visualizations and exports data
- **Methods**: `spectraAndFluxes()`, `plotAll()`, `exportData()`
- **Usage**: Comprehensive output generation

### Main Controller

#### `CoagulationSimulation`
- **Purpose**: Orchestrates the entire simulation process
- **Methods**: `run()`, `generateOutputs()`, `computeDiagnostics()`, `exportResults()`
- **Usage**: Main interface for running simulations

## Usage Examples

### Basic Simulation
```matlab
% Create and run simulation with default parameters
sim = CoagulationSimulation();
result = sim.run();

% Generate outputs and plots
sim.generateOutputs(true);
```

### Custom Configuration
```matlab
% Create custom configuration
config = SimulationConfig(...
    'n_sections', 15, ...
    't_final', 10.0, ...
    'growth', 0.1, ...
    'gamma', 0.05
);

% Run simulation with custom config
sim = CoagulationSimulation(config);
result = sim.run();
```

### Custom Initial Conditions
```matlab
% Create exponential initial spectrum
v0 = InitialSpectrumBuilder.exponentialSpectrum(config, grid, 1000, 0.2);

% Run with custom initial conditions
result = sim.run('v0', v0);
```

### Custom Solver Options
```matlab
% Create custom solver
solver = ODESolver('ode23s');
solver.setOptions('RelTol', 1e-10, 'AbsTol', 1e-20);

% Run with custom solver
result = sim.run('solver_options', solver.options);
```

## Key Benefits

1. **Modularity**: Each class has a single, well-defined responsibility
2. **Reusability**: Components can be used independently or in different combinations
3. **Maintainability**: Clear separation of concerns makes code easier to understand and modify
4. **Flexibility**: Easy to modify parameters, initial conditions, and solver options
5. **Extensibility**: New kernel types, initial conditions, or analysis methods can be added easily
6. **Validation**: Built-in validation ensures consistency and catches errors early

## Migration Notes

- **Mathematical Formulations**: All original mathematical formulas are preserved exactly
- **Matrix Structures**: Original matrix shapes and operations are maintained
- **Integration Logic**: ODE solver interface and Jacobian computation are identical
- **Output Formats**: Generated plots and data structures match the original output

## Future Enhancements

1. **Tracer Integration**: Extend framework to support radioactive tracer calculations
2. **Additional Kernels**: Implement more sophisticated coagulation kernels
3. **Parallel Processing**: Add support for parallel computation of beta matrices
4. **Advanced Diagnostics**: Implement more sophisticated mass balance and conservation checks
5. **Parameter Estimation**: Add framework for parameter optimization and sensitivity analysis

## File Organization

```
particle_aggregation/
├── SimulationConfig.m          # Configuration management
├── DerivedGrid.m               # Derived constants and grids
├── KernelLibrary.m             # Coagulation kernels
├── SettlingVelocityService.m   # Settling velocity calculations
├── LinearProcessBuilder.m      # Linear process operators
├── InitialSpectrumBuilder.m    # Initial condition generation
├── BetaAssembler.m             # Coagulation matrix computation
├── BetaMatrices.m              # Beta matrix container
├── CoagulationRHS.m            # ODE right-hand side
├── ODESolver.m                 # Time integration
├── MassBalanceAnalyzer.m       # Diagnostics and analysis
├── OutputGenerator.m            # Visualization and export
├── CoagulationSimulation.m     # Main simulation controller
├── example_usage.m             # Usage examples
├── README_OOP.md               # This documentation

```

## Dependencies

The OOP framework requires MATLAB with the following toolboxes:
- **Base MATLAB**: Core functionality
- **Optimization Toolbox**: For potential future steady-state calculations
- **Statistics and Machine Learning Toolbox**: For potential future analysis features

## Testing

Run the example script to verify the framework:
```matlab
example_usage
```

This will demonstrate all major features and validate the implementation.
