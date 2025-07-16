# Examples Directory

This directory contains example scripts demonstrating how to use the Shilling ODE solver library. Each example is self-contained and can be run independently.

## Available Examples

### Basic Usage Examples

- **`example_simple_ode.m`** - Simple first-order ODE with analytical solution
  - Demonstrates basic usage of different solvers
  - Compares Euler, RK4, and RK5 methods
  - Shows error analysis and performance comparison

- **`example_validation.m`** - Validation suite for the solver library
  - Tests various solvers with known analytical solutions
  - Validates both direct solver calls and odeManager interface
  - Includes error handling tests

### Advanced System Examples

- **`example_van_der_pol.m`** - Van der Pol oscillator (system of ODEs)
  - Demonstrates solving systems of equations
  - Shows phase portrait analysis
  - Compares different solver performance
  - Parameter study with different damping values

- **`example_lorenz.m`** - Lorenz chaotic system
  - Complex 3D dynamical system
  - Demonstrates high-order solver usage
  - Shows sensitivity to initial conditions
  - 3D visualization and parameter studies

## How to Run Examples

1. **Add the library to your MATLAB path:**
   ```matlab
   addpath('../src');  % Add the src directory to path
   ```

2. **Run any example:**
   ```matlab
   run('example_simple_ode.m');
   ```

3. **Or run individual sections:**
   ```matlab
   open('example_van_der_pol.m');  % Open in editor and run sections
   ```

## Example Categories

### Solver Comparison
- `example_simple_ode.m` - Basic comparison of different methods
- `example_validation.m` - Systematic validation of all solvers

### Physical Systems
- `example_van_der_pol.m` - Nonlinear oscillator
- `example_lorenz.m` - Chaotic dynamics

### Mathematical Problems
- Simple ODEs with analytical solutions
- Systems of ODEs
- Nonlinear dynamics
- Chaotic systems

## Usage Patterns

### Basic ODE Solving
```matlab
% Define ODE
f = @(t, y) -2*y;

% Set up problem
tspan = [0, 5];
y0 = 1;

% Solve with default settings
[t, y] = odeManager(f, tspan, y0);

% Or with custom options
options = odeOptions('Solver', 'rk5', 'StepSize', 0.01);
[t, y] = odeManager(f, tspan, y0, options);
```

### System of ODEs
```matlab
% Define system (Van der Pol oscillator)
f = @(t, x) [x(2); (1-x(1)^2)*x(2) - x(1)];

% Initial conditions
y0 = [2; 0];

% Solve system
[t, y] = odeManager(f, tspan, y0, 'Solver', 'rk4');
```

### Direct Solver Calls
```matlab
% Call solver directly
[t, y] = odeRK4(f, tspan, y0, 0.01);
```

## Tips for Using Examples

1. **Start with simple examples** - Begin with `example_simple_ode.m`
2. **Check system requirements** - Some examples require specific MATLAB versions
3. **Modify parameters** - Change step sizes, time spans, and initial conditions
4. **Add your own problems** - Use examples as templates for your own ODEs
5. **Compare solvers** - Try different methods to see performance differences

## Common Issues

- **Path not set**: Make sure to add the `src/` directory to your MATLAB path
- **Step size too large**: Reduce step size if solution appears unstable
- **Memory issues**: Use smaller time spans or larger step sizes for long simulations
- **Plotting issues**: Some examples create multiple figures - close unused ones

## Further Reading

- Main documentation: `../README.md`
- API reference: Check function help using `help odeManager`
- Solver-specific help: `help odeRK4`, `help odeEuler`, etc.