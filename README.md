# Shilling - MATLAB ODE Solver Collection

A comprehensive collection of fixed-step ODE solvers for nth-order non-linear ordinary differential equations, all implemented in MATLAB. This library provides a unified interface to access various Runge-Kutta and other numerical integration methods.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Basic Usage](#basic-usage)
  - [Advanced Usage](#advanced-usage)
  - [Options Configuration](#options-configuration)
- [Available Solvers](#available-solvers)
- [Examples](#examples)
- [API Reference](#api-reference)
- [Contributing](#contributing)
- [License](#license)

## Installation

1. Clone or download this repository
2. Add the `src/` directory to your MATLAB path:
   ```matlab
   addpath('path/to/shilling/src');
   ```
3. You're ready to use the solvers!

## Quick Start

```matlab
% Define a simple ODE: dy/dt = -2*y + 1
f = @(t, y) -2*y + 1;

% Set up the problem
tspan = [0, 5];     % Time interval
y0 = 0;             % Initial condition
options = odeOptions('Solver', 'rk4', 'StepSize', 0.01);

% Solve the ODE
[t, y] = odeManager(f, tspan, y0, options);

% Plot results
plot(t, y);
xlabel('Time');
ylabel('y(t)');
title('Solution of dy/dt = -2*y + 1');
```

## Usage

### Basic Usage

The library provides two main interfaces:

1. **High-level interface** via `odeManager` (recommended)
2. **Direct solver calls** for specific methods

#### Using odeManager (Recommended)

```matlab
% Basic usage with default settings (RK4 solver)
[t, y] = odeManager(f, tspan, y0);

% With custom options
options = odeOptions('Solver', 'rk5', 'StepSize', 0.001);
[t, y] = odeManager(f, tspan, y0, options);

% With inline options
[t, y] = odeManager(f, tspan, y0, 'Solver', 'rk4', 'StepSize', 0.01);
```

#### Direct Solver Calls

```matlab
% Call solver directly
[t, y] = odeRK4(f, tspan, y0, 0.01);
```

### Advanced Usage

#### System of ODEs

```matlab
% Van der Pol oscillator: [dx/dt; dy/dt] = [y; (1-x^2)*y - x]
f = @(t, x) [x(2); (1-x(1)^2)*x(2) - x(1)];
tspan = [0, 20];
y0 = [2; 0];

options = odeOptions('Solver', 'rk4', 'StepSize', 0.01);
[t, y] = odeManager(f, tspan, y0, options);

% Plot phase portrait
figure;
plot(y(:,1), y(:,2));
xlabel('x');
ylabel('dx/dt');
title('Van der Pol Oscillator - Phase Portrait');
```

#### Higher-Order ODEs

```matlab
% Second-order ODE: d²y/dt² + 2*dy/dt + y = 0
% Convert to system: [y; dy/dt] = [y2; -2*y2 - y1]
f = @(t, x) [x(2); -2*x(2) - x(1)];
tspan = [0, 10];
y0 = [1; 0];  % y(0) = 1, dy/dt(0) = 0

[t, y] = odeManager(f, tspan, y0, 'Solver', 'rk5', 'StepSize', 0.01);
```

### Options Configuration

Use `odeOptions` to configure solver behavior:

```matlab
% Create options with default values
options = odeOptions();

% Create options with specific values
options = odeOptions('Solver', 'rk4', 'StepSize', 0.001);

% Modify existing options
options = odeOptions(options, 'Solver', 'rk5');

% Available options:
% - 'Solver': Solver name (string) or function handle [default: 'rk4']
% - 'StepSize': Fixed time-step size [default: 0.01]
% - 'RelTol': Relative tolerance [default: 1e-6]
% - 'AbsTol': Absolute tolerance [default: 1e-9]
% - 'MaxStep': Maximum step size [default: []]
% - 'UserData': Arbitrary user data [default: []]
```

## Available Solvers

### Fixed Step Solvers

| Name | Description | Properties | Order | Stages | Registers | CFL | Reference |
|------|-------------|-----------|--------|--------|-----------|-----|-----------|
| `euler` | Euler | | 1 | 1 | 1 | 1.0 | [Euler (1768)](https://archive.org/details/institutionescal020326mbp) |
| `midpoint` | Midpoint | | 2 | 2 | 2 | | Classic method |
| `heun` | Heun | | 2 | 2 | 2 | | Classic method |
| `rkssp22` | 2-stage, 2nd order TVD Runge-Kutta Shu-Osher | Strong stability preserving | 2 | 2 | 1 | 1.0 | [Shu & Osher (1988)](https://ntrs.nasa.gov/api/citations/19880014833/downloads/19880014833.pdf) |
| `rk3` | 3rd order Runge-Kutta | | 3 | 3 | 3 | | Classic method |
| `rkssp33` | 3-stage, 3rd order TVD Runge-Kutta Shu-Osher | Strong stability preserving | 3 | 3 | 1 | 1.0 | [Shu & Osher (1988)](https://ntrs.nasa.gov/api/citations/19880014833/downloads/19880014833.pdf) |
| `rkssp53` | 5-stage, 3rd order SSP Runge-Kutta Spiteri-Ruuth | Strong stability preserving | 3 | 5 | 2 | 2.65 | [Ruuth (2006)](https://www.ams.org/journals/mcom/2006-75-253/S0025-5718-05-01772-2/S0025-5718-05-01772-2.pdf) |
| `rk4` | Classic 4th order Runge-Kutta | | 4 | 4 | 4 | | [Kutta (1901)](https://archive.org/stream/zeitschriftfrma12runggoog#page/n449/mode/2up) |
| `rks4` | 4th order Runge-Kutta Shanks | | 4 | 4 | 4 | | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf) |
| `rkr4` | 4th order Runge-Kutta Ralston | | 4 | 4 | 4 | | [Ralston (1962)](https://doi.org/10.1090%2FS0025-5718-1962-0150954-0) |
| `rkls44` | 4-stage, 4th order low storage non-TVD Runge-Kutta Jiang-Shu | Low storage | 4 | 4 | 2 | | [Jiang and Shu (1988)](https://ntrs.nasa.gov/api/citations/19960007052/downloads/19960007052.pdf) |
| `rkls54` | 5-stage, 4th order low storage Runge-Kutta Carpenter-Kennedy | Low storage | 4 | 5 | 2 | 0.32 | [Carpenter & Kennedy (1994)](https://ntrs.nasa.gov/api/citations/19940028444/downloads/19940028444.pdf) |
| `rkssp54` | 5-stage, 4th order SSP Runge-Kutta Spiteri-Ruuth | Strong stability preserving | 4 | 5 | 4 | 1.51 | [Ruuth (2006)](https://www.ams.org/journals/mcom/2006-75-253/S0025-5718-05-01772-2/S0025-5718-05-01772-2.pdf) |
| `rks5` | 5th order Runge-Kutta Shanks | | 5 | 5 | 5 | | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf) |
| `rk5` | 5th order Runge-Kutta | | 5 | 6 | 6 | | Classic method |
| `rkc5` | 5th order Runge-Kutta Cassity | | 5 | 6 | 6 | | [Cassity (1966)](https://epubs.siam.org/doi/10.1137/0703052) |
| `rkl5` | 5th order Runge-Kutta Lawson | | 5 | 6 | 6 | | [Lawson (1966)](https://epubs.siam.org/doi/abs/10.1137/0703051) |
| `rklk5a` | 5th order Runge-Kutta Luther-Konen 1 | | 5 | 6 | 6 | | [Luther & Konen (1965)](https://epubs.siam.org/doi/abs/10.1137/1007112) |
| `rklk5b` | 5th order Runge-Kutta Luther-Konen 2 | | 5 | 6 | 6 | | [Luther & Konen (1965)](https://epubs.siam.org/doi/abs/10.1137/1007112) |
| `rkb6` | 6th order Runge-Kutta Butcher | | 6 | 7 | 7 | | [Butcher (1963)](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/40DFE501CAB781C9AAE1439B6B8F481A/S1446788700023387a.pdf/div-class-title-on-runge-kutta-processes-of-high-order-div.pdf) |
| `rk7` | 7th order Runge-Kutta Shanks | | 7 | 9 | 9 | | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf) |
| `rk8_10` | 10-stage, 8th order Runge-Kutta Shanks | | 8 | 10 | 10 | | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf) |
| `rkcv8` | 11-stage, 8th order Runge-Kutta Cooper-Verner | | 8 | 11 | 11 | | [Cooper & Verner (1972)](https://epubs.siam.org/doi/abs/10.1137/0709037) |
| `rk8_12` | 12-stage, 8th order Runge-Kutta Shanks | | 8 | 12 | 12 | | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf) |
| `rkz10` | 10th order Runge-Kutta Zhang | | 10 | 16 | 16 | | [Zhang (2019)](https://arxiv.org/abs/1911.00318) |
| `rko10` | 10th order Runge-Kutta Ono | | 10 | 17 | 17 | | [Ono (2003)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK10/RKcoeff10f_1.pdf) |
| `rkh10` | 10th order Runge-Kutta Hairer | | 10 | 17 | 17 | | [Hairer (1978)](https://www.researchgate.net/publication/31221486_A_Runge-Kutta_Method_of_Order_10) |

## Examples

### Example 1: Simple First-Order ODE

```matlab
% Solve dy/dt = -2*y + 1 with y(0) = 0
f = @(t, y) -2*y + 1;
tspan = [0, 5];
y0 = 0;

% Solve with different methods
[t1, y1] = odeManager(f, tspan, y0, 'Solver', 'euler', 'StepSize', 0.1);
[t2, y2] = odeManager(f, tspan, y0, 'Solver', 'rk4', 'StepSize', 0.1);
[t3, y3] = odeManager(f, tspan, y0, 'Solver', 'rk5', 'StepSize', 0.1);

% Analytical solution for comparison
t_exact = 0:0.01:5;
y_exact = 0.5 * (1 - exp(-2*t_exact));

% Plot comparison
figure;
plot(t_exact, y_exact, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t1, y1, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Euler');
plot(t2, y2, 'b:', 'LineWidth', 1.5, 'DisplayName', 'RK4');
plot(t3, y3, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'RK5');
xlabel('Time');
ylabel('y(t)');
title('Comparison of ODE Solvers');
legend;
grid on;
```

### Example 2: Van der Pol Oscillator

```matlab
% Van der Pol oscillator: d²x/dt² - μ(1-x²)dx/dt + x = 0
% Convert to system: [x; dx/dt] = [x2; μ(1-x1²)x2 - x1]
mu = 1;
f = @(t, x) [x(2); mu*(1-x(1)^2)*x(2) - x(1)];
tspan = [0, 20];
y0 = [2; 0];

% Solve with high-order method
options = odeOptions('Solver', 'rk5', 'StepSize', 0.01);
[t, y] = odeManager(f, tspan, y0, options);

% Plot results
figure;
subplot(2,1,1);
plot(t, y(:,1), 'b-', 'LineWidth', 2);
xlabel('Time');
ylabel('x(t)');
title('Van der Pol Oscillator - Time Series');
grid on;

subplot(2,1,2);
plot(y(:,1), y(:,2), 'r-', 'LineWidth', 2);
xlabel('x');
ylabel('dx/dt');
title('Van der Pol Oscillator - Phase Portrait');
grid on;
```

### Example 3: Lorenz System

```matlab
% Lorenz system: chaotic dynamics
sigma = 10;
rho = 28;
beta = 8/3;

f = @(t, x) [sigma*(x(2)-x(1)); 
             x(1)*(rho-x(3))-x(2); 
             x(1)*x(2)-beta*x(3)];

tspan = [0, 30];
y0 = [1; 1; 1];

% Use high-order solver for accuracy
options = odeOptions('Solver', 'rk8_10', 'StepSize', 0.001);
[t, y] = odeManager(f, tspan, y0, options);

% 3D phase space plot
figure;
plot3(y(:,1), y(:,2), y(:,3), 'b-', 'LineWidth', 1);
xlabel('x');
ylabel('y');
zlabel('z');
title('Lorenz Attractor');
grid on;
```

## API Reference

### Main Functions

#### `odeManager(f, tspan, y0, options)`
High-level interface to all ODE solvers.

**Parameters:**
- `f`: Function handle `@(t, y)` representing dy/dt = f(t, y)
- `tspan`: Time span as `[t0, tf]`
- `y0`: Initial condition (column vector)
- `options`: Options structure from `odeOptions()` (optional)

**Returns:**
- `t`: Time vector
- `y`: Solution matrix (each row is a time point, each column is a state variable)

#### `odeOptions(varargin)`
Create or modify options structure for ODE solvers.

**Parameters:**
- Name-value pairs or existing options structure

**Returns:**
- `options`: Options structure

### Individual Solver Functions

All solver functions have the same interface:
```matlab
[t, y] = solverName(f, tspan, y0, h)
```

Where:
- `f`: Function handle for the ODE
- `tspan`: Time span `[t0, tf]`
- `y0`: Initial condition
- `h`: Step size

## Contributing

1. Fork the repository
2. Create a feature branch
3. Add your solver following the existing pattern
4. Update documentation
5. Submit a pull request

### Adding a New Solver

1. Create a new file `src/odeNewSolver.m`
2. Follow the existing function signature and help format
3. Add the solver to the `solverMap` in `odeManager.m`
4. Update the README solver table
5. Add tests and examples

## License

This project is open source. Please check the license file for details.

## Citation

If you use this library in academic work, please cite:

```
@software{shilling_ode_solvers,
  title = {Shilling: A Collection of ODE Solvers for MATLAB},
  author = {[Author Name]},
  year = {2024},
  url = {https://github.com/[username]/shilling}
}
```

