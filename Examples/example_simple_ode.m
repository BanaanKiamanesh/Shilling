%% Simple ODE Example
% This example demonstrates how to solve a simple first-order ODE
% using the Shilling ODE solver library.

clear; close all; clc;

%% Define the ODE
% Solve dy/dt = -2*y + 1 with y(0) = 0
% Analytical solution: y(t) = 0.5*(1 - exp(-2*t))

f = @(t, y) -2*y + 1;

%% Set up the problem
tspan = [0, 5];     % Time interval
y0 = 0;             % Initial condition

%% Solve with different methods
fprintf('Solving dy/dt = -2*y + 1 with y(0) = 0\n');
fprintf('Time interval: [0, 5]\n');
fprintf('Step size: 0.1\n\n');

% Euler method
[t1, y1] = odeManager(f, tspan, y0, 'Solver', 'euler', 'StepSize', 0.1);
fprintf('Euler method: %d time steps\n', length(t1));

% RK4 method
[t2, y2] = odeManager(f, tspan, y0, 'Solver', 'rk4', 'StepSize', 0.1);
fprintf('RK4 method: %d time steps\n', length(t2));

% RK5 method
[t3, y3] = odeManager(f, tspan, y0, 'Solver', 'rk5', 'StepSize', 0.1);
fprintf('RK5 method: %d time steps\n', length(t3));

%% Analytical solution for comparison
t_exact = 0:0.01:5;
y_exact = 0.5 * (1 - exp(-2*t_exact));

%% Plot comparison
figure('Name', 'Simple ODE - Method Comparison', 'Position', [100, 100, 800, 600]);

subplot(2,1,1);
plot(t_exact, y_exact, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t1, y1, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Euler', 'MarkerSize', 4);
plot(t2, y2, 'b:', 'LineWidth', 1.5, 'DisplayName', 'RK4', 'MarkerSize', 4);
plot(t3, y3, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'RK5', 'MarkerSize', 4);
xlabel('Time');
ylabel('y(t)');
title('Solution Comparison: dy/dt = -2*y + 1');
legend('Location', 'best');
grid on;

% Error plot
subplot(2,1,2);
% Interpolate exact solution at solver time points
y_exact_interp1 = interp1(t_exact, y_exact, t1, 'linear');
y_exact_interp2 = interp1(t_exact, y_exact, t2, 'linear');
y_exact_interp3 = interp1(t_exact, y_exact, t3, 'linear');

semilogy(t1, abs(y1 - y_exact_interp1), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Euler Error');
hold on;
semilogy(t2, abs(y2 - y_exact_interp2), 'b:', 'LineWidth', 1.5, 'DisplayName', 'RK4 Error');
semilogy(t3, abs(y3 - y_exact_interp3), 'g-.', 'LineWidth', 1.5, 'DisplayName', 'RK5 Error');
xlabel('Time');
ylabel('Absolute Error');
title('Error Analysis');
legend('Location', 'best');
grid on;

%% Print final values
fprintf('\nFinal values at t = 5:\n');
fprintf('Analytical: %.6f\n', y_exact(end));
fprintf('Euler:      %.6f (error: %.2e)\n', y1(end), abs(y1(end) - y_exact(end)));
fprintf('RK4:        %.6f (error: %.2e)\n', y2(end), abs(y2(end) - y_exact(end)));
fprintf('RK5:        %.6f (error: %.2e)\n', y3(end), abs(y3(end) - y_exact(end)));