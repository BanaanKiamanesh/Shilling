%% Van der Pol Oscillator Example
% This example demonstrates how to solve a system of ODEs
% using the Van der Pol oscillator as a test case.

clear; close all; clc;

%% Define the Van der Pol Oscillator
% Van der Pol oscillator: d²x/dt² - μ(1-x²)dx/dt + x = 0
% Convert to system: [x; dx/dt] = [x2; μ(1-x1²)x2 - x1]

mu = 1;  % Damping parameter
f = @(t, x) [x(2); mu*(1-x(1)^2)*x(2) - x(1)];

%% Set up the problem
tspan = [0, 20];    % Time interval
y0 = [2; 0];        % Initial condition: x(0) = 2, dx/dt(0) = 0

%% Solve with high-order method
fprintf('Solving Van der Pol oscillator with μ = %.1f\n', mu);
fprintf('Initial condition: x(0) = %.1f, dx/dt(0) = %.1f\n', y0(1), y0(2));
fprintf('Time interval: [0, 20]\n');
fprintf('Step size: 0.01\n\n');

options = odeOptions('Solver', 'rk5', 'StepSize', 0.01);
tic;
[t, y] = odeManager(f, tspan, y0, options);
solve_time = toc;

fprintf('Solution completed in %.3f seconds\n', solve_time);
fprintf('Number of time steps: %d\n', length(t));

%% Plot results
figure('Name', 'Van der Pol Oscillator', 'Position', [100, 100, 1200, 800]);

% Time series plot
subplot(2,2,1);
plot(t, y(:,1), 'b-', 'LineWidth', 2);
xlabel('Time');
ylabel('x(t)');
title('Position vs Time');
grid on;

subplot(2,2,2);
plot(t, y(:,2), 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('dx/dt');
title('Velocity vs Time');
grid on;

% Phase portrait
subplot(2,2,3);
plot(y(:,1), y(:,2), 'g-', 'LineWidth', 2);
xlabel('x');
ylabel('dx/dt');
title('Phase Portrait');
grid on;
axis equal;

% Combined plot
subplot(2,2,4);
plot(t, y(:,1), 'b-', 'LineWidth', 2, 'DisplayName', 'x(t)');
hold on;
plot(t, y(:,2), 'r-', 'LineWidth', 2, 'DisplayName', 'dx/dt');
xlabel('Time');
ylabel('State Variables');
title('Combined Time Series');
legend('Location', 'best');
grid on;

%% Compare different solvers
fprintf('\nComparing different solvers:\n');
solvers = {'rk4', 'rk5', 'rk7', 'rk8_10'};
colors = {'b', 'r', 'g', 'm'};
figure('Name', 'Solver Comparison', 'Position', [150, 150, 1000, 400]);

for i = 1:length(solvers)
    fprintf('Solving with %s... ', solvers{i});
    tic;
    [t_comp, y_comp] = odeManager(f, tspan, y0, 'Solver', solvers{i}, 'StepSize', 0.01);
    time_comp = toc;
    fprintf('Done in %.3f seconds\n', time_comp);
    
    subplot(1,2,1);
    plot(t_comp, y_comp(:,1), colors{i}, 'LineWidth', 1.5, 'DisplayName', upper(solvers{i}));
    hold on;
    
    subplot(1,2,2);
    plot(y_comp(:,1), y_comp(:,2), colors{i}, 'LineWidth', 1.5, 'DisplayName', upper(solvers{i}));
    hold on;
end

subplot(1,2,1);
xlabel('Time');
ylabel('x(t)');
title('Position: Solver Comparison');
legend('Location', 'best');
grid on;

subplot(1,2,2);
xlabel('x');
ylabel('dx/dt');
title('Phase Portrait: Solver Comparison');
legend('Location', 'best');
grid on;
axis equal;

%% Different damping parameters
fprintf('\nAnalyzing different damping parameters:\n');
mu_values = [0.1, 0.5, 1.0, 2.0];
figure('Name', 'Damping Parameter Study', 'Position', [200, 200, 1000, 600]);

for i = 1:length(mu_values)
    mu_current = mu_values(i);
    f_current = @(t, x) [x(2); mu_current*(1-x(1)^2)*x(2) - x(1)];
    
    fprintf('μ = %.1f... ', mu_current);
    [t_mu, y_mu] = odeManager(f_current, tspan, y0, 'Solver', 'rk5', 'StepSize', 0.01);
    fprintf('Done\n');
    
    subplot(2,2,i);
    plot(y_mu(:,1), y_mu(:,2), 'b-', 'LineWidth', 2);
    xlabel('x');
    ylabel('dx/dt');
    title(sprintf('μ = %.1f', mu_current));
    grid on;
    axis equal;
end