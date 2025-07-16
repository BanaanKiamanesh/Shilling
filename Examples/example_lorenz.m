%% Lorenz System Example
% This example demonstrates solving the Lorenz system,
% a famous chaotic dynamical system.

clear; close all; clc;

%% Define the Lorenz System
% dx/dt = σ(y - x)
% dy/dt = x(ρ - z) - y  
% dz/dt = xy - βz

% Standard parameters for chaotic behavior
sigma = 10;
rho = 28;
beta = 8/3;

f = @(t, x) [sigma*(x(2)-x(1)); 
             x(1)*(rho-x(3))-x(2); 
             x(1)*x(2)-beta*x(3)];

%% Set up the problem
tspan = [0, 30];    % Time interval
y0 = [1; 1; 1];     % Initial condition

%% Solve with high-order solver for accuracy
fprintf('Solving Lorenz system with parameters:\n');
fprintf('σ = %.1f, ρ = %.1f, β = %.3f\n', sigma, rho, beta);
fprintf('Initial condition: [%.1f, %.1f, %.1f]\n', y0(1), y0(2), y0(3));
fprintf('Time interval: [0, 30]\n');
fprintf('Step size: 0.001\n\n');

options = odeOptions('Solver', 'rk8_10', 'StepSize', 0.001);
tic;
[t, y] = odeManager(f, tspan, y0, options);
solve_time = toc;

fprintf('Solution completed in %.3f seconds\n', solve_time);
fprintf('Number of time steps: %d\n', length(t));

%% 3D Phase Space Plot
figure('Name', 'Lorenz Attractor', 'Position', [100, 100, 1200, 800]);

% Main 3D plot
subplot(2,2,[1,3]);
plot3(y(:,1), y(:,2), y(:,3), 'b-', 'LineWidth', 0.8);
xlabel('x');
ylabel('y');
zlabel('z');
title('Lorenz Attractor - 3D Phase Space');
grid on;
view(45, 30);

% 2D projections
subplot(2,2,2);
plot(y(:,1), y(:,2), 'r-', 'LineWidth', 0.8);
xlabel('x');
ylabel('y');
title('x-y Projection');
grid on;
axis equal;

subplot(2,2,4);
plot(y(:,1), y(:,3), 'g-', 'LineWidth', 0.8);
xlabel('x');
ylabel('z');
title('x-z Projection');
grid on;
axis equal;

%% Time Series Analysis
figure('Name', 'Lorenz Time Series', 'Position', [150, 150, 1000, 600]);

subplot(3,1,1);
plot(t, y(:,1), 'b-', 'LineWidth', 1);
xlabel('Time');
ylabel('x(t)');
title('x component');
grid on;

subplot(3,1,2);
plot(t, y(:,2), 'r-', 'LineWidth', 1);
xlabel('Time');
ylabel('y(t)');
title('y component');
grid on;

subplot(3,1,3);
plot(t, y(:,3), 'g-', 'LineWidth', 1);
xlabel('Time');
ylabel('z(t)');
title('z component');
grid on;

%% Sensitivity to Initial Conditions
fprintf('\nAnalyzing sensitivity to initial conditions:\n');
epsilon = 1e-6;  % Small perturbation
y0_pert = y0 + epsilon * [1; 0; 0];  % Perturb x coordinate

fprintf('Perturbing initial condition by ε = %.1e\n', epsilon);
fprintf('Original IC: [%.6f, %.6f, %.6f]\n', y0(1), y0(2), y0(3));
fprintf('Perturbed IC: [%.6f, %.6f, %.6f]\n', y0_pert(1), y0_pert(2), y0_pert(3));

% Solve perturbed system
[t_pert, y_pert] = odeManager(f, tspan, y0_pert, options);

% Calculate distance between trajectories
distance = sqrt(sum((y - y_pert).^2, 2));

figure('Name', 'Sensitivity Analysis', 'Position', [200, 200, 1000, 600]);

subplot(2,1,1);
semilogy(t, distance, 'k-', 'LineWidth', 2);
xlabel('Time');
ylabel('Distance');
title('Distance Between Trajectories');
grid on;

subplot(2,1,2);
plot3(y(:,1), y(:,2), y(:,3), 'b-', 'LineWidth', 1, 'DisplayName', 'Original');
hold on;
plot3(y_pert(:,1), y_pert(:,2), y_pert(:,3), 'r-', 'LineWidth', 1, 'DisplayName', 'Perturbed');
xlabel('x');
ylabel('y');
zlabel('z');
title('Trajectory Comparison');
legend('Location', 'best');
grid on;
view(45, 30);

%% Parameter Study
fprintf('\nParameter study - varying ρ:\n');
rho_values = [24, 28, 32, 36];
figure('Name', 'Parameter Study', 'Position', [250, 250, 1000, 600]);

for i = 1:length(rho_values)
    rho_current = rho_values(i);
    f_current = @(t, x) [sigma*(x(2)-x(1)); 
                         x(1)*(rho_current-x(3))-x(2); 
                         x(1)*x(2)-beta*x(3)];
    
    fprintf('ρ = %.1f... ', rho_current);
    [t_param, y_param] = odeManager(f_current, [0, 15], y0, 'Solver', 'rk5', 'StepSize', 0.001);
    fprintf('Done\n');
    
    subplot(2,2,i);
    plot3(y_param(:,1), y_param(:,2), y_param(:,3), 'b-', 'LineWidth', 1);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(sprintf('ρ = %.1f', rho_current));
    grid on;
    view(45, 30);
end

%% Statistics
fprintf('\nSolution statistics:\n');
fprintf('Final point: [%.3f, %.3f, %.3f]\n', y(end,1), y(end,2), y(end,3));
fprintf('Range of x: [%.3f, %.3f]\n', min(y(:,1)), max(y(:,1)));
fprintf('Range of y: [%.3f, %.3f]\n', min(y(:,2)), max(y(:,2)));
fprintf('Range of z: [%.3f, %.3f]\n', min(y(:,3)), max(y(:,3)));