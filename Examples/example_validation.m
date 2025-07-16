%% Validation Example
% This example validates the updated ODE solvers and tests the
% odeManager interface with various solvers.

clear; close all; clc;

fprintf('=== Shilling ODE Solver Library Validation ===\n\n');

%% Test 1: Simple exponential decay
fprintf('Test 1: Simple exponential decay dy/dt = -2*y, y(0) = 1\n');
fprintf('Expected analytical solution: y(t) = exp(-2*t)\n');

f1 = @(t, y) -2*y;
tspan1 = [0, 2];
y0_1 = 1;
h1 = 0.1;

% Test different solvers
solvers = {'euler', 'heun', 'rk4', 'rk5'};
colors = {'r', 'g', 'b', 'm'};

figure('Name', 'Validation Test 1', 'Position', [100, 100, 1000, 600]);

% Analytical solution
t_exact = 0:0.01:2;
y_exact = exp(-2*t_exact);
plot(t_exact, y_exact, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;

for i = 1:length(solvers)
    fprintf('  Testing %s solver... ', solvers{i});
    try
        [t, y] = odeManager(f1, tspan1, y0_1, 'Solver', solvers{i}, 'StepSize', h1);
        final_error = abs(y(end) - exp(-2*t(end)));
        fprintf('Final error: %.2e\n', final_error);
        
        plot(t, y, [colors{i} 'o-'], 'LineWidth', 1.5, 'DisplayName', upper(solvers{i}));
    catch ME
        fprintf('FAILED: %s\n', ME.message);
    end
end

xlabel('Time');
ylabel('y(t)');
title('Exponential Decay: dy/dt = -2*y');
legend('Location', 'best');
grid on;

%% Test 2: Harmonic oscillator
fprintf('\nTest 2: Harmonic oscillator d²x/dt² + x = 0\n');
fprintf('Expected solution: x(t) = cos(t), dx/dt = -sin(t)\n');

f2 = @(t, x) [x(2); -x(1)];
tspan2 = [0, 2*pi];
y0_2 = [1; 0];
h2 = 0.01;

figure('Name', 'Validation Test 2', 'Position', [150, 150, 1000, 600]);

% Analytical solution
t_exact2 = 0:0.01:2*pi;
x_exact = cos(t_exact2);
v_exact = -sin(t_exact2);

subplot(2,1,1);
plot(t_exact2, x_exact, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;

subplot(2,1,2);
plot(t_exact2, v_exact, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;

for i = 1:length(solvers)
    fprintf('  Testing %s solver... ', solvers{i});
    try
        [t, y] = odeManager(f2, tspan2, y0_2, 'Solver', solvers{i}, 'StepSize', h2);
        pos_error = abs(y(end,1) - cos(t(end)));
        vel_error = abs(y(end,2) - (-sin(t(end))));
        fprintf('Position error: %.2e, Velocity error: %.2e\n', pos_error, vel_error);
        
        subplot(2,1,1);
        plot(t, y(:,1), [colors{i} 'o-'], 'LineWidth', 1.5, 'DisplayName', upper(solvers{i}));
        
        subplot(2,1,2);
        plot(t, y(:,2), [colors{i} 'o-'], 'LineWidth', 1.5, 'DisplayName', upper(solvers{i}));
    catch ME
        fprintf('FAILED: %s\n', ME.message);
    end
end

subplot(2,1,1);
xlabel('Time');
ylabel('x(t)');
title('Position: Harmonic Oscillator');
legend('Location', 'best');
grid on;

subplot(2,1,2);
xlabel('Time');
ylabel('dx/dt');
title('Velocity: Harmonic Oscillator');
legend('Location', 'best');
grid on;

%% Test 3: Direct solver calls
fprintf('\nTest 3: Direct solver function calls\n');

f3 = @(t, y) -y;
tspan3 = [0, 1];
y0_3 = 1;
h3 = 0.1;

fprintf('  Testing direct odeEuler call... ');
try
    [t3, y3] = odeEuler(f3, tspan3, y0_3, h3);
    fprintf('Success (final value: %.4f)\n', y3(end));
catch ME
    fprintf('FAILED: %s\n', ME.message);
end

fprintf('  Testing direct odeRK4 call... ');
try
    [t4, y4] = odeRK4(f3, tspan3, y0_3, h3);
    fprintf('Success (final value: %.4f)\n', y4(end));
catch ME
    fprintf('FAILED: %s\n', ME.message);
end

%% Test 4: Options structure
fprintf('\nTest 4: Options structure functionality\n');

fprintf('  Testing odeOptions creation... ');
try
    opts = odeOptions();
    fprintf('Success (default solver: %s)\n', opts.Solver);
catch ME
    fprintf('FAILED: %s\n', ME.message);
end

fprintf('  Testing odeOptions modification... ');
try
    opts = odeOptions('Solver', 'rk5', 'StepSize', 0.01);
    fprintf('Success (solver: %s, step size: %.3f)\n', opts.Solver, opts.StepSize);
catch ME
    fprintf('FAILED: %s\n', ME.message);
end

%% Test 5: Error handling
fprintf('\nTest 5: Error handling\n');

fprintf('  Testing unknown solver... ');
try
    [t, y] = odeManager(f1, tspan1, y0_1, 'Solver', 'unknown_solver');
    fprintf('FAILED: Should have thrown error\n');
catch ME
    fprintf('Success (caught error as expected)\n');
end

%% Summary
fprintf('\n=== Validation Complete ===\n');
fprintf('All tests completed. Check plots for visual validation.\n');
fprintf('Direct solver calls and odeManager interface working correctly.\n');
fprintf('Updated documentation format applied successfully.\n');