% VanDerPolExample.m — Shilling demonstration script
%
% Solves the Van-der-Pol oscillator using the high-level helper `odeManager`
% with a 5-th order Runge–Kutta solver and plots the result.

addpath(genpath('..'));   % ensure the Shilling src/ folder is on the MATLAB path

% Van-der-Pol oscillator
f      = @(t,x) [x(2); (1-x(1).^2).*x(2) - x(1)];
tspan  = [0 20];
y0     = [2; 0];

% Choose solver and step size
opt = odeOptions('Solver','rk5', 'StepSize', 1e-2);

% Integrate
[T, Y] = odeManager(f, tspan, y0, opt);

% Plot the two state variables
figure;
plot(T, Y, 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('States');
legend({'x_1','x_2'});

title('Van-der-Pol oscillator — 5^{th}-order Runge–Kutta');
grid on;