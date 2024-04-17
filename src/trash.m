clear
close all
clc

f = @(t, x) [x(2); (1 - x(1)^2)*x(2) - x(1)]; % Van der Pol ODE
TSpan= [0, 20];    % Time span
Y0= [2; 2];        % Initial condition
h = 0.01;          % Step size
[T, Y] = odeRKLS44(f, TSpan, Y0, h);

% Plot results
figure;
plot(T, Y, 'LineWidth', 2);
xlabel('Time [s]');
ylabel('State Variables');
title('Van der Pol Equation');
grid on;
axis([0, 20, -3, 3]);
