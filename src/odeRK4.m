function [Time, Y] = odeRK4(f, TSpan, Y0, h)
    %ODERK4 4th Order Runge-Kutta ODE Solver Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         4th order Runge-Kutta
    %     Order:
    %                         4
    %     Number of Stages:
    %                         4
    %     Number of Registers:
    %                         4
    %     Links:
    %         1. https://archive.org/stream/zeitschriftfrma12runggoog#page/n449/mode/2up
    %
    %
    % Inputs:
    %   f: function handle for the ODE
    %   TSpan: time span as [t0, tf]
    %   Y0: initial condition as a column vector
    %   h: step size (default: 0.01)
    %
    %
    % Outputs:
    %   Time: time vector associated with the integration
    %   Y: solved ode state evolution matrix
    %
    %
    % Example Usage:
    %   f = @(t, x) [x(2); (1 - x(1)^2)*x(2) - x(1)]; % Van der Pol ODE
    %   TSpan= [0, 20];    % Time span
    %   Y0= [2; 0];        % Initial condition
    %   h = 0.01;          % Step size
    %   [T, Y] = odeRK4(f, TSpan, Y0, h);
    %
    %   % Plot results
    %   figure;
    %   plot(T, Y, 'LineWidth', 2);
    %   xlabel('Time [s]');
    %   ylabel('State Variables');
    %   title('Van der Pol Equation');
    %   grid on;
    %   axis([0, 20, -3, 3]);

    % Set default values if not provided
    if nargin < 4
        h = 0.01;
    end

    t0 = TSpan(1);
    tf = TSpan(2);

    % Initial setup
    t = t0;
    y = Y0;

    % Preallocate arrays to store values
    StepNum = ceil((tf - t0) / h);
    Time = zeros(StepNum, 1);
    Y = zeros(StepNum, length(Y0));
    Time(1) = t;
    Y(1, :) = y';

    % Main loop
    idx = 1;
    while t < tf
        % Runge-Kutta method
        k1 = h * f(t, y);
        k2 = h * f(t + h/2, y + k1/2);
        k3 = h * f(t + h/2, y + k2/2);
        k4 = h * f(t + h, y + k3);

        % Update values
        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6;
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
