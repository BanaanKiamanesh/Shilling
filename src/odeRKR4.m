function [Time, Y] = odeRKR4(f, TSpan, Y0, h)
    %ODERKR4 4th Order Runge-Kutta With Minimum Truncation Error Method Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         4th order Runge-Kutta Ralston
    %     Order:
    %                         4
    %     Number of Stages:
    %                         4
    %     Number of Registers:
    %                         4
    %     Links:
    %         1. https://doi.org/10.1090%2FS0025-5718-1962-0150954-0
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
    %   [T, Y] = odeRKR4(f, TSpan, Y0, h);
    %
    %   % Plot results
    %   figure;
    %   plot(T, Y, 'LineWidth', 2);
    %   xlabel('Time [s]');
    %   ylabel('State Variables');
    %   title('Van der Pol Equation');
    %   grid on;
    %   axis([0, 20, -3, 3]);
    %
    % Reference
    %   * Ralston, Anthony (1962).
    %       "[Runge-Kutta Methods with Minimum Error Bounds](https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf)".
    %       Math. Comput. 16 (80): 431-437.

    % Set default values if not provided
    if nargin < 4
        h = 0.01;
    end

    t0 = TSpan(1);
    tf = TSpan(2);

    % Initial setup
    t = t0;
    y = Y0;

    % Method Params
    a2  = 4.0 / 10.0;
    a3  = (14.0 - 3.0 * sqrt(5)) / 16.0;
    b21 = 4.0 / 10.0;
    b31 = (-2889.0 + 1428.0 * sqrt(5)) / 1024.0;
    b32 = (3785.0 - 1620.0 * sqrt(5)) / 1024.0;
    b41 = (-3365.0 + 2094.0 * sqrt(5)) / 6040.0;
    b42 = (-975.0 - 3046.0 * sqrt(5)) / 2552.0;
    b43 = (467040.0 + 203968.0 * sqrt(5)) / 240845.0;
    c1  = (263.0 + 24.0 * sqrt(5)) / 1812.0;
    c2  = (125.0 - 1000.0 * sqrt(5)) / 3828.0;
    c3  = 1024.0 * (3346.0 + 1623.0 * sqrt(5)) / 5924787.0;
    c4  = (30.0 - 4.0 * sqrt(5)) / 123.0;

    % Preallocate arrays to store values
    StepNum = ceil((tf - t0) / h);
    Time = zeros(StepNum, 1);
    Y = zeros(StepNum, length(Y0));
    Time(1) = t;
    Y(1, :) = y';

    % Main loop
    idx = 1;
    while t < tf
        % Method
        k1 = f(       t,                                y);
        k2 = f(t + a2*h,                   y + h*(b21*k1));
        k3 = f(t + a3*h,          y + h*(b31*k1 + b32*k2));
        k4 = f(     t+h, y + h*(b41*k1 + b42*k2 + b43*k3));

        % Update values
        y = y + h*(c1*k1 + c2*k2 + c3*k3 + c4*k4);
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
