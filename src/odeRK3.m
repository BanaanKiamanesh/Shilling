function [Time, Y] = odeRK3(f, TSpan, Y0, h)
    %ODERK3, 3rd Order, 3-Stages Runge-Kutta ODE Solver Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         3th order Runge-Kutta
    %     Method Type:
    %                         Fixed Time Step, Explicit
    %     Order:
    %                         3
    %     Number of Stages:
    %                         3
    %     Number of Registers:
    %                         3
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
    %   [T, Y] = odeRK3(f, TSpan, Y0, h);
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

    % RK3 Params
    a1  = 1/6;
    a2  = 2/3;
    a3  = 1/6;
    b2  = 1/2;
    c21 = 1/2;
    c32 = 2;

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
        k1 = f(       t,                   y);
        k2 = f(t + b2*h,        y + h*c21*k1);
        k3 = f(   t + h, y + h*(-k1 + c32*k2));

        y = y + h * (a1*k1 + a2*k2 + a3*k3);
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
