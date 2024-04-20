function [Time, Y] = odeRK5(f, TSpan, Y0, h)
    %ODERK5 5th Order Runge-Kutta Method Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         5th order Runge-Kutta
    %     Introduced in Year:
    %                         1965
    %     Method Type:
    %                         Fixed Time Step, Explicit
    %     Order:
    %                         5
    %     Number of Stages:
    %                         6
    %     Number of Registers:
    %                         6
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
    %   [T, Y] = odeRK5(f, TSpan, Y0, h);
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

    % Method Params
    a2 = 1.0 / 5.0;
    a3 = 2.0 / 5.0;
    a5 = 3.0 / 5.0;
    a6 = 4.0 / 5.0;

    b21 = 1.0 / 5.0;
    b32 = 2.0 / 5.0;
    b41 = 9.0 / 4.0;
    b42 = -5.0;
    b43 = 15.0 / 4.0;
    b51 = -63.0 / 100.0;
    b52 = 9.0 / 5.0;
    b53 = -13.0 / 20.0;
    b54 = 2.0 / 25.0;
    b61 = -6.0 / 25.0;
    b62 = 4.0 / 5.0;
    b63 = 2.0 / 15.0;
    b64 = 8.0 / 75.0;

    c1 = 17.0 / 144.0;
    c3 = 25.0 / 36.0;
    c4 = 1.0 / 72.0;
    c5 = -25.0 / 72.0;
    c6 = 25.0 / 48.0;


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
        k1 = f(       t,                                         y);
        k2 = f(t + a2*h,                            y + h*(b21*k1));
        k3 = f(t + a3*h,                            y + h*(b32*k2));
        k4 = f(   t + h,          y + h*(b41*k1 + b42*k2 + b43*k3));
        k5 = f(t + a5*h, y + h*(b51*k1 + b52*k2 + b53*k3 + b54*k4));
        k6 = f(t + a6*h, y + h*(b61*k1 + b62*k2 + b63*k3 + b64*k4));

        % Update values
        y = y + h * (c1*k1 + c3*k3 + c4*k4 + c5*k5 + c6*k6);
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
