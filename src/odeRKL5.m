function [Time, Y] = odeRKL5(f, TSpan, Y0, h)
    %ODERKL5 5th Order Runge-Kutta Lawson Method Implementation
    % Inputs:
    %   f: function handle for the ODE
    %   TSpan: time span as [t0, tf]
    %   Y0: initial condition as a column vector (default: zeros)
    %   h: step size (default: 0.01)
    %
    % Example Usage:
    %   f = @(t, x) [x(2); (1 - x(1)^2)*x(2) - x(1)]; % Van der Pol ODE
    %   TSpan= [0, 20];    % Time span
    %   Y0= [2; 0];        % Initial condition
    %   h = 0.01;          % Step size
    %   [T, Y] = odeRKL5(f, TSpan, Y0, h);
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
    % Reference:
    %   * An Order Five Runge Kutta Process with Extended Region of Stability,
    %       J. Douglas Lawson, Siam Journal on Numerical Analysis,
    %       Vol. 3, No. 4, (Dec., 1966) pages 593-597

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
    a2 = 1.0 / 12.0;
    a3 = 1.0 / 4.0;
    a4 = 1.0 / 2.0;
    a5 = 3.0 / 4.0;

    b21 =  1.0   / 12.0;
    b31 = -1.0   / 8.0;
    b32 =  3.0   / 8.0;
    b41 =  3.0   / 5.0;
    b42 = -9.0   / 10.0;
    b43 =  4.0   / 5.0;
    b51 =  39.0  / 80.0;
    b52 = -9.0   / 20.0;
    b53 =  3.0   / 20.0;
    b54 =  9.0   / 16.0;
    b61 = -59.0  / 35.0;
    b62 =  66.0  / 35.0;
    b63 =  48.0  / 35.0;
    b64 = -12.0  / 7.0;
    b65 =  8.0   / 7.0;

    c1 = 7.0  / 90.0;
    c3 = 16.0 / 45.0;
    c4 = 2.0  / 15.0;
    c5 = 16.0 / 45.0;
    c6 = 7.0  / 90.0;

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
        k1 = f(       t,                                                  y);
        k2 = f(t + a2*h,                                     y + h*(b21*k1));
        k3 = f(t + a3*h,                            y + h*(b31*k1 + b32*k2));
        k4 = f(t + a4*h,                   y + h*(b41*k1 + b42*k2 + b43*k3));
        k5 = f(t + a5*h,          y + h*(b51*k1 + b52*k2 + b53*k3 + b54*k4));
        k6 = f(   t + h, y + h*(b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5));

        % Update values
        y = y + h * (c1*k1 + c3*k3 + c4*k4 + c5*k5 + c6*k6);
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
