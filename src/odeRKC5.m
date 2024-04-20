function [Time, Y] = odeRKC5(f, TSpan, Y0, h)
    %ODERKC5 Cassity's 5th Order Runge-Kutta Method Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         5th order Runge-Kutta Cassity
    %     Introduced in Year:
    %                         1966
    %     Method Type:
    %                         Fixed Time Step, Explicit
    %     Order:
    %                         5
    %     Number of Stages:
    %                         6
    %     Number of Registers:
    %                         6
    %     Links:
    %         1. https://epubs.siam.org/doi/10.1137/0703052
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
    %   [T, Y] = odeRKC5(f, TSpan, Y0, h);
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
    %   * C.R. Cassity, Solutions of the fifth order Runge-Kutta equations,
    %       SIAM J. Numer. Anal., 3, (1966), pp. 598-606
    %   * [Coefficients](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK5/RKcoeff5a_3.pdf)

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
    a2 = 1.0 / 7.0;
    a3 = 5.0 / 14.0;
    a4 = 9.0 / 14.0;
    a5 = 6.0 / 7.0;

    b21 =  1.0      / 7.0;
    b31 = -367.0    / 4088.0;
    b32 =  261.0    / 584.0;
    b41 =  41991.0  / 2044.0;
    b42 = -2493.0   / 73.0;
    b43 =  57.0     / 4.0;
    b51 = -108413.0 / 196224.0;
    b52 =  58865.0  / 65408.0;
    b53 =  5.0      / 16.0;
    b54 =  265.0    / 1344.0;
    b61 = -204419.0 / 58984.0;
    b62 =  143829.0 / 58984.0;
    b63 =  171.0    / 202.0;
    b64 =  2205.0   / 404.0;
    b65 = -432.0    / 101.0;

    c1 =  1.0   / 9.0;
    c2 =  7.0   / 2700.0;
    c3 =  413.0 / 810.0;
    c4 =  7.0   / 450.0;
    c5 =  28.0  / 75.0;
    c6 = -101.0 / 8100.0;


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
