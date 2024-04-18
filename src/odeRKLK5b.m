function [Time, Y] = odeRKLK5b(f, TSpan, Y0, h)
    %ODERKLK5B 5th Order Runge-Kutta Luther and Konen Second Method Implementation
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
    %   [T, Y] = odeRKLK5b(f, TSpan, Y0, h);
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
    %   * H.A.Luther and H.P.Konen,
    %       "Some Fifth-Order Classical Runge Kutta Formulas",
    %       Siam Review, Vol. 3, No. 7, (Oct., 1965) pages 551-558.

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
    a2 = 2.0 / 5.0;
    a3 = 1.0 / 2.0;
    a5 = 1.0 / 2.0 - 1.0 / 10.0 * sqrt(15);
    a6 = 1.0 / 2.0 + 1.0 / 10.0 * sqrt(15);

    b21 =  2.0 / 5.0;
    b31 =  3.0 / 16.0;
    b32 =  5.0 / 16.0;
    b41 =  1.0 / 4.0;
    b42 = -5.0 / 4.0;
    b43 =  2.0;
    b51 =  3.0 / 20.0 - 1.0 / 100.0 * sqrt(15);
    b52 = -1.0 / 4.0;
    b53 =  3.0 / 5.0 - 2.0 / 25.0 * sqrt(15);
    b54 = -1.0 / 100.0 * sqrt(15);
    b61 = -3.0 / 20.0 - 1.0 / 20.0 * sqrt(15);
    b62 = -1.0 / 4.0;
    b63 =  3.0 / 5.0;
    b64 =  3.0 / 10.0 - 1.0 / 20.0 * sqrt(15);
    b65 =  1.0 / 5.0 * sqrt(15);

    c3 = 4.0 / 9.0;
    c5 = 5.0 / 18.0;
    c6 = 5.0 / 18.0;

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
        k4 = f(   t + h,                   y + h*(b41*k1 + b42*k2 + b43*k3));
        k5 = f(t + a5*h,          y + h*(b51*k1 + b52*k2 + b53*k3 + b54*k4));
        k6 = f(t + a6*h, y + h*(b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5));

        % Update values
        y = y + h * (c3*k3 + c5*k5 + c6*k6);
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
