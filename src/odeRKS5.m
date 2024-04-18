function [Time, Y] = odeRKS5(f, TSpan, Y0, h)
    %ODERKS5 5th Order Runge-Kutta Shanks Method Implementation
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
    %   [T, Y] = odeRKS5(f, TSpan, Y0, h);
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
    %   * E. B. Shanks, "Solutions of Differential Equations by Evaluations of Functions"
    %       Math. Comp. 20 (1966).

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
    a1  =  1.0 / 9000.0;
    a2  =  3.0 / 10.0;
    a3  =  3.0 / 4.0;
    c   =  1.0 / 1134.0;
    c0  =  105.0;
    c2  =  500.0;
    c3  =  448.0;
    c4  =  81.0;
    aa1 =  1.0 / 9000.0;
    aa2 =  1.0 / 10.0;
    aa3 =  1.0 / 8.0;
    aa4 =  1.0 / 81.0;
    b20 = -4047.0;
    b21 =  4050.0;
    b30 =  20241.0;
    b31 = -20250.0;
    b32 =  15.0;
    b40 = -931041.0;
    b41 =  931500.0;
    b42 = -490.0;
    b43 =  112.0;

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
        k1 = f(       t,                                             y);
        k2 = f(t + a1*h,                                  y + aa1*h*k1);
        k3 = f(t + a2*h,                   y + aa2*h*(b20*k1 + b21*k2));
        k4 = f(t + a3*h,          y + aa3*h*(b30*k1 + b31*k2 + b32*k3));
        k5 = f(       t, y + aa4*h*(b40*k1 + b41*k2 + b42*k3 + b43*k4));

        % Update values
        y = y + h * c * (c0*k1 + c2*k3 + c3*k4 + c4*k5);
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
