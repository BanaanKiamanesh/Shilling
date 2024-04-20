function [Time, Y] = odeRKB6(f, TSpan, Y0, h)
    %ODERKB6 6th Order Runge-Kutta Butcher's Method Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         6th order Runge-Kutta Butcher
    %     Introduced in Year:
    %                         1963
    %     Method Type:
    %                         Fixed Time Step, Explicit
    %     Order:
    %                         6
    %     Number of Stages:
    %                         7
    %     Number of Registers:
    %                         7
    %     Links:
    %         1. https://www.cambridge.org/core/services/aop-cambridge-core/content/view/40DFE501CAB781C9AAE1439B6B8F481A/S1446788700023387a.pdf/div-class-title-on-runge-kutta-processes-of-high-order-div.pdf
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
    %   [T, Y] = odeRKB6(f, TSpan, Y0, h);
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
    %   * Butcher, J. (1964). On Runge-Kutta processes of high order.
    %       Journal of the Australian Mathematical Society, 4(2), 179-194.

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
    a2 = 1.0 / 3.0;
    a3 = 2.0 / 3.0;
    a4 = 1.0 / 3.0;
    a5 = 1.0 / 2.0;
    a6 = 1.0 / 2.0;

    b21 =  1.0  / 3.0;
    b32 =  2.0  / 3.0;
    b41 =  1.0  / 12.0;
    b42 =  1.0  / 3.0;
    b43 = -1.0  / 12.0;
    b51 = -1.0  / 16.0;
    b52 =  9.0  / 8.0;
    b53 = -3.0  / 16.0;
    b54 = -3.0  / 8.0;
    b62 =  9.0  / 8.0;
    b63 = -3.0  / 8.0;
    b64 = -3.0  / 4.0;
    b65 =  1.0  / 2.0;
    b71 =  9.0  / 44.0;
    b72 = -9.0  / 11.0;
    b73 =  63.0 / 44.0;
    b74 =  18.0 / 11.0;
    b76 = -16.0 / 11.0;

    c1 =  11.0 / 120.0;
    c3 =  27.0 / 40.0;
    c4 =  27.0 / 40.0;
    c5 = -4.0  / 15.0;
    c6 = -4.0  / 15.0;
    c7 =  11.0 / 120.0;

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
        k3 = f(t + a3*h,                                     y + h*(b32*k2));
        k4 = f(t + a4*h,                   y + h*(b41*k1 + b42*k2 + b43*k3));
        k5 = f(t + a5*h,          y + h*(b51*k1 + b52*k2 + b53*k3 + b54*k4));
        k6 = f(t + a6*h,          y + h*(b62*k2 + b63*k3 + b64*k4 + b65*k5));
        k7 = f(   t + h, y + h*(b71*k1 + b72*k2 + b73*k3 + b74*k4 + b76*k6));

        % Update values
        y = y + h * (c1*k1 + c3*k3 + c4*k4 + c5*k5 + c6*k6 + c7*k7);
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
