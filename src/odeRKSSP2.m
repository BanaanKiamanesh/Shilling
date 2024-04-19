function [Time, Y] = odeRKSSP2(ODEfun, TSpan, Y0, h)
    %ODERKSSP2, 2-stage, 2nd order TVD Runge-Kutta Shu-Osher ODE Solver Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         2-stage, 2nd Order TVD Runge-Kutta Shu-Osher
    %     Order:
    %                         2
    %     Number of Stages:
    %                         2
    %     Number of Registers:
    %                         1
    %     CFL:
    %                         1.0
    %     Strong Stability Preserving:
    %                         true
    %     Links:
    %         1. https://ntrs.nasa.gov/api/citations/19880014833/downloads/19880014833.pdf
    %
    %
    % Inputs:
    %   ODEfun: function handle for the ODE
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
    %   [T, Y] = odeRKSSP2(f, TSpan, Y0, h);
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
    % References
    %   * C.-W. Shu, S. Osher, "Efficient implementation of essentially non-oscillatory
    %       shock-capturing schemes", Journal of Computational Physics, 77, 1988, 439-471.
    %       https://doi.org/10.1016/0021-9991(88)90177-5.

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
        % Method
        k = ODEfun(    t,   y); Tmp = y + h * k;
        k = ODEfun(t + h, Tmp);

        y = (y + Tmp + h*k) / 2;
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
