function [Time, Y] = odeRKLS54(f, TSpan, Y0, h)
    %ODERKLS54 5-Stage, 4th Order Low Storage Runge-Kutta Method Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         5-stage, 4th Order Low Storage Runge-Kutta Carpenter-Kennedy
    %     Order:
    %                         4
    %     Number of Stages:
    %                         5
    %     Number of Registers:
    %                         2
    %     Low Storage:
    %                         true
    %     CFL:
    %                         0.32
    %     Links:
    %         1. https://ntrs.nasa.gov/api/citations/19940028444/downloads/19940028444.pdf
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
    %   [T, Y] = odeRKLS54(f, TSpan, Y0, h);
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
    %   * Carpenter, Mark H., and Christopher A. Kennedy. Fourth-order 2N-storage Runge-Kutta
    %      schemes. No. NASA-TM-109112. 1994.
    %      https://ntrs.nasa.gov/api/citations/19940028444/downloads/19940028444.pdf

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
    a = zeros(1, 5);
    b = a;
    c = a;

    a(2:5) = [ -567301805773.0  / 1357537059087.0, ...
        -2404267990393.0 / 2016746695238.0, ...
        -3550918686646.0 / 2091501179385.0, ...
        -1275806237668.0 / 842570457699.0];

    b(1:5) = [1432997174477.0 / 9575080441755.0,  ...
        5161836677717.0 / 13612068292357.0, ...
        1720146321549.0 / 2090206949498.0,  ...
        3134564353537.0 / 4481467310338.0,  ...
        2277821191437.0 / 14882151754819.0];

    c(2:5) = [1432997174477.0 / 9575080441755.0, ...
        2526269341429.0 / 6820363962896.0, ...
        2006345519317.0 / 3224310063776.0, ...
        2802321613138.0 / 2924317926251.0];

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
        k   = f(t, y);
        ds  = h*k;
        Tmp = y + b(1)*ds;

        for i = 2:5
            k   = f(t + c(i)*h, Tmp);
            ds  = a(i)*ds + h*k;
            Tmp = Tmp + b(i)*ds;
        end

        % Update values
        y = Tmp;
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
