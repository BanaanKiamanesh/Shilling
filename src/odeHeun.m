function [Time, Y] = odeHeun(f, TSpan, Y0, h)
    %ODEHEUN Heun's method for solving ordinary differential equations.
    %
    %   [T, Y] = ODEHEUN(F, TSPAN, Y0, H) solves the system of differential
    %   equations dy/dt = F(t,y) using Heun's method (improved Euler)
    %   with fixed step size H.
    %
    %   METHOD PROPERTIES:
    %       Name:           Heun's method (improved Euler)
    %       Order:          2
    %       Type:           Fixed step, explicit
    %       Stages:         2
    %       Registers:      2
    %       CFL:            2.0
    %       Year:           1900
    %       Reference:      Heun (1900)
    %
    %   INPUTS:
    %       F           Function handle @(t,y) that returns dy/dt
    %       TSPAN       Time span [t0, tf] as 2-element vector
    %       Y0          Initial condition as column vector
    %       H           Step size (default: 0.01)
    %
    %   OUTPUTS:
    %       T           Time vector
    %       Y           Solution matrix (each row is a time point)
    %
    %   EXAMPLE:
    %       % Solve logistic equation: dy/dt = r*y*(1-y/K)
    %       r = 0.5; K = 10;
    %       f = @(t, y) r*y*(1-y/K);
    %       tspan = [0, 10];
    %       y0 = 1;
    %       h = 0.1;
    %       [t, y] = odeHeun(f, tspan, y0, h);
    %       plot(t, y, 'o-');
    %       xlabel('Time'); ylabel('Population'); title('Logistic Growth');
    %
    %   SEE ALSO:
    %       ODEMANAGER, ODEOPTIONS, ODEEULER, ODERK4
    %
    %   REFERENCE:
    %       Heun, K. (1900). Neue Methode zur approximativen Integration
    %       der Differentialgleichungen einer unabhängigen Veränderlichen.
    %       Zeitschrift für Mathematik und Physik, 45, 23-38.

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
        k1 = f(    t,        y);
        k2 = f(t + h, y + h*k1);

        y = y + h * (k1 + k2) / 2;
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
