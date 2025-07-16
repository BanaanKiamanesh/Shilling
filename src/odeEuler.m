function [Time, Y] = odeEuler(f, TSpan, Y0, h)
    %ODEEULER Euler's method for solving ordinary differential equations.
    %
    %   [T, Y] = ODEEULER(F, TSPAN, Y0, H) solves the system of differential
    %   equations dy/dt = F(t,y) using Euler's method with fixed step size H.
    %
    %   METHOD PROPERTIES:
    %       Name:           Euler's method
    %       Order:          1
    %       Type:           Fixed step, explicit
    %       Stages:         1
    %       Registers:      1
    %       CFL:            1.0
    %       Year:           1768
    %       Reference:      Euler (1768)
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
    %       % Solve simple exponential decay: dy/dt = -2*y, y(0) = 1
    %       f = @(t, y) -2*y;
    %       tspan = [0, 2];
    %       y0 = 1;
    %       h = 0.1;
    %       [t, y] = odeEuler(f, tspan, y0, h);
    %       plot(t, y, 'o-');
    %       xlabel('Time'); ylabel('y(t)'); title('Exponential Decay');
    %
    %   SEE ALSO:
    %       ODEMANAGER, ODEOPTIONS, ODERK4, ODEHEUN
    %
    %   REFERENCE:
    %       Euler, L. (1768). Institutionum calculi integralis.
    %       https://archive.org/details/institutionescal020326mbp

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
        % Euler's method
        y = y + h * f(t, y);
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
