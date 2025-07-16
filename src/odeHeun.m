function [Time, Y] = odeHeun(f, TSpan, Y0, h)
%ODEHEUN Heun's ODE Solver Implementation
%
%   Implements the explicit Heun (improved Euler) method for solving ODEs.
%
%   Method Properties:
%     Method Name:        Heun
%     Method Type:        Fixed Time Step, Explicit
%     Order:              2
%     Number of Stages:   2
%     Number of Registers:2
%
%   Inputs:
%     f     - Function handle for the ODE, f(t, y)
%     TSpan - Time span as [t0, tf]
%     Y0    - Initial condition (column vector)
%     h     - Step size (default: 0.01)
%
%   Outputs:
%     Time  - Time vector
%     Y     - Solution matrix (each row corresponds to a time in Time)
%
%   Example Usage:
%     f = @(t, x) [x(2); (1 - x(1)^2)*x(2) - x(1)];
%     TSpan = [0, 20];
%     Y0 = [2; 0];
%     h = 0.01;
%     [T, Y] = odeHeun(f, TSpan, Y0, h);
%     plot(T, Y);
%
%   See also: odeManager, odeOptions, README.md

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
