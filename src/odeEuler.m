function [Time, Y] = odeEuler(f, TSpan, Y0, h)
%ODEEULER Euler's ODE Solver Implementation
%
%   Implements the classic explicit Euler method for solving ODEs.
%
%   Method Properties:
%     Method Name:        Euler
%     Introduced in Year: 1768
%     Method Type:        Fixed Time Step, Explicit
%     Order:              1
%     Number of Stages:   1
%     Number of Registers:1
%     CFL:                1.0
%     Reference:          https://archive.org/details/institutionescal020326mbp
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
%     [T, Y] = odeEuler(f, TSpan, Y0, h);
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
        % Euler's method
        y = y + h * f(t, y);
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
