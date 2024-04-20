function [Time, Y] = odeRKS4(f, TSpan, Y0, h)
    %ODERKS4 Fourth-Order Runge-Kutta Shanks Method Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         4th order Runge-Kutta Shanks
    %     Introduced in Year:
    %                         1965
    %     Method Type:
    %                         Fixed Time Step, Explicit
    %     Order:
    %                         4
    %     Number of Stages:
    %                         4
    %     Number of Registers:
    %                         4
    %     Links:
    %         http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf
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
    %   [T, Y] = odeRKS4(f, TSpan, Y0, h);
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
    %    * E. B. Shanks, "Solutions of Differential Equations by Evaluations of Functions"
    %        Math. Comp. 20 (1966).

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
    a1  =  1.0 / 100.0;
    a2  =  3.0 / 5.0;
    c   =  1.0 / 70092.0;
    c0  = -179124.0;
    c1  =  200000.0;
    c2  =  40425.0;
    c3  =  8791.0;
    aa1 =  1.0 / 100.0;
    aa2 =  1.0 / 245.0;
    aa3 =  1.0 / 8791.0;
    b20 = -4278.0;
    b21 =  4425.0;
    b30 =  524746.0;
    b31 = -532125.0;
    b32 =  16170.0;

    % Preallocate arrays to store values
    StepNum = ceil((tf - t0) / h);
    Time = zeros(StepNum, 1);
    Y = zeros(StepNum, length(Y0));
    Time(1) = t;
    Y(1, :) = y';

    % Main loop
    idx = 1;
    while t < tf
        % Runge-Kutta method
        f0 = f(       t,                                    y);
        f1 = f(t + a1*h,                         y + aa1*h*f0);
        f2 = f(t + a2*h,          y + aa2*h*(b20*f0 + b21*f1));
        f3 = f(     t+h, y + aa3*h*(b30*f0 + b31*f1 + b32*f2));

        % Update values
        y = y + h*c*(c0*f0 + c1*f1 + c2*f2 + c3*f3);
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
