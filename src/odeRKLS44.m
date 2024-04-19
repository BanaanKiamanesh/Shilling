function [Time, Y] = odeRKLS44(f, TSpan, Y0, h)
    %ODERKLS44 4th Order Runge-Kutta Low Storage non-TVD Method Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         4-stage, 4th order low storage non-TVD Runge-Kutta Jiang-Shu
    %     Order:
    %                         4
    %     Number of Stages:
    %                         4
    %     Number of Registers:
    %                         2
    %     Low Storage:
    %                         true
    %     Links:
    %         1. https://ntrs.nasa.gov/api/citations/19960007052/downloads/19960007052.pdf
    %         2. https://iopscience.iop.org/article/10.3847/1538-4365/ab09fb
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
    %   [T, Y] = odeRKLS44(f, TSpan, Y0, h);
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
    %   * Method: Jiang, Guang-Shan, and Chi-Wang Shu. "Efficient implementation of weighted ENO
    %       schemes." Journal of computational physics 126.1 (1996): 202-228.
    %       https://ntrs.nasa.gov/api/citations/19960007052/downloads/19960007052.pdf
    %   * Implementation: J. M. F. Donnert et al 2019 ApJS 241 23.
    %       https://iopscience.iop.org/article/10.3847/1538-4365/ab09fb

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
        Tmp1 = y;           Tmp2 = -4*y/3;
        k = f(      t, Tmp1); Tmp1 = y + h*k/2;   Tmp2 = Tmp2 + Tmp1/3;
        k = f(t + h/2, Tmp1); Tmp1 = y + h*k/2.0; Tmp2 = Tmp2 + 2*Tmp1/3;
        k = f(t + h/2, Tmp1); Tmp1 = y + h*k;     Tmp2 = Tmp2 + Tmp1/3;
        k = f(  t + h, Tmp1); Tmp1 = y + h*k/6;

        % Update values
        y = Tmp1 + Tmp2;
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
