function [Time, Y] = odeRKSSP53(f, TSpan, Y0, h)
    %ODERKSSP53, 3-stage, 3rd Order TVD Runge-Kutta Method of Spiteri-Ruuth Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         5-stage, 3rd order SSP Runge-Kutta Spiteri-Ruuth
    %     Order:
    %                         3
    %     Number of Stages:
    %                         5
    %     Number of Registers:
    %                         2
    %     CFL:
    %                         2.65
    %     Strong Stability Preserving:
    %                         true
    %     Links:
    %         1. https://www.ams.org/journals/mcom/2006-75-253/S0025-5718-05-01772-2/S0025-5718-05-01772-2.pdf
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
    %   [T, Y] = odeRKSSP53(f, TSpan, Y0, h);
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
    %    * Ruuth, Steven. "Global optimization of explicit strong-stability-preserving Runge-Kutta
    %       methods." Mathematics of Computation 75.253 (2006): 183-207.
    %       https://www.ams.org/journals/mcom/2006-75-253/S0025-5718-05-01772-2/S0025-5718-05-01772-2.pdf

    % Set default values if not provided
    if nargin < 4
        h = 0.01;
    end

    t0 = TSpan(1);
    tf = TSpan(2);

    % Initial setup
    t = t0;
    y = Y0;

    % Solver Params
    a30 = 0.355909775063327;
    a32 = 0.644090224936674;
    a40 = 0.367933791638137;
    a43 = 0.632066208361863;
    a52 = 0.237593836598569;
    a54 = 0.762406163401431;
    b10 = 0.377268915331368;
    b21 = 0.377268915331368;
    b32 = 0.242995220537396;
    b43 = 0.238458932846290;
    b54 = 0.287632146308408;
    c1  = 0.377268915331368;
    c2  = 0.754537830662736;
    c3  = 0.728985661612188;
    c4  = 0.699226135931670;

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
        k = f(     t,    y); Tmp1 = y + b10*h*k;
        k = f(t+c1*h, Tmp1); Tmp2 = Tmp1 + b21*h*k;
        k = f(t+c2*h, Tmp2); Tmp1 = a30*y + a32*Tmp2 + b32*h*k;
        k = f(t+c3*h, Tmp1); Tmp1 = a40*y + a43*Tmp1 + b43*h*k;
        k = f(t+c4*h, Tmp1);

        y = a52*Tmp2 + a54*Tmp1 + b54*h*k;
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
