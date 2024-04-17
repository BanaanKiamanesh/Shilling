function [Time, Y] = odeRKSSP54(f, TSpan, Y0, h)
    %ODERKSSP54 5-Stage Fourth-Order Runge-Kutta Method Implementation (CFL=1.508)
    % Inputs:
    %   f: function handle for the ODE
    %   TSpan: time span as [t0, tf]
    %   Y0: initial condition as a column vector (default: zeros)
    %   h: step size (default: 0.01)
    %
    % Example Usage:
    %   f = @(t, x) [x(2); (1 - x(1)^2)*x(2) - x(1)]; % Van der Pol ODE
    %   TSpan= [0, 20];    % Time span
    %   Y0= [2; 0];        % Initial condition
    %   h = 0.01;          % Step size
    %   [T, Y] = odeRKSSP54(f, TSpan, Y0, h);
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
    %   * Ruuth, Steven. "Global optimization of explicit strong-stability-preserving Runge-Kutta
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

    % Method Params
    b10 = 0.391752226571890;
    a20 = 0.444370493651235;
    a21 = 0.555629506348765;
    b21 = 0.368410593050371;
    a30 = 0.620101851488403;
    a32 = 0.379898148511597;
    b32 = 0.251891774271694;
    a40 = 0.178079954393132;
    a43 = 0.821920045606868;
    b43 = 0.544974750228521;
    a52 = 0.517231671970585;
    a53 = 0.096059710526147;
    b53 = 0.063692468666290;
    a54 = 0.386708617503269;
    b54 = 0.226007483236906;
    c1  = 0.391752226571890;
    c2  = 0.586079689311540;
    c3  = 0.474542363121400;
    c4  = 0.935010630967653;

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
        k1 = f(       t,  y);   Tmp2 = y + b10*h*k1;
        k1 = f(t + c1*h, Tmp2); Tmp2 = a20*y + a21*Tmp2 + b21*h*k1;
        k3 = f(t + c2*h, Tmp2); Tmp3 = a30*y + a32*Tmp2 + b32*h*k1;
                                Tmp4 = a40*y + a43*Tmp3 + b43*h*k3;
        k1 = f(t + c4*h, Tmp4);

        % Update values
        y = a52*Tmp2 + a53*Tmp3 + b53*h*k3 + a54*Tmp4 + b54*h*k1;
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
