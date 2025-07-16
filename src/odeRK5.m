function [Time, Y] = odeRK5(f, TSpan, Y0, h)
    %ODERK5 5th-order Runge-Kutta method for solving ODEs.
    %
    %   [T, Y] = ODERK5(F, TSPAN, Y0, H) solves the system of differential
    %   equations dy/dt = F(t,y) using a 5th-order Runge-Kutta method
    %   with fixed step size H.
    %
    %   METHOD PROPERTIES:
    %       Name:           5th-order Runge-Kutta
    %       Order:          5
    %       Type:           Fixed step, explicit
    %       Stages:         6
    %       Registers:      6
    %       CFL:            ~3.2
    %       Year:           1965
    %       Reference:      Classical method
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
    %       % Solve pendulum equation: d²θ/dt² + (g/L)*sin(θ) = 0
    %       % Convert to system: [θ; dθ/dt] = [θ2; -(g/L)*sin(θ1)]
    %       g = 9.81; L = 1;
    %       f = @(t, x) [x(2); -(g/L)*sin(x(1))];
    %       tspan = [0, 10];
    %       y0 = [pi/4; 0];  % θ(0) = π/4, dθ/dt(0) = 0
    %       h = 0.01;
    %       [t, y] = odeRK5(f, tspan, y0, h);
    %       plot(t, y(:,1), 'o-');
    %       xlabel('Time'); ylabel('θ(t)'); title('Pendulum Oscillation');
    %
    %   SEE ALSO:
    %       ODEMANAGER, ODEOPTIONS, ODERK4, ODERK7
    %
    %   REFERENCE:
    %       Classical 5th-order Runge-Kutta method.

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
    a2 = 1.0 / 5.0;
    a3 = 2.0 / 5.0;
    a5 = 3.0 / 5.0;
    a6 = 4.0 / 5.0;

    b21 = 1.0 / 5.0;
    b32 = 2.0 / 5.0;
    b41 = 9.0 / 4.0;
    b42 = -5.0;
    b43 = 15.0 / 4.0;
    b51 = -63.0 / 100.0;
    b52 = 9.0 / 5.0;
    b53 = -13.0 / 20.0;
    b54 = 2.0 / 25.0;
    b61 = -6.0 / 25.0;
    b62 = 4.0 / 5.0;
    b63 = 2.0 / 15.0;
    b64 = 8.0 / 75.0;

    c1 = 17.0 / 144.0;
    c3 = 25.0 / 36.0;
    c4 = 1.0 / 72.0;
    c5 = -25.0 / 72.0;
    c6 = 25.0 / 48.0;


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
        k1 = f(       t,                                         y);
        k2 = f(t + a2*h,                            y + h*(b21*k1));
        k3 = f(t + a3*h,                            y + h*(b32*k2));
        k4 = f(   t + h,          y + h*(b41*k1 + b42*k2 + b43*k3));
        k5 = f(t + a5*h, y + h*(b51*k1 + b52*k2 + b53*k3 + b54*k4));
        k6 = f(t + a6*h, y + h*(b61*k1 + b62*k2 + b63*k3 + b64*k4));

        % Update values
        y = y + h * (c1*k1 + c3*k3 + c4*k4 + c5*k5 + c6*k6);
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
