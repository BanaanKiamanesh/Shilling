function [Time, Y] = odeRK4(f, TSpan, Y0, h)
    %ODERK4 Classic 4th-order Runge-Kutta method for solving ODEs.
    %
    %   [T, Y] = ODERK4(F, TSPAN, Y0, H) solves the system of differential
    %   equations dy/dt = F(t,y) using the classic 4th-order Runge-Kutta
    %   method with fixed step size H.
    %
    %   METHOD PROPERTIES:
    %       Name:           Classic 4th-order Runge-Kutta
    %       Order:          4
    %       Type:           Fixed step, explicit
    %       Stages:         4
    %       Registers:      4
    %       CFL:            ~2.8
    %       Year:           1901
    %       Reference:      Kutta (1901)
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
    %       % Solve harmonic oscillator: d²x/dt² + x = 0
    %       % Convert to system: [x; dx/dt] = [x2; -x1]
    %       f = @(t, x) [x(2); -x(1)];
    %       tspan = [0, 2*pi];
    %       y0 = [1; 0];  % x(0) = 1, dx/dt(0) = 0
    %       h = 0.1;
    %       [t, y] = odeRK4(f, tspan, y0, h);
    %       plot(t, y(:,1), 'o-');
    %       xlabel('Time'); ylabel('x(t)'); title('Harmonic Oscillator');
    %
    %   SEE ALSO:
    %       ODEMANAGER, ODEOPTIONS, ODEEULER, ODERK5
    %
    %   REFERENCE:
    %       Kutta, W. (1901). Beitrag zur näherungsweisen Integration
    %       totaler Differentialgleichungen. Zeitschrift für Mathematik
    %       und Physik, 46, 435-453.

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
        % Runge-Kutta method
        k1 = h * f(t, y);
        k2 = h * f(t + h/2, y + k1/2);
        k3 = h * f(t + h/2, y + k2/2);
        k4 = h * f(t + h, y + k3);

        % Update values
        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6;
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
