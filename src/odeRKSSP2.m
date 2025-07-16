function [Time, Y] = odeRKSSP2(f, TSpan, Y0, h)
    %ODERKSSP2 2-stage, 2nd-order TVD Runge-Kutta Shu-Osher method.
    %
    %   [T, Y] = ODERKSSP2(F, TSPAN, Y0, H) solves the system of differential
    %   equations dy/dt = F(t,y) using the 2-stage, 2nd-order TVD Runge-Kutta
    %   Shu-Osher method with fixed step size H.
    %
    %   METHOD PROPERTIES:
    %       Name:           2-stage, 2nd-order TVD Runge-Kutta Shu-Osher
    %       Order:          2
    %       Type:           Fixed step, explicit, SSP
    %       Stages:         2
    %       Registers:      1
    %       CFL:            1.0
    %       SSP:            Strong stability preserving
    %       Year:           1988
    %       Reference:      Shu & Osher (1988)
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
    %       % Solve advection equation with shock: du/dt + u*du/dx = 0
    %       % (demonstrating strong stability preservation)
    %       f = @(t, u) -u.*(circshift(u, -1) - circshift(u, 1)) / (2*0.01);
    %       tspan = [0, 1];
    %       x = linspace(0, 1, 100);
    %       u0 = (x > 0.3 & x < 0.7)';  % Initial step function
    %       h = 0.001;
    %       [t, u] = odeRKSSP2(f, tspan, u0, h);
    %       plot(x, u(end,:), 'o-');
    %       xlabel('x'); ylabel('u(x,t)'); title('Advection with SSP-RK2');
    %
    %   SEE ALSO:
    %       ODEMANAGER, ODEOPTIONS, ODERKSSP33, ODERKSSP53
    %
    %   REFERENCE:
    %       Shu, C.W. & Osher, S. (1988). Efficient implementation of
    %       essentially non-oscillatory shock-capturing schemes.
    %       Journal of Computational Physics, 77(2), 439-471.
    %   axis([0, 20, -3, 3]);
    %
    % References
    %   * C.-W. Shu, S. Osher, "Efficient implementation of essentially non-oscillatory
    %       shock-capturing schemes", Journal of Computational Physics, 77, 1988, 439-471.
    %       https://doi.org/10.1016/0021-9991(88)90177-5.

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
        k = f(    t,   y); Tmp = y + h * k;
        k = f(t + h, Tmp);

        y = (y + Tmp + h*k) / 2;
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
