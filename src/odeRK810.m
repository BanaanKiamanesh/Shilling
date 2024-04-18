function [Time, Y] = odeRK810(ODEfun, TSpan, Y0, h)
    %ODERK810 8th Order Take-One Runge-Kutta Method Implementation
    % Inputs:
    %   ODEfun: function handle for the ODE
    %   TSpan: time span as [t0, tf]
    %   Y0: initial condition as a column vector (default: zeros)
    %   h: initial step size (default: 0.01)
    %
    % Example Usage:
    %   f = @(t, x) [x(2); (1 - x(1)^2)*x(2) - x(1)]; % Van der Pol ODE
    %   TSpan= [0, 20];    % Time span
    %   Y0= [2; 0];        % Initial condition
    %   h = 0.01;          % Initial step size
    %   [T, Y] = odeRK810(f, TSpan, Y0, h);
    %
    %   % Plot results
    %   figure;
    %   plot(T, Y, 'LineWidth', 2);
    %   xlabel('Time [s]');
    %   ylabel('State Variables');
    %   title('Van der Pol Equation Solution using Runge Kutta 8(10) Method');
    %   grid on;
    %   axis([0, 20, -3, 3]);
    %
    % Reference
    %   1. E. B. Shanks, "[Higher Order Approximations of Runge-Kutta Type](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)",
    %       NASA Technical Note, NASA TN D-2920, Sept. 1965.

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
    StepNum = ceil((tf - t0)  /  h);
    Time    = zeros(StepNum, 1);
    Y       = zeros(StepNum, length(Y0));
    Time(1) = t;
    Y(1, :) = y';

    % Method Params
    a1  =  4.0 / 27.0;
    a2  =  2.0 / 9.0;
    a3  =  1.0 / 3.0;
    a4  =  1.0 / 2.0;
    a5  =  2.0 / 3.0;
    a6  =  1.0 / 6.0;
    a8  =  5.0 / 6.0;

    aa1 =  4.0 / 27.0;
    aa2 =  1.0 / 18.0;
    aa3 =  1.0 / 12.0;
    aa4 =  1.0 / 8.0;
    aa5 =  1.0 / 54.0;
    aa6 =  1.0 / 4320.0;
    aa7 =  1.0 / 20.0;
    aa8 =  1.0 / 288.0;
    aa9 =  1.0 / 820.0;

    b21 =  3.0;
    b32 =  3.0;
    b43 =  3.0;
    b50 =  13.0;
    b52 = -27.0;
    b53 =  42.0;
    b54 =  8.0;
    b60 =  389.0;
    b62 = -54.0;
    b63 =  966.0;
    b64 = -824.0;
    b65 =  243.0;
    b70 = -231.0;
    b72 =  81.0;
    b73 = -1164.0;
    b74 =  656.0;
    b75 = -122.0;
    b76 =  800.0;
    b80 = -127.0;
    b82 =  18.0;
    b83 = -678.0;
    b84 =  456.0;
    b85 = -9.0;
    b86 =  576.0;
    b87 =  4.0;
    b90 =  1481.0;
    b92 = -81.0;
    b93 =  7104.0;
    b94 = -3376.0;
    b95 =  72.0;
    b96 = -5040.0;
    b97 = -60.0;
    b98 =  720.0;

    c   =  1.0 / 840.0;
    c0  =  41.0;
    c3  =  27.0;
    c4  =  272.0;
    c5  =  27.0;
    c6  =  216.0;
    c8  =  216.0;
    c9  =  41.0;

    % Main Loop
    Idx = 1;
    while t < tf
        % Method
        k1  = ODEfun(       t,                                                                                 y);
        k2  = ODEfun(t + a1*h,                                                                      y + aa1*h*k1);
        k3  = ODEfun(t + a2*h,                                                           y + aa2*h*(k1 + b21*k2));
        k4  = ODEfun(t + a3*h,                                                           y + aa3*h*(k1 + b32*k3));
        k5  = ODEfun(t + a4*h,                                                           y + aa4*h*(k1 + b43*k4));
        k6  = ODEfun(t + a5*h,                                     y + aa5*h*(b50*k1 + b52*k3 + b53*k4 + b54*k5));
        k7  = ODEfun(t + a6*h,                            y + aa6*h*(b60*k1 + b62*k3 + b63*k4 + b64*k5 + b65*k6));
        k8  = ODEfun(   t + h,                   y + aa7*h*(b70*k1 + b72*k3 + b73*k4 + b74*k5 + b75*k6 + b76*k7));
        k9  = ODEfun(t + a8*h,          y + aa8*h*(b80*k1 + b82*k3 + b83*k4 + b84*k5 + b85*k6 + b86*k7 + b87*k8));
        k10 = ODEfun(   t + h, y + aa9*h*(b90*k1 + b92*k3 + b93*k4 + b94*k5 + b95*k6 + b96*k7 + b97*k8 + b98*k9));

        % Update Values
        y = y + c * h * (c0*k1 + c3*k4 + c4*k5 + c5*k6 + c6*k7 + c8*k9 + c9*k10);
        t = t + h;

        % Store Values
        Idx       = Idx + 1;
        Time(Idx) = t;
        Y(Idx, :) = y';
    end
end
