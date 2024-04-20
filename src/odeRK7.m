function [Time, Y] = odeRK7(f, TSpan, Y0, h)
    %ODERK7 7th Order Take-One Runge-Kutta Method Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         7th Order Runge-Kutta Shanks
    %     Introduced in Year:
    %                         1965
    %     Method Type:
    %                         Fixed Time Step, Explicit
    %     Order:
    %                         7
    %     Number of Stages:
    %                         9
    %     Number of Registers:
    %                         9
    %     Links:
    %         1. http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf
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
    %   [T, Y] = odeRK7(f, TSpan, Y0, h);
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
    % Reference:
    %   * E. B. Shanks, "Solutions of Differential Equations by Evaluations of Functions"
    %       Math. Comp. 20 (1966).

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
    a1 = 2.0 / 9.0;
    a2 = 1.0 / 3.0;
    a3 = 1.0 / 2.0;
    a4 = 1.0 / 6.0;
    a5 = 8.0 / 9.0;
    a6 = 1.0 / 9.0;
    a7 = 5.0 / 6.0;

    aa1 = 2.0 / 9.0;
    aa2 = 1.0 / 12.0;
    aa3 = 1.0 / 8.0;
    aa4 = 1.0 / 216.0;
    aa5 = 1.0 / 729.0;
    aa6 = 1.0 / 151632.0;
    aa7 = 1.0 / 1375920.0;
    aa8 = 1.0 / 251888.0;

    b21 =  3.0;
    b32 =  3.0;
    b40 =  23.0;
    b42 =  21.0;
    b43 = -8.0;
    b50 = -4136.0;
    b52 = -13584.0;
    b53 =  5264.0;
    b54 =  13104.0;
    b60 =  105131.0;
    b62 =  302016.0;
    b63 = -107744.0;
    b64 = -284256.0;
    b65 =  1701.0;
    b70 = -775229.0;
    b72 = -2770950.0;
    b73 =  1735136.0;
    b74 =  2547216.0;
    b75 =  81891.0;
    b76 =  328536.0;
    b80 =  23569.0;
    b82 = -122304.0;
    b83 = -20384.0;
    b84 =  695520.0;
    b85 = -99873.0;
    b86 = -466560.0;
    b87 =  241920.0;

    c = 1.0 / 2140320.0;

    c0 =  110201.0;
    c3 =  767936.0;
    c4 =  635040.0;
    c5 = -59049.0;
    c6 = -59049.0;
    c7 =  635040.0;
    c8 =  110201.0;

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
        k1 = f(       t,                                                                        y);
        k2 = f(t + a1*h,                                                             y + aa1*h*k1);
        k3 = f(t + a2*h,                                                  y + aa2*h*(k1 + b21*k2));
        k4 = f(t + a3*h,                                                  y + aa3*h*(k1 + b32*k3));
        k5 = f(t + a4*h,                                     y + aa4*h*(b40*k1 + b42*k3 + b43*k4));
        k6 = f(t + a5*h,                            y + aa5*h*(b50*k1 + b52*k3 + b53*k4 + b54*k5));
        k7 = f(t + a6*h,                   y + aa6*h*(b60*k1 + b62*k3 + b63*k4 + b64*k5 + b65*k6));
        k8 = f(t + a7*h,          y + aa7*h*(b70*k1 + b72*k3 + b73*k4 + b74*k5 + b75*k6 + b76*k7));
        k9 = f(   t + h, y + aa8*h*(b80*k1 + b82*k3 + b83*k4 + b84*k5 + b85*k6 + b86*k7 + b87*k8));

        % Update values
        y = y + h * c * (c0*k1 + c3*k4 + c4*k5 + c5*k6 + c6*k7 + c7*k8 + c8*k9);
        t = t + h;

        % Store values
        idx = idx + 1;
        Time(idx) = t;
        Y(idx, :) = y';
    end
end
