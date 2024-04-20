function [Time, Y] = odeRK812(ODEfun, TSpan, Y0, h)
    %ODERK812 8th Order Runge-Kutta Shanks Method Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         12-stage, 8th order Runge-Kutta Shanks
    %     Introduced in Year:
    %                         1965
    %     Method Type:
    %                         Fixed Time Step, Explicit
    %     Order:
    %                         8
    %     Number of Stages:
    %                         12
    %     Number of Registers:
    %                         12
    %     Links:
    %         1. http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf
    %
    %
    % Inputs:
    %   ODEfun: function handle for the ODE
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
    %   h = 0.01;          % Initial step size
    %   [T, Y] = odeRK812(f, TSpan, Y0, h);
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
    a1    =  1.0 / 9.0;
    a2    =  1.0 / 6.0;
    a3    =  1.0 / 4.0;
    a4    =  1.0 / 10.0;
    a5    =  1.0 / 6.0;
    a6    =  1.0 / 2.0;
    a7    =  2.0 / 3.0;
    a8    =  1.0 / 3.0;
    a9    =  5.0 / 6.0;
    a10   =  5.0 / 6.0;

    aa1   =  1.0 / 9.0;
    aa2   =  1.0 / 24.0;
    aa3   =  1.0 / 16.0;
    aa4   =  1.0 / 500.0;
    aa5   =  1.0 / 972.0;
    aa6   =  1.0 / 36.0;
    aa7   =  1.0 / 243.0;
    aa8   =  1.0 / 324.0;
    aa9   =  1.0 / 324.0;
    aa10  =  1.0 / 1620.0;
    aa11  =  1.0 / 4428.0;

    b21   =  3.0;
    b32   =  3.0;
    b40   =  29.0;
    b42   =  33.0;
    b43   = -12.0;
    b50   =  33.0;
    b53   =  4.0;
    b54   =  125.0;
    b60   = -21.0;
    b63   =  76.0;
    b64   =  125.0;
    b65   = -162.0;
    b70   = -30.0;
    b73   = -32.0;
    b74   =  125.0;
    b76   =  99.0;
    b80   =  1175.0;
    b83   = -3456.0;
    b84   = -6250.0;
    b85   =  8424.0;
    b86   =  242.0;
    b87   = -27.0;
    b90   =  293.0;
    b93   = -852.0;
    b94   = -1375.0;
    b95   =  1836.0;
    b96   = -118.0;
    b97   =  162.0;
    b98   =  324.0;
    b100  =  1303.0;
    b103  = -4260.0;
    b104  = -6875.0;
    b105  =  9990.0;
    b106  =  1030.0;
    b109  =  162.0;
    b110  = -8595.0;
    b113  =  30720.0;
    b114  =  48750.0;
    b115  = -66096.0;
    b116  =  378.0;
    b117  = -729.0;
    b118  = -1944.0;
    b119  = -1296.0;
    b1110 =  3240.0;

    c     =  1.0 / 840.0;
    c0    =  41.0;
    c5    =  216.0;
    c6    =  272.0;
    c7    =  27.0;
    c8    =  27.0;
    c9    =  36.0;
    c10   =  180.0;
    c11   =  41.0;


    % Preallocate arrays to store values
    StepNum = ceil((tf - t0)  /  h);
    Time    = zeros(StepNum, 1);
    Y       = zeros(StepNum, length(Y0));
    Time(1) = t;
    Y(1, :) = y';

    % Main Loop
    Idx = 1;
    while t < tf
        % Method
        k1   = ODEfun(        t,                                                                                                       y);
        k2   = ODEfun( t + a1*h,                                                                                            y + aa1*h*k1);
        k3   = ODEfun( t + a2*h,                                                                                 y + aa2*h*(k1 + b21*k2));
        k4   = ODEfun( t + a3*h,                                                                                 y + aa3*h*(k1 + b32*k3));
        k5   = ODEfun( t + a4*h,                                                                    y + aa4*h*(b40*k1 + b42*k3 + b43*k4));
        k6   = ODEfun( t + a5*h,                                                                    y + aa5*h*(b50*k1 + b53*k4 + b54*k5));
        k7   = ODEfun( t + a6*h,                                                           y + aa6*h*(b60*k1 + b63*k4 + b64*k5 + b65*k6));
        k8   = ODEfun( t + a7*h,                                                           y + aa7*h*(b70*k1 + b73*k4 + b74*k5 + b76*k7));
        k9   = ODEfun( t + a8*h,                                         y + aa8*h*(b80*k1 + b83*k4 + b84*k5 + b85*k6 + b86*k7 + b87*k8));
        k10  = ODEfun( t + a9*h,                                y + aa9*h*(b90*k1 + b93*k4 + b94*k5 + b95*k6 + b96*k7 + b97*k8 + b98*k9));
        k11  = ODEfun(t + a10*h,                                 y + aa10*h*(b100*k1 + b103*k4 + b104*k5 + b105*k6 + b106*k7 + b109*k10));
        k12  = ODEfun(    t + h, y + aa11*h*(b110*k1 + b113*k4 + b114*k5 + b115*k6 + b116*k7 + b117*k8 + b118*k9 + b119*k10 + b1110*k11));

        % Update Values
        y = y + h * c * (c0*k1 + c5*k6 + c6*k7 + c7*k8 + c8*k9 + c9*k10 + c10*k11 + c11*k12);
        t = t + h;

        % Store Values
        Idx       = Idx + 1;
        Time(Idx) = t;
        Y(Idx, :) = y';
    end
end
