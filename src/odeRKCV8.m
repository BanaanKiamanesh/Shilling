function [Time, Y] = odeRKCV8(f, TSpan, Y0, h)
    %ODERKCV8 11-Stage 8th Order Runge-Kutta Cooper-Verner Method Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         11-stage, 8th order Runge-Kutta Cooper-Verner
    %     Introduced in Year:
    %                         1972
    %     Method Type:
    %                         Fixed Time Step, Explicit
    %     Order:
    %                         8
    %     Number of Stages:
    %                         11
    %     Number of Registers:
    %                         11
    %     Links:
    %         1. https://epubs.siam.org/doi/abs/10.1137/0709037
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
    %   h = 0.01;          % Initial step size
    %   [T, Y] = odeRKCV8(f, TSpan, Y0, h);
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
    %   * Some Explicit Runge-Kutta Methods of High Order, by G. J. Cooper and J. H. Verner,
    %       SIAM Journal on Numerical Analysis, Vol. 9, No. 3, (September 1972), pages 389 to 405
    %   * http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK8/RKcoeff8b_1.pdf

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
    a2  = 1.0 / 2.0;
    a3  = 1.0 / 2.0;
    a4  = 1.0 / 2.0 - 1.0 / 14.0*sqrt(21);
    a5  = 1.0 / 2.0 - 1.0 / 14.0*sqrt(21);
    a6  = 1.0 / 2.0;
    a7  = 1.0 / 2.0 + 1.0 / 14.0*sqrt(21);
    a8  = 1.0 / 2.0 + 1.0 / 14.0*sqrt(21);
    a9  = 1.0 / 2.0;
    a10 = 1.0 / 2.0 - 1.0 / 14.0*sqrt(21);

    b21   =  1.0   / 2.0;
    b31   =  1.0   / 4.0;
    b32   =  1.0   / 4.0;
    b41   =  1.0   / 7.0;
    b42   = -1.0   / 14.0   + 3.0  / 98.0  * sqrt(21);
    b43   =  3.0   / 7.0    - 5.0  / 49.0  * sqrt(21);
    b51   =  11.0  / 84.0   - 1.0  / 84.0  * sqrt(21);
    b53   =  2.0   / 7.0    - 4.0  / 63.0  * sqrt(21);
    b54   =  1.0   / 12.0   + 1.0  / 252.0 * sqrt(21);
    b61   =  5.0   / 48.0   - 1.0  / 48.0  * sqrt(21);
    b63   =  1.0   / 4.0    - 1.0  / 36.0  * sqrt(21);
    b64   = -77.0  / 120.0  -7.0   / 180.0 * sqrt(21);
    b65   =  63.0  / 80.0   + 7.0  / 80.0  * sqrt(21);
    b71   =  5.0   / 21.0   + 1.0  / 42.0  * sqrt(21);
    b73   = -48.0  / 35.0   - 92.0 / 315.0 * sqrt(21);
    b74   =  211.0 / 30.0   + 29.0 / 18.0  * sqrt(21);
    b75   = -36.0  / 5.0    - 23.0 / 14.0  * sqrt(21);
    b76   =  9.0   / 5.0    + 13.0 / 35.0  * sqrt(21);
    b81   =  1.0   / 14.0;
    b85   =  1.0   / 9.0    + 1.0  / 42.0  * sqrt(21);
    b86   =  13.0  / 63.0   + 1.0  / 21.0  * sqrt(21);
    b87   =  1.0   / 9.0;
    b91   =  1.0   / 32.0;
    b95   =  91.0  / 576.0  + 7.0  / 192.0 * sqrt(21);
    b96   =  11.0  / 72.0;
    b97   = -385.0 / 1152.0 + 25.0 / 384.0 * sqrt(21);
    b98   =  63.0  / 128.0  - 13.0 / 128.0 * sqrt(21);
    b101  =  1.0   / 14.0;
    b105  =  1.0   / 9.0;
    b106  = -733.0 / 2205.0 + 1.0  / 15.0  * sqrt(21);
    b107  =  515.0 / 504.0  - 37.0 / 168.0 * sqrt(21);
    b108  = -51.0  / 56.0   + 11.0 / 56.0  * sqrt(21);
    b109  =  132.0 / 245.0  - 4.0  / 35.0  * sqrt(21);
    b115  = -7.0   / 3.0    - 7.0  / 18.0  * sqrt(21);
    b116  = -2.0   / 5.0    - 28.0 / 45.0  * sqrt(21);
    b117  = -91.0  / 24.0   + 53.0 / 72.0  * sqrt(21);
    b118  =  301.0 / 72.0   - 53.0 / 72.0  * sqrt(21);
    b119  =  28.0  / 45.0   + 28.0 / 45.0  * sqrt(21);
    b1110 =  49.0  / 18     + 7.0  / 18.0  * sqrt(21);

    c1  = 1.0  / 20.0;
    c8  = 49.0 / 180.0;
    c9  = 16.0 / 45.0;
    c10 = 49.0 / 180.0;
    c11 = 1.0  / 20.0;


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
        k1  = f(    t + h,                                                                   y);
        k2  = f( t + a2*h,                                                      y + h*(b21*k1));
        k3  = f( t + a3*h,                                             y + h*(b31*k1 + b32*k2));
        k4  = f( t + a4*h,                                   y + h*(b41*k1 + b42*k2  + b43*k3));
        k5  = f( t + a5*h,                                   y + h*(b51*k1 + b53*k3  + b54*k4));
        k6  = f( t + a6*h,                         y + h*(b61*k1 + b63*k3  + b64*k4  + b65*k5));
        k7  = f( t + a7*h,                y + h*(b71*k1 + b73*k3  + b74*k4  + b75*k5 + b76*k6));
        k8  = f( t + a8*h,                         y + h*(b81*k1 + b85*k5  + b86*k6  + b87*k7));
        k9  = f( t + a9*h,                y + h*(b91*k1 + b95*k5  + b96*k6  + b97*k7 + b98*k8));
        k10 = f(t + a10*h,   y + h*(b101*k1 + b105*k5 + b106*k6 + b107*k7 + b108*k8 + b109*k9));
        k11 = f(    t + h, y + h*(b115*k5 + b116*k6 + b117*k7 + b118*k8 + b119*k9 + b1110*k10));

        % Update Values
        y = y + h * (c1*k1 + c8*k8 + c9*k9 + c10*k10 + c11*k11);
        t = t + h;

        % Store Values
        Idx       = Idx + 1;
        Time(Idx) = t;
        Y(Idx, :) = y';
    end
end
