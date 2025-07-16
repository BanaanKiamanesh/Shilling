function [Time, Y] = odeManager(f, TSpan, Y0, varargin)
%ODEMANAGER High-level interface for all ODE solvers in the library.
%
%   [T, Y] = ODEMANAGER(F, TSPAN, Y0) solves the system defined by the
%   function handle F with the default solver (classic RK4) on the time
%   interval TSPAN and with initial condition Y0.
%
%   [T, Y] = ODEMANAGER(F, TSPAN, Y0, OPTIONS) allows the caller to pass
%   an options structure generated with ODEOPTIONS.  Important fields are:
%       OPTIONS.Solver   - Name of the solver (string) or a function handle.
%       OPTIONS.StepSize - Fixed time-step size (scalar) used by fixed-step
%                          solvers.
%       OPTIONS.RelTol   - Relative tolerance (not used by fixed-step solvers).
%       OPTIONS.AbsTol   - Absolute tolerance (not used by fixed-step solvers).
%
%   Inputs:
%     F      - Function handle for the ODE system, f(t, y)
%     TSPAN  - Time interval as [t0, tf]
%     Y0     - Initial condition (column vector)
%     OPTIONS- (optional) struct from odeOptions
%
%   Outputs:
%     T      - Time vector
%     Y      - Solution matrix (each row corresponds to a time in T)
%
%   Example (Van der Pol ODE):
%       f  = @(t, x) [x(2); (1-x(1)^2).*x(2) - x(1)];
%       ts = [0 20];
%       y0 = [2; 0];
%       opt = odeOptions('Solver','rk5','StepSize',1e-2);
%       [T, Y] = odeManager(f, ts, y0, opt);
%       plot(T, Y);
%
%   See also: odeOptions, README.md for more details and examples.

    % Merge varargin into options
    if nargin < 4 || isempty(varargin)
        options = odeOptions();
    elseif isstruct(varargin{1})
        options = odeOptions(varargin{1}, varargin{2:end});
    else
        options = odeOptions(varargin{:});
    end

    % Ensure mandatory option fields exist (fall back to defaults)
    def = odeOptions();
    fields = fieldnames(def);
    for k = 1:numel(fields)
        if ~isfield(options, fields{k})
            options.(fields{k}) = def.(fields{k});
        end
    end

    % Retrieve solver function handle
    solverFunc = getSolverFunction(options.Solver);

    % Determine common parameters
    h = options.StepSize;

    % Call solver (variable number of outputs supported)
    if nargout > 1
        [Time, Y] = solverFunc(f, TSpan, Y0, h);
    else
        Time = solverFunc(f, TSpan, Y0, h);
    end
end

function solverFunc = getSolverFunction(solver)
%GETSOLVERFUNCTION Return function handle corresponding to solver name.
%
%   SOLVER can be either a string (case-insensitive) or a function handle.

    if isa(solver, 'function_handle')
        solverFunc = solver;
        return;
    end

    if ~ischar(solver)
        error('odeManager:InvalidSolver', 'Solver must be a string or function handle.');
    end

    solverKey = lower(strtrim(solver));

    solverMap = struct( ...
        'euler',      @odeEuler, ...
        'midpoint',   @odeMidpoint, ...
        'heun',       @odeHeun, ...
        'rk3',        @odeRK3, ...
        'rk4',        @odeRK4, ...
        'rkr4',       @odeRKR4, ...
        'rks4',       @odeRKS4, ...
        'rkls44',     @odeRKLS44, ...
        'rkls54',     @odeRKLS54, ...
        'rkssp2',     @odeRKSSP2, ...
        'rkssp22',    @odeRKSSP2, ...
        'rkssp33',    @odeRKSSP33, ...
        'rkssp53',    @odeRKSSP53, ...
        'rkssp54',    @odeRKSSP54, ...
        'rks5',       @odeRKS5, ...
        'rk5',        @odeRK5, ...
        'rkc5',       @odeRKC5, ...
        'rkl5',       @odeRKL5, ...
        'rklk5a',     @odeRKLK5a, ...
        'rklk5b',     @odeRKLK5b, ...
        'rkb6',       @odeRKB6, ...
        'rk7',        @odeRK7, ...
        'rk8_10',     @odeRK810, ...
        'rkcv8',      @odeRKCV8, ...
        'rk8_12',     @odeRK812, ...
        'rkz10',      @odeRKZ10, ...
        'rko10',      @odeRKO10, ...
        'rkh10',      @odeRKH10 ...
    );

    if ~isfield(solverMap, solverKey)
        % Offer friendly error with list of available solvers
        names = sort(fieldnames(solverMap));
        error('odeManager:UnknownSolver', ...
            'Unknown solver ''%s''. Available solvers are:\n  * %s', ...
            solverKey, strjoin(names, '\n  * '));
    end

    solverFunc = solverMap.(solverKey);
end