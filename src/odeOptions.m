function options = odeOptions(varargin)
%ODEOPTIONS Create/modify option structure for ODE solvers.
%
%   OPTIONS = ODEOPTIONS() creates an options structure OPTIONS where
%   each field is set to the default value.
%
%   OPTIONS = ODEOPTIONS('NAME1',VALUE1,'NAME2',VALUE2,...) creates an
%   options structure OPTIONS in which the named properties have the
%   specified values.  Any unspecified properties have default values.
%
%   OPTIONS = ODEOPTIONS(OLDOPTIONS,'NAME',VALUE,...) modifies an existing
%   options structure OLDOPTIONS by changing the specified properties and
%   returns the modified structure OPTIONS.
%
%   Recognized property names:
%       'Solver'    - Solver to use (string or function handle) [default: 'rk4']
%       'StepSize'  - Fixed time-step size (scalar)             [default: 0.01]
%       'RelTol'    - Relative tolerance (scalar)               [default: 1e-6]
%       'AbsTol'    - Absolute tolerance (scalar)               [default: 1e-9]
%       'MaxStep'   - Maximum step size (scalar)                []
%       'UserData'  - Arbitrary user data                       []
%
%   Example:
%       opts = odeOptions('Solver','rk4','StepSize',1e-2,'RelTol',1e-8);
%
%   See also ODEMANAGER.

    % Default option values
    defaultOptions = struct( ...
        'Solver',   'rk4', ...
        'StepSize', 0.01, ...
        'RelTol',   1e-6, ...
        'AbsTol',   1e-9, ...
        'MaxStep',  [],   ...
        'UserData', []    ...
    );

    % Handle no-input case
    if nargin == 0
        options = defaultOptions;
        return;
    end

    argIdx = 1;

    % If the first argument is an existing options structure, start from it
    if isstruct(varargin{1})
        options = varargin{1};
        argIdx = 2;
    else
        options = defaultOptions;
    end

    % The rest of the arguments must come in name/value pairs
    if mod(nargin - argIdx + 1, 2) ~= 0
        error('odeOptions:NameValuePairs', ...
              'Properties and values must be supplied in pairs.');
    end

    % Parse name/value pairs (case-insensitive)
    for k = argIdx:2:nargin
        name  = varargin{k};
        value = varargin{k+1};

        if ~ischar(name)
            error('odeOptions:InvalidPropertyName', ...
                  'Property names must be character vectors.');
        end

        % Allow case-insensitive property names
        switch lower(strtrim(name))
            case 'solver'
                options.Solver = value;
            case 'stepsize'
                options.StepSize = value;
            case 'reltol'
                options.RelTol = value;
            case 'abstol'
                options.AbsTol = value;
            case 'maxstep'
                options.MaxStep = value;
            case 'userdata'
                options.UserData = value;
            otherwise
                warning('odeOptions:UnknownProperty', ...
                        'Unknown property ''%s'' ignored.', name);
        end
    end
end