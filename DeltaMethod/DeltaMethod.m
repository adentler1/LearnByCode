function strc = DeltaMethod(varargin)
%% DELTA METHOD
%  Variance estimation for transformations of estimated parameters.
%
%  Given an estimator theta_hat with asymptotic distribution
%       sqrt(n)(theta_hat - theta_0) --> N(0, V)
%  and a continuously differentiable function g: R^k -> R^q, the delta
%  method yields (Hansen 2021, Thm 6.8; Robinson 2008, Lemma 1):
%       sqrt(n)(g(theta_hat) - g(theta_0)) --> N(0, G * V * G')
%  where G = dg/dtheta' is the q x k Jacobian evaluated at theta_0.
%
%  In finite-sample practice, we replace theta_0 with theta_hat:
%       para    = g(coef)
%       varpara = F * varcoef * F'
%  where F = JacobianEst(g, coef) is computed numerically.
%
%  Usage:
%   strc = DeltaMethod(coef, varcoef, g)       % nonlinear g(coef)
%   strc = DeltaMethod(coef, varcoef, R, r)    % linear R*coef - r
%
%  Inputs:
%   coef    - k x 1 vector of estimated coefficients
%   varcoef - k x k estimated covariance matrix
%   g       - function handle g: R^k -> R^q (nonlinear case)
%   R, r    - q x k matrix and q x 1 vector (linear case: phi = R*coef - r)
%
%  Optional name-value pairs (after positional args):
%   'lbl'            - string array of output parameter labels
%   'lblstem'        - stem for auto-generated labels (default: "delta")
%   'CriticalValues' - significance levels (default: [.9 .95 .99])
%   'cellfrmt'       - number format (default: '%- 12.4f')
%   'DerivArg'       - cell array of extra arguments passed to JacobianEst
%
%  Outputs (struct):
%   strc.para       - q x 1 transformed point estimates
%   strc.varpara    - q x q covariance matrix of transformed estimates
%   strc.stderrpara - q x 1 standard errors
%   strc.FullTable  - table with estimates, SEs, t-stats, p-values, stars
%   strc.SlimTable  - compact table with formatted estimate strings
%
%  References:
%   Hansen, B. E. (2021). Econometrics. Princeton. Theorem 6.8, Section 7.10.
%   Robinson, P. M. (2008). EC484 Lecture Notes, LSE. Section 5.2.1.
%
%  See also: JacobianEst, CreateCoefficientTable, MinimalExample

%% SET OPTIONS
CriticalValues = [.9 .95 .99];
lbl = [];
lblstem=[];
cellfrmt = '%- 12.4f';
DerivArg = {};

%% DYNAMIC READ IN
% Handle backward compatibility with older versions or different argument patterns
if nargin < 1
    help DeltaMethod
    return
elseif nargin > 1 && isa(varargin{1}, 'function_handle') && isa(varargin{2}, 'struct')
    warning("The 'NonLinearDeltaMethod' approach is deprecated and will be removed in future releases.")
    strc = NonLinearDeltaMethod(varargin{:});
    return
elseif nargin > 1 && isa(varargin{2}, 'function_handle') && isa(varargin{1}, 'struct')
    % OldFunctionalDeltaMethod format
    warning("The 'OldFunctionalDeltaMethod' approach is deprecated and will be removed in future releases.")
    strc = OldFunctionalDeltaMethod(varargin{:});
    return
elseif nargin == 3 && isnumeric(varargin{1}) && isnumeric(varargin{2}) && isa(varargin{3}, 'struct')
    warning("The 'LinearDeltaMethod' approach is deprecated and will be removed in future releases.")
    strc = LinearDeltaMethod(varargin{:});
    return
elseif nargin > 2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
    % Read in coefficient data and arguments
    coef = varargin{1};
    varcoef = varargin{2};
    arg1 = varargin{3};
    %   Determine how many positional arguments were consumed
    if isa(arg1, 'function_handle')
        % Nonlinear case: DeltaMethod(coef, varcoef, g, ...)
        arg2 = [];
        nPosArgs = 3;
    elseif nargin > 3 && isnumeric(varargin{4})
        % Linear case: DeltaMethod(coef, varcoef, R, r, ...)
        arg2 = varargin{4};
        nPosArgs = 4;
    else
        arg2 = [];
        nPosArgs = 3;
    end
    %   Shift `varargin` to exclude the positional arguments
    varargin = varargin(nPosArgs+1:end);
else
    error("Invalid input format.")
end

% Parse input arguments
[opt, argList] = parseInputs(varargin, who);
% Set options from the structure
setOptionsFromStruct(opt, who);
setOptionsFromNameValuePairs(argList, who);

%% READ IN
% Two use cases:
% 1. A linear transformation using a matrix "R" and vector "r" (from `arg1` and `arg2`).
% 2. A nonlinear function `fun` that only uses `arg1`.
if isa(arg1, 'function_handle') && isempty(arg2)
    fun = arg1;
elseif ~isempty(arg1) && ~isempty(arg2) && isnumeric(arg1) && isnumeric(arg2)
    fun = @(coef) arg1 * coef - arg2;
else
    error("Check input formats.")
end

%% CORE CALCULATIONS
% Convert point estimates
para = fun(coef);
% Compute the derivative (Jacobian)
F = JacobianEst(fun, coef, DerivArg{:});
% Convert dispersion estimates
varpara = F * varcoef * F';
stderrpara = sqrt(diag(varpara));

%% CREATE TABLE
% Check for variable names
if isempty(lblstem)
    lblstem="delta";
end
if isempty(lbl)
    lbl = lblstem + string((1:numel(para))');
end

% Create tables (formatted for output)
[SlimTable, FullTable] = CreateCoefficientTable(para, stderrpara, ...
    'VarNames', lbl, 'cellfrmt', cellfrmt, 'CriticalValues', CriticalValues);

%% COLLECT OUTPUT
strc.lbl = lbl;
strc.para = para(:);
strc.varpara = varpara;
strc.stderrpara = stderrpara;
strc.FullTable = FullTable;
strc.SlimTable = SlimTable;
end
