function J = JacobianEst(fun, x, varargin)
%JACOBIANEST  Numerical Jacobian via central finite differences.
%
%   J = JacobianEst(fun, x)
%
%   Inputs:
%       fun - Function handle mapping R^n -> R^m
%       x   - Point at which to evaluate the Jacobian (n x 1)
%
%   Optional name-value pairs:
%       'f0'  - Pre-computed fun(x) to avoid redundant evaluation
%
%   Output:
%       J - Jacobian matrix (m x n), where J(i,j) = d fun_i / d x_j

%% Parse optional f0
f0 = [];
for i = 1:2:numel(varargin)
    if strcmpi(varargin{i}, 'f0')
        f0 = varargin{i+1};
    end
end
if isempty(f0)
    f0 = fun(x);
end

%% Central finite differences
n = numel(x);
m = numel(f0);
J = zeros(m, n);
h = eps^(1/3) * max(abs(x), 1);

for j = 1:n
    xp = x;  xp(j) = xp(j) + h(j);
    xm = x;  xm(j) = xm(j) - h(j);
    J(:, j) = (fun(xp) - fun(xm)) / (2 * h(j));
end
end
