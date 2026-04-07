%% Delta Method — Minimal Example
%  Given sqrt(n)(theta_hat - theta_0) ~ N(0, V), the delta method gives
%  the asymptotic distribution of any smooth transformation g(theta_hat):
%
%       sqrt(n)(g(theta_hat) - g(theta_0)) ~ N(0, G * V * G')
%
%  where G = dg/dtheta' is the Jacobian (Hansen 2021, Thm 6.8).
%  In practice:  varpara = F * varcoef * F'  with F = JacobianEst(g, coef).
%  See THEORY.md for derivation and assumptions.
%
%  When is this useful? Whenever you need the standard error of a
%  quantity that is a function of estimated coefficients:
%    1. Ratio of coefficients (e.g., relative effect of education vs experience)
%    2. Elasticities from a log-linear model (e.g., price elasticity = b_price * p/q)
%    3. Marginal effects in nonlinear models (probit/logit: b * phi(X'b))
%    4. Long-run multiplier in a dynamic model (sum of lag coefficients / (1 - rho))
%    5. Implied structural parameters from reduced-form estimates (IV/SEM)
%    6. Predicted values at a specific point (e.g., wage at education=16, experience=5)
%    7. Differences or sums of coefficients (e.g., testing beta1 = beta2)
ccc

%% Simulated data: OLS regression with two regressors
rng(42);
n = 500;  k = 2;
X = randn(n, k);
y = X * [3; 1] + randn(n, 1);          % true beta = [3; 1]
coef    = (X'*X) \ (X'*y);              % OLS point estimates
varcoef = inv(X'*X) * (sum((y - X*coef).^2) / (n-k));  % OLS covariance

%% Case 1: Linear transformation  R*coef - r
%  Test whether beta1 - beta2 = 0.
%  This is a special case where we know the analytical solution:
%  the Jacobian is simply R, so Var(phi) = R * varcoef * R' exactly.
R = [1 -1];  r = 0;
lin = DeltaMethod(coef, varcoef, R, r);
fprintf('\n--- Linear: beta1 - beta2 ---\n');
disp(lin.FullTable)

%  Verify: analytical variance must match
var_analytical = R * varcoef * R';
fprintf('  Analytical SE:  %.6f\n', sqrt(var_analytical));
fprintf('  Delta method SE: %.6f  (must match exactly)\n\n', lin.stderrpara);

%% Case 2: Nonlinear transformation  g(coef)
%  Estimate the ratio beta1 / beta2.
%  Here the Jacobian G = [1/b2, -b1/b2^2] is evaluated numerically.
g   = @(b) b(1) / b(2);
nln = DeltaMethod(coef, varcoef, g);
fprintf('\n--- Nonlinear: beta1 / beta2 ---\n');
disp(nln.FullTable)

%  Verify: compare numerical Jacobian against analytical Jacobian
G_analytical = [1/coef(2), -coef(1)/coef(2)^2];
se_analytical = sqrt(G_analytical * varcoef * G_analytical');
fprintf('  Analytical SE:  %.6f\n', se_analytical);
fprintf('  Delta method SE: %.6f  (should match to ~10 digits)\n', nln.stderrpara);
