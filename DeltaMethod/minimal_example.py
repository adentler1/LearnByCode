"""
Delta Method -- Minimal Example (Python)

Given sqrt(n)(theta_hat - theta_0) ~ N(0, V), the delta method gives
the asymptotic distribution of any smooth transformation g(theta_hat):

    sqrt(n)(g(theta_hat) - g(theta_0)) ~ N(0, G @ V @ G')

where G = dg/dtheta' is the Jacobian (Hansen 2021, Thm 6.8).
In practice:  varpara = F @ varcoef @ F.T  with F = jacobian_est(g, coef).
See THEORY.md for derivation and assumptions.

When is this useful? Whenever you need the standard error of a
quantity that is a function of estimated coefficients:
  1. Ratio of coefficients (e.g., relative effect of education vs experience)
  2. Elasticities from a log-linear model (e.g., price elasticity = b_price * p/q)
  3. Marginal effects in nonlinear models (probit/logit: b * phi(X'b))
  4. Long-run multiplier in a dynamic model (sum of lag coefficients / (1 - rho))
  5. Implied structural parameters from reduced-form estimates (IV/SEM)
  6. Predicted values at a specific point (e.g., wage at education=16, experience=5)
  7. Differences or sums of coefficients (e.g., testing beta1 = beta2)
"""

import numpy as np
from scipy import stats


# ── Core functions ────────────────────────────────────────────────

def jacobian_est(fun, x):
    """Numerical Jacobian via central finite differences.

    Parameters
    ----------
    fun : callable, R^k -> R^q
    x   : (k,) array, point at which to evaluate

    Returns
    -------
    J : (q, k) Jacobian matrix
    """
    x = np.asarray(x, dtype=float)
    f0 = np.atleast_1d(fun(x))
    n = len(x)
    m = len(f0)
    J = np.zeros((m, n))
    h = np.finfo(float).eps ** (1/3) * np.maximum(np.abs(x), 1.0)
    for j in range(n):
        xp = x.copy(); xp[j] += h[j]
        xm = x.copy(); xm[j] -= h[j]
        J[:, j] = (fun(xp) - fun(xm)) / (2 * h[j])
    return J


def delta_method(coef, varcoef, g):
    """Delta method for a nonlinear transformation g(coef).

    Parameters
    ----------
    coef    : (k,) estimated coefficients
    varcoef : (k, k) estimated covariance matrix
    g       : callable, R^k -> R^q

    Returns
    -------
    dict with keys: para, varpara, se, t_stat, p_value
    """
    coef = np.asarray(coef, dtype=float)
    para = np.atleast_1d(g(coef))
    F = jacobian_est(g, coef)
    varpara = F @ varcoef @ F.T
    se = np.sqrt(np.diag(varpara))
    t_stat = para / se
    p_value = 2 * stats.norm.sf(np.abs(t_stat))
    return dict(para=para, varpara=varpara, se=se,
                t_stat=t_stat, p_value=p_value)


def delta_method_linear(coef, varcoef, R, r):
    """Delta method for a linear transformation R @ coef - r.

    This is a special case where the Jacobian is simply R,
    so Var(phi) = R @ varcoef @ R.T exactly (no approximation).
    """
    coef = np.asarray(coef, dtype=float)
    R = np.atleast_2d(R)
    r = np.atleast_1d(r)
    return delta_method(coef, varcoef, lambda b: R @ b - r)


def print_results(label, res):
    """Print a formatted table of delta method results."""
    print(f"\n{'─' * 60}")
    print(f"  {label}")
    print(f"{'─' * 60}")
    print(f"  {'':12s} {'Estimate':>12s} {'Std.Err':>12s} "
          f"{'t-stat':>10s} {'p-value':>10s}")
    for i in range(len(res['para'])):
        stars = ('***' if res['p_value'][i] < 0.01 else
                 '**'  if res['p_value'][i] < 0.05 else
                 '*'   if res['p_value'][i] < 0.10 else '')
        print(f"  {'delta'+str(i+1):12s} {res['para'][i]:12.6f} "
              f"{res['se'][i]:12.6f} {res['t_stat'][i]:10.4f} "
              f"{res['p_value'][i]:10.4f} {stars}")
    print()


# ── Simulated data: OLS regression with two regressors ────────────

if __name__ == '__main__':
    rng = np.random.default_rng(42)
    n, k = 500, 2
    X = rng.standard_normal((n, k))
    beta_true = np.array([3.0, 1.0])
    y = X @ beta_true + rng.standard_normal(n)

    # OLS
    coef = np.linalg.solve(X.T @ X, X.T @ y)
    resid = y - X @ coef
    s2 = resid @ resid / (n - k)
    varcoef = s2 * np.linalg.inv(X.T @ X)

    print(f"OLS estimates: beta = [{coef[0]:.4f}, {coef[1]:.4f}]")
    print(f"True values:   beta = [{beta_true[0]:.1f}, {beta_true[1]:.1f}]")

    # ── Case 1: Linear transformation  R @ coef - r ──────────────
    #    Test whether beta1 - beta2 = 0.
    #    Special case: Jacobian is R, so Var(phi) = R @ V @ R' exactly.
    R = np.array([[1, -1]])
    r = np.array([0.0])
    lin = delta_method_linear(coef, varcoef, R, r)
    print_results("Linear: beta1 - beta2", lin)

    # Verify: analytical variance must match
    var_analytical = R @ varcoef @ R.T
    print(f"  Analytical SE:   {np.sqrt(var_analytical[0,0]):.10f}")
    print(f"  Delta method SE: {lin['se'][0]:.10f}  (must match exactly)")

    # ── Case 2: Nonlinear transformation  g(coef) ────────────────
    #    Estimate the ratio beta1 / beta2.
    #    Jacobian G = [1/b2, -b1/b2^2] is evaluated numerically.
    g = lambda b: np.array([b[0] / b[1]])
    nln = delta_method(coef, varcoef, g)
    print_results("Nonlinear: beta1 / beta2", nln)

    # Verify: compare numerical Jacobian against analytical Jacobian
    G_analytical = np.array([1/coef[1], -coef[0]/coef[1]**2])
    se_analytical = np.sqrt(G_analytical @ varcoef @ G_analytical)
    print(f"  Analytical SE:   {se_analytical:.10f}")
    print(f"  Delta method SE: {nln['se'][0]:.10f}  "
          f"(should match to ~10 digits)")
