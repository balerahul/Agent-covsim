# HCD model population package

**Generated:** 2025-08-14T07:42:58.548156Z

## Purpose & background

This package contains population-level parameter summaries derived from human-challenge viral-load time-series using an HCD within-host ODE model. It is intended to be shipped to an ABM so that each newly infected agent can be assigned a realistic within-host parameter set. The workflow used here: fit each subject independently (on log-parameter scale), compute population mean and covariance of log-parameters, regularize covariance for safe sampling, and save results. Model fits were performed by least-squares on observed log10 viral load vs time.

## What is inside the JSON package

`mu_log` — dictionary with means of log-parameters (log_beta, log_gamma, log_p).

`sigma_log_marginal` — marginal standard deviations of log-parameters (if >=2 successful fits).

`Sigma_log_matrix` — 3x3 covariance matrix on log-parameter scale (preferred for joint sampling). If `null`, use independent marginals.

`Sigma_log_regularization_jitter` — tiny jitter added to the diagonal to ensure PSD (if any).

`natural_scale_log_normal_moments` — log-normal mean and sd on the natural parameter scale for quick checks.

`model_constants` — U0, I0, V0, c, V_m used during fitting (ABM should use the same or explicitly document differences).

`fitted_subjects` — list of fitted per-subject parameters (log and natural) for diagnostics.


## HCD model (within-host) — equations and meaning

The HCD (simple viral kinetics) ODE used:

```

dU/dt = -beta * U * V

dI/dt = beta * U * V - gamma * I

dV/dt = (p / V_m) * I - c * V

```

- **U**: susceptible/uninfected target cells (cells available to infect)

- **I**: infected cells

- **V**: free virions (viral load)

- **beta**: infection rate constant (per virion per target cell per time)

- **gamma**: infected cell clearance rate (1/time)

- **p**: virion production scale (units matched to V and V_m)

- **c**: virion clearance / apoptosis (given constant)

- **V_m**: scaling constant used in the p*I/V_m term (given)


## How the ABM should use this package (recommended procedure)

1. Load the JSON package.

2. Prefer joint multivariate sampling on the log-scale: sample `log_theta ~ N(mu_log, Sigma_log)` and then set `theta = exp(log_theta)`.

   - This preserves correlations between parameters.

3. If `Sigma_log_matrix` is null or unreliable, sample independent marginals using `mu_log` and `sigma_log_marginal`.

4. For each newly infected agent in the ABM: sample `theta_agent`, store them as that agent's fixed within-host params, and integrate the HCD ODE once for that agent to produce `V_agent(t)` used to compute infectiousness.


## Example ABM sampling code (python)

```python
import numpy as np
# load mu_log_dict, Sigma_log_matrix (or sigma_log_dict)
mu_log_arr = np.array([mu_log['beta'], mu_log['gamma'], mu_log['p']])
Sigma = np.array(Sigma_log_matrix)  # or None
# joint sampling (preferred)
log_theta = np.random.multivariate_normal(mu_log_arr, Sigma)
beta, gamma, p = np.exp(log_theta)
# now integrate HCD with these params for the agent
```

## Notes & caveats

- Parameters are modeled on the **log** scale. The ABM must exponentiate to produce positive parameters.

- If you wish to propagate uncertainty in the population estimates (mu_log, Sigma_log), provide multiple draws of (mu_log, Sigma_log) from the posterior; otherwise simulations treat the package values as fixed population parameters.

- Use the `model_constants` in the JSON for initial conditions and ODE constants to reproduce the fitted curves.


## Contact / provenance

Package created by fit_hcd_package.py. For reproduction, run the script with the same subject CSV files in the data folder. See per-subject fits included in `fitted_subjects` for diagnostics.
