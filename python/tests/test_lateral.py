import numpy as np
import pytest

from helpers import NG_WEDGE
from pydelivery import Inversion
from pydelivery.lateral import LateralCoupling, correlated_chain, invert_lateral


def _dense(lik_J, lik_h, Sigma, rho):
    n, k = lik_h.shape
    C0 = np.zeros((n * k, n * k))
    for i in range(n):
        for j in range(n):
            C0[i * k:(i + 1) * k, j * k:(j + 1) * k] = rho ** abs(i - j) * Sigma
    J = np.linalg.inv(C0)
    h = np.zeros(n * k)
    for i in range(n):
        J[i * k:(i + 1) * k, i * k:(i + 1) * k] += lik_J[i]
        h[i * k:(i + 1) * k] += lik_h[i]
    C = np.linalg.inv(J)
    return C @ h, C


def test_chain_smoother_is_exact():
    rng = np.random.default_rng(1)
    n, k, rho = 9, 3, 0.7
    A = rng.normal(size=(k, k))
    Sigma = A @ A.T + k * np.eye(k)
    lik_J = np.array([np.diag(rng.uniform(0, 2, k)) for _ in range(n)])
    lik_h = rng.normal(size=(n, k))
    mu, C = _dense(lik_J, lik_h, Sigma, rho)
    oJ, oh, pm, pc = correlated_chain(lik_J, lik_h, Sigma, rho)
    np.testing.assert_allclose(pm.ravel(), mu, atol=1e-10)
    np.testing.assert_allclose(pc, [C[i * k:(i + 1) * k, i * k:(i + 1) * k] for i in range(n)], atol=1e-10)
    # own prior + own data + information from the others = smoothed posterior
    rec = np.linalg.inv(np.linalg.inv(Sigma) + lik_J + oJ)
    np.testing.assert_allclose(rec, pc, atol=1e-10)


def test_no_data_gives_no_lateral_information():
    Sigma = np.diag([4.0, 1.0])
    oJ, oh, _, pc = correlated_chain(np.zeros((6, 2, 2)), np.zeros((6, 2)), Sigma, 0.9)
    np.testing.assert_allclose(oJ, 0, atol=1e-9)
    np.testing.assert_allclose(oh, 0, atol=1e-12)
    np.testing.assert_allclose(pc, np.broadcast_to(Sigma, (6, 2, 2)))  # prior unchanged


def test_invert_lateral_runs_and_reduces_jitter():
    inv = Inversion.from_files(NG_WEDGE / "ModelDescription.xml", NG_WEDGE / "prior_traces.su")
    kw = dict(burn_in=200, thin=5, n_chains=16)
    base = inv.run(80, seed=0, **kw)
    run = invert_lateral(inv, LateralCoupling(("vp_m", "rho_m"), correlation_length=30),
                         n_samples=80, seed=1, start=base, **kw)
    assert len(run.results) == inv.n_traces and ("vp_m", 0) in run.keys
    k = run.keys.index(("vp_m", 0))
    assert np.std(np.diff(run.smoothed_mean[:, k])) < np.std(np.diff(run.local_mean[:, k]))
    assert np.all(run.smoothed_sd[:, k] < run.local_sd[:, k] + 1e-9)
    assert run.results[0].modes[0].model.lateral_keys  # the extra prior term is present
