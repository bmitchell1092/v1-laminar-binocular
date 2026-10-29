"""Dynamic (delayed) binocular normalization model, after the delayed-normalization model of Zhou, Benson, Kay &
Winawer (2019, PLoS Comput Biol), extended to two eyes.

For dominant-eye contrast c_D and non-dominant-eye contrast c_N switched on at t = 0:
    E(t) = (c_D + a c_N) u(t)          excitatory drive: eye-weighted contrasts summed linearly
    P(t) = (c_D + w c_N) v(t)          normalization pool: its own non-dominant weight w
    R(t) = g E(t)^n / (sigma^n + P(t)^n)
u(t) is the step response of a gamma impulse response (shape K_GAMMA, time constant tau1) delayed by delta;
v(t) is u(t) low-pass filtered by an exponential with time constant tau2 (the normalization delay).
Because the drive is linear in contrast, every condition is a scaled copy of the same u and v, so the model is fast.

Reduced models: static normalization (tau2 -> 0, so v = u) and w = a (the pool weighted like excitation).
Extended model 'tau_c': the pool time constant shortens with pool contrast, tau2 / (1 + beta (c_D + w c_N)).
Derived consequences (not fitted assumptions): the response is linear in the eye-weighted contrast at onset,
before v catches up; the pool reaches sigma sooner at higher contrast; and at high contrast the binocular response
falls below the dominant eye alone only if w > a.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, stats

K_GAMMA = 3.0
PARAMS = ["g", "sigma", "n", "delta", "tau1", "tau2", "a", "w", "beta", "S", "omega", "tau_s"]
BOUNDS = {"g": (0.05, 50.0), "sigma": (0.005, 3.0), "n": (0.5, 6.0), "delta": (5.0, 70.0), "tau1": (0.5, 25.0),
          "tau2": (0.5, 300.0), "a": (0.0, 1.5), "w": (0.0, 3.0),
          "beta": (0.0, 50.0), "S": (0.001, 5.0), "omega": (0.0, 3.0), "tau_s": (0.5, 300.0)}
TWO_STAGE = ["g", "sigma", "n", "delta", "tau1", "tau2", "a", "S", "omega", "tau_s"]


def lowpass(x: np.ndarray, tau: float, dt: float = 1.0) -> np.ndarray:
    """First-order low-pass along the last axis (zero initial state)."""
    alpha = 1 - np.exp(-dt / tau)
    y = np.zeros_like(x, dtype=float)
    for i in range(1, x.shape[-1]):
        y[..., i] = y[..., i - 1] + alpha * (x[..., i] - y[..., i - 1])
    return y


def predict_two_stage(p: dict, t: np.ndarray, cd: np.ndarray, cn: np.ndarray) -> np.ndarray:
    """Dynamic two-stage model (after Meese, Georgeson & Baker 2006, made dynamic).

    Stage 1, interocular gain control per eye, with its own delay tau_s:
        b_D = c_D u / (S + v_s (c_D + omega c_N)),   b_N = a c_N u / (S + v_s (c_N + omega c_D))
    Stage 2, delayed binocular normalization of the sum B = b_D + b_N:
        R = g B^n / (sigma^n + lowpass(B, tau2)^n)
    As S grows, stage 1 becomes linear and the model reduces to single-stage normalization with w = a."""
    u, _ = step_responses(t, p["delta"], p["tau1"], p["tau2"], static=True)
    dt = t[1] - t[0]
    vs = lowpass(u[None, :], p["tau_s"], dt)[0]
    bd = cd[:, None] * u / (p["S"] + vs[None, :] * (cd + p["omega"] * cn)[:, None])
    bn = p["a"] * cn[:, None] * u / (p["S"] + vs[None, :] * (cn + p["omega"] * cd)[:, None])
    b = bd + bn
    pool = lowpass(b, p["tau2"], dt)
    with np.errstate(over="ignore", invalid="ignore"):
        return p["g"] * b ** p["n"] / (p["sigma"] ** p["n"] + pool ** p["n"])


def step_responses(t: np.ndarray, delta: float, tau1: float, tau2: float, static: bool = False):
    """u(t): gamma-CDF step response delayed by delta; v(t): u low-passed by exp(-t / tau2) (v = u if static)."""
    u = stats.gamma.cdf(np.maximum(t - delta, 0.0), K_GAMMA, scale=tau1)
    if static:
        return u, u
    dt = t[1] - t[0]
    k = np.exp(-np.arange(0, 6 * tau2 + dt, dt) / tau2)
    k /= k.sum()
    v = np.convolve(u, k)[: len(t)]
    return u, v


def predict(p: dict, t: np.ndarray, cd: np.ndarray, cn: np.ndarray, static: bool = False) -> np.ndarray:
    """Responses [condition, time] for condition contrasts cd, cn (arrays of equal length).

    If p has "beta" (> 0), the pool time constant shortens with pool contrast, tau2 / (1 + beta (c_D + w c_N)),
    so normalization is faster at higher contrast ('tau_c' model). If p has "S", the two-stage model is used."""
    if "S" in p:
        return predict_two_stage(p, t, cd, cn)
    u, v = step_responses(t, p["delta"], p["tau1"], p["tau2"], static)
    e = (cd + p["a"] * cn)[:, None] * u[None, :]
    pool_c = cd + p["w"] * cn
    if p.get("beta", 0) > 0 and not static:
        dt = t[1] - t[0]
        tau = p["tau2"] / (1 + p["beta"] * pool_c)                      # per condition
        alpha = 1 - np.exp(-dt / tau)
        vv = np.zeros((len(cd), len(t)))
        for i in range(1, len(t)):                                      # first-order low-pass, per condition
            vv[:, i] = vv[:, i - 1] + alpha * (u[i] - vv[:, i - 1])
        q = pool_c[:, None] * vv
    else:
        q = pool_c[:, None] * v[None, :]
    with np.errstate(over="ignore", invalid="ignore"):
        return p["g"] * e ** p["n"] / (p["sigma"] ** p["n"] + q ** p["n"])


@dataclass
class Fit:
    params: dict
    sse: float
    n_points: int
    model: str


def fit(t: np.ndarray, cd: np.ndarray, cn: np.ndarray, y: np.ndarray, model: str = "full",
        starts: int = 6, seed: int = 0) -> Fit:
    """Least-squares fit to responses y [condition, time]. model: 'full', 'static' (no pool delay), 'w=a', or
    'tau_c' (full + a pool time constant that shortens with contrast; the others fix beta = 0), or '2stage'
    (predict_two_stage)."""
    if model == "2stage":
        free = list(TWO_STAGE)
    else:
        free = [k for k in PARAMS if not (model == "static" and k == "tau2") and not (model == "w=a" and k == "w")
                and not (model != "tau_c" and k == "beta") and k not in ("S", "omega", "tau_s")]
    lo = np.array([BOUNDS[k][0] for k in free])
    hi = np.array([BOUNDS[k][1] for k in free])
    ok = np.isfinite(y)
    rng = np.random.default_rng(seed)

    def unpack(x):
        p = dict(zip(free, x))
        if model == "static":
            p["tau2"] = BOUNDS["tau2"][0]
        if model == "w=a":
            p["w"] = p["a"]
        if model != "2stage":
            p.setdefault("beta", 0.0)
        return p

    def resid(x):
        r = predict(unpack(x), t, cd, cn, static=(model == "static")) - y
        return np.where(ok, np.nan_to_num(r, nan=1e3), 0.0).ravel()

    base = {"g": 1.5, "sigma": 0.15, "n": 2.0, "delta": 30.0, "tau1": 4.0, "tau2": 40.0, "a": 0.7, "w": 1.0,
            "beta": 2.0, "S": 0.2, "omega": 1.0, "tau_s": 20.0}
    best = None
    for s in range(starts):
        x0 = np.array([base[k] for k in free], float)
        if s:
            x0 = np.exp(np.log(np.clip(x0, lo + 1e-9, hi)) + rng.normal(0, 0.4, len(x0)))
        x0 = np.clip(x0, lo + 1e-6, hi - 1e-6)
        try:
            res = optimize.least_squares(resid, x0, bounds=(lo, hi), x_scale="jac", max_nfev=3000)
        except ValueError:
            continue
        if best is None or res.cost < best.cost:
            best = res
    return Fit(unpack(best.x), float(2 * best.cost), int(ok.sum()), model)


def sse(p: dict, t, cd, cn, y, model: str) -> float:
    r = predict(p, t, cd, cn, static=(model == "static")) - y
    return float(np.nansum(r ** 2))
