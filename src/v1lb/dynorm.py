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
from scipy import optimize, signal, stats

K_GAMMA = 3.0
PARAMS = ["g", "sigma", "n", "delta", "tau1", "tau2", "a", "w", "beta", "S", "omega", "tau_s"]
BOUNDS = {"g": (0.05, 50.0), "sigma": (0.005, 3.0), "n": (0.5, 6.0), "delta": (5.0, 70.0), "tau1": (0.5, 25.0),
          "tau2": (0.5, 300.0), "a": (0.0, 1.5), "w": (0.0, 3.0),
          "beta": (0.0, 50.0), "S": (0.001, 5.0), "omega": (0.0, 3.0), "tau_s": (0.5, 300.0)}
TWO_STAGE = ["g", "sigma", "n", "delta", "tau1", "tau2", "a", "S", "omega", "tau_s"]
TWO_STAGE_FULL = TWO_STAGE + ["w"]               # 'full2': stage-2 pool with its own non-dominant weight
FULL2_BOUNDS = {"sigma": (1e-5, 3.0), "S": (0.001, 500.0)}   # wide enough to reach the single-stage limit
DS = ["g", "sigma", "n", "delta", "tau1", "tau2", "a", "rho", "gamma", "tau_g"]   # dynamic Ding-Sperling
DS_BOUNDS = {"rho": (0.01, 200.0), "gamma": (0.5, 4.0), "tau_g": (1.0, 400.0)}
ORG = ["g", "b", "sigma", "n", "delta", "tau1", "tau_y", "tau_a", "a", "w"]        # binocular ORGaNICs
ORG_BOUNDS = {"b": (0.01, 20.0), "sigma": (0.005, 1.0), "n": (0.5, 4.0), "tau_y": (0.3, 30.0), "tau_a": (0.5, 200.0)}


def lowpass(x: np.ndarray, tau: float, dt: float = 1.0) -> np.ndarray:
    """First-order low-pass along the last axis (zero initial state)."""
    alpha = 1 - np.exp(-dt / tau)
    return signal.lfilter([alpha], [1.0, -(1.0 - alpha)], np.asarray(x, float), axis=-1)


def predict_two_stage(p: dict, t: np.ndarray, cd: np.ndarray, cn: np.ndarray) -> np.ndarray:
    """Dynamic two-stage model (after Meese, Georgeson & Baker 2006, made dynamic).

    Stage 1, interocular gain control per eye, with its own delay tau_s:
        b_D = c_D u / (S + v_s (c_D + omega c_N)),   b_N = a c_N u / (S + v_s (c_N + omega c_D))
    Stage 2, delayed binocular normalization of the sum B = b_D + b_N:
        R = g B^n / (sigma^n + lowpass(B, tau2)^n)
    As S grows, stage 1 becomes linear and the model reduces to single-stage normalization with w = a; with a
    pool weight "w" in p ('full2'), the pool is lowpass(b_D + w b'_N) and the limit is single-stage 'full'."""
    u, _ = step_responses(t, p["delta"], p["tau1"], p["tau2"], static=True)
    dt = t[1] - t[0]
    vs = lowpass(u[None, :], p["tau_s"], dt)[0]
    bd = cd[:, None] * u / (p["S"] + vs[None, :] * (cd + p["omega"] * cn)[:, None])
    bn1 = cn[:, None] * u / (p["S"] + vs[None, :] * (cn + p["omega"] * cd)[:, None])   # unweighted non-dominant
    b = bd + p["a"] * bn1
    # 'full2' (p has "w"): the pool weights the non-dominant stage-1 output by w; otherwise the pool is b itself
    pool = lowpass(bd + p["w"] * bn1 if "w" in p else b, p["tau2"], dt)
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


def predict_ds(p: dict, t: np.ndarray, cd: np.ndarray, cn: np.ndarray) -> np.ndarray:
    """Dynamic Ding-Sperling (2006) gain-control combination followed by delayed normalization.

    x_D = c_D u, x_N = a c_N u; gain-control energy eps_e = (rho lowpass(x_e, tau_g))^gamma (Ding & Sperling filter
    the energy temporally); Gamma = sum_e x_e (1 + eps_e) / (1 + eps_D + eps_N), which gives Gamma(I, 0) = I and
    Gamma(I, I) = I (2 + 2 eps) / (1 + 2 eps) as in the paper; R = g Gamma^n / (sigma^n + lowpass(Gamma, tau2)^n).
    Linear summation while eps is small (onset), contrast-weighted averaging once eps >> 1 (late), so with a < 1
    the binocular response can fall below the dominant eye alone."""
    u, _ = step_responses(t, p["delta"], p["tau1"], p["tau2"], static=True)
    dt = t[1] - t[0]
    xd = cd[:, None] * u[None, :]
    xn = p["a"] * cn[:, None] * u[None, :]
    ed = (p["rho"] * lowpass(xd, p["tau_g"], dt)) ** p["gamma"]
    en = (p["rho"] * lowpass(xn, p["tau_g"], dt)) ** p["gamma"]
    gam = (xd * (1 + ed) + xn * (1 + en)) / (1 + ed + en)
    pool = lowpass(gam, p["tau2"], dt)
    with np.errstate(over="ignore", invalid="ignore"):
        return p["g"] * gam ** p["n"] / (p["sigma"] ** p["n"] + pool ** p["n"])


def predict_org(p: dict, t: np.ndarray, cd: np.ndarray, cn: np.ndarray, dur: float | None = None) -> np.ndarray:
    """Binocular ORGaNICs (Heeger & Zemlianova 2020, in the form of Rawat, Heeger & Martiniani 2024).

    Two populations, D (recorded) driven by z_D = (c_D + a c_N) u(t) and N driven by z_N = (c_N + a c_D) u(t):
        tau_y dy_i/dt = -sqrt(a_i) y_i + b z_i
        tau_a da_i/dt = -a_i + sigma^2 + y_i^2 a_i + w y_j^2 a_j
    R = g max(y_D, 0)^n. The effective time constant tau_y / sqrt(a_i) is long at rest and shortens as the
    modulator grows with the pooled (both-eye) response, so weak inputs rise slowly and late, and normalization
    builds up after onset. At steady state with w = 0, y^2 = (b z)^2 / (sigma^2 + (b z)^2). u(t) is the gamma step
    of 'full' (delta, tau1), switched off after `dur` ms if given. Integrated exactly for coefficients held
    constant over each time step."""
    dt = t[1] - t[0]
    u = stats.gamma.cdf(np.maximum(t - p["delta"], 0.0), K_GAMMA, scale=p["tau1"])
    if dur is not None:
        u = u - stats.gamma.cdf(np.maximum(t - p["delta"] - dur, 0.0), K_GAMMA, scale=p["tau1"])
    zd, zn = (cd + p["a"] * cn)[:, None] * u, (cn + p["a"] * cd)[:, None] * u
    s2, b, ty, ta, w = p["sigma"] ** 2, p["b"], p["tau_y"], p["tau_a"], p["w"]
    yd, yn = np.zeros(len(cd)), np.zeros(len(cd))
    ad, an = np.full(len(cd), s2), np.full(len(cd), s2)
    out = np.empty((len(cd), len(t)))

    def relax(x, target, rate):                      # exact step of dx/dt = rate (target - x), rate > 0
        return target + (x - target) * np.exp(-rate * dt)

    def modulator(a, y, src):                        # tau_a da/dt = -a (1 - y^2) + src
        k = (1 - y ** 2) / ta
        grow = k <= 1e-9
        safe_k = np.where(grow, 1.0, k)
        nxt = np.where(grow, a + dt * (src / ta - k * a), relax(a, src / (safe_k * ta), safe_k))
        return np.clip(nxt, s2 * 1e-3, 1e8)

    for i in range(len(t)):
        sd, sn = np.sqrt(ad), np.sqrt(an)
        yd_new = relax(yd, b * zd[:, i] / sd, sd / ty)
        yn_new = relax(yn, b * zn[:, i] / sn, sn / ty)
        ad, an = modulator(ad, yd_new, s2 + w * yn_new ** 2 * an), modulator(an, yn_new, s2 + w * yd_new ** 2 * ad)
        yd, yn = yd_new, yn_new
        out[:, i] = yd
    with np.errstate(over="ignore", invalid="ignore"):
        return p["g"] * np.maximum(out, 0.0) ** p["n"]


def interaction_residuals(pred: np.ndarray, y: np.ndarray, cd: np.ndarray, cn: np.ndarray) -> np.ndarray:
    """Residuals for the interaction-targeted loss: monocular time courses (one eye at 0) and the binocular
    interaction I(c_D, c_N, t) = R(c_D, c_N, t) - R(c_D, 0, t), each part divided by sqrt(data variance x number
    of points), so both parts weigh equally. Sum of squares = the loss."""
    mono = (cd == 0) | (cn == 0)
    bino = np.flatnonzero((cd > 0) & (cn > 0))
    ref = np.array([np.flatnonzero((cd == cd[i]) & (cn == 0))[0] for i in bino])
    ym, pm = y[mono], pred[mono]
    yi, pi = y[bino] - y[ref], pred[bino] - pred[ref]
    out = []
    for yy, pp in ((ym, pm), (yi, pi)):
        ok = np.isfinite(yy)
        scale = np.sqrt(np.nanvar(yy) * ok.sum()) + 1e-12
        out.append(np.where(ok, np.nan_to_num(pp - yy, nan=1e3), 0.0).ravel() / scale)
    return np.concatenate(out)


def loss_value(p: dict, t, cd, cn, y, loss: str = "interaction") -> float:
    pred = predict(p, t, cd, cn)
    if loss == "interaction":
        return float(np.sum(interaction_residuals(pred, y, cd, cn) ** 2))
    return float(np.nansum((pred - y) ** 2))


def predict(p: dict, t: np.ndarray, cd: np.ndarray, cn: np.ndarray, static: bool = False) -> np.ndarray:
    """Responses [condition, time] for condition contrasts cd, cn (arrays of equal length).

    If p has "beta" (> 0), the pool time constant shortens with pool contrast, tau2 / (1 + beta (c_D + w c_N)),
    so normalization is faster at higher contrast ('tau_c' model). If p has "S", the two-stage model is used; if it
    has "tau_a", binocular ORGaNICs (predict_org)."""
    if "tau_a" in p:
        return predict_org(p, t, cd, cn)
    if "S" in p:
        return predict_two_stage(p, t, cd, cn)
    if "rho" in p:
        return predict_ds(p, t, cd, cn)
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
        starts: int = 6, seed: int = 0, init: dict | None = None, loss: str = "sse") -> Fit:
    """Least-squares fit to responses y [condition, time]. model: 'full', 'static' (no pool delay), 'w=a', or
    'tau_c' (full + a pool time constant that shortens with contrast; the others fix beta = 0), or '2stage'
    (predict_two_stage), or 'full2' (two-stage with a stage-2 pool weight w; nests 'full'). `init` adds a first
    start at the given parameters (e.g. a single-stage solution embedded at the S bound)."""
    if model == "2stage":
        free = list(TWO_STAGE)
    elif model == "full2":
        free = list(TWO_STAGE_FULL)
    elif model == "ds":
        free = list(DS)
    elif model == "org":
        free = list(ORG)
    else:
        free = [k for k in PARAMS if not (model == "static" and k == "tau2") and not (model == "w=a" and k == "w")
                and not (model != "tau_c" and k == "beta") and k not in ("S", "omega", "tau_s")]
    bounds = (BOUNDS | (FULL2_BOUNDS if model == "full2" else {}) | (DS_BOUNDS if model == "ds" else {})
              | (ORG_BOUNDS if model == "org" else {}))
    lo = np.array([bounds[k][0] for k in free])
    hi = np.array([bounds[k][1] for k in free])
    ok = np.isfinite(y)
    rng = np.random.default_rng(seed)

    def unpack(x):
        p = dict(zip(free, x))
        if model == "static":
            p["tau2"] = BOUNDS["tau2"][0]
        if model == "w=a":
            p["w"] = p["a"]
        if model not in ("2stage", "full2", "ds", "org"):
            p.setdefault("beta", 0.0)
        return p

    def resid(x):
        pred = predict(unpack(x), t, cd, cn, static=(model == "static"))
        if loss == "interaction":
            return interaction_residuals(pred, y, cd, cn)
        r = pred - y
        return np.where(ok, np.nan_to_num(r, nan=1e3), 0.0).ravel()

    base = {"g": 1.5, "sigma": 0.15, "n": 2.0, "delta": 30.0, "tau1": 4.0, "tau2": 40.0, "a": 0.7, "w": 1.0,
            "beta": 2.0, "S": 0.2, "omega": 1.0, "tau_s": 20.0, "rho": 5.0, "gamma": 1.2, "tau_g": 50.0,
            "b": 0.5, "tau_y": 2.0, "tau_a": 10.0}
    if model == "org":
        base = base | {"sigma": 0.1, "delta": 25.0}
    best = None
    for s in range(starts + (1 if init else 0)):
        first_base = 1 if init else 0                     # with init: start 0 = init, start 1 = base, then random
        if init and s == 0:
            x0 = np.array([init.get(k, base[k]) for k in free], float)
        else:
            x0 = np.array([base[k] for k in free], float)
            if s > first_base:
                x0 = np.exp(np.log(np.clip(x0, lo + 1e-9, hi)) + rng.normal(0, 0.4, len(x0)))
        x0 = np.clip(x0, lo + 1e-6, hi - 1e-6)
        try:
            res = optimize.least_squares(resid, x0, bounds=(lo, hi), x_scale="jac", max_nfev=3000)
        except ValueError:
            continue
        if best is None or res.cost < best.cost:
            best = res
    return Fit(unpack(best.x), float(2 * best.cost), int(ok.sum()), model)   # sse = loss value for loss='interaction'


def sse(p: dict, t, cd, cn, y, model: str) -> float:
    r = predict(p, t, cd, cn, static=(model == "static")) - y
    return float(np.nansum(r ** 2))
