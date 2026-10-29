# Derived toy check (not a fit): 10-90% rise time and post-offset half-decay vs contrast,
# for DN (fixed low-pass pool), feedforward shunting (Sit-like), recurrent ORGaNICs (Rawat form),
# and a 2-population SSN. Illustrative parameters only.
import numpy as np
dt = 0.05; T = 400.0; t = np.arange(0, T, dt); TOFF = 200.0
def lp(x, tau):
    y = np.zeros_like(x); a = dt / tau
    for i in range(1, len(x)): y[i] = y[i-1] + a * (x[i-1] - y[i-1])
    return y
stim = ((t >= 0) & (t < TOFF)).astype(float)
def front(tau_in):
    return lp(lp(stim, tau_in), tau_in)

def DN(z, tau2=30.0, n=2.0, sig=0.1):
    P = lp(np.abs(z), tau2); return np.abs(z)**n / (sig**n + P**n)
def SHUNT_FF(z, tau0=40.0, k=None, kz=30.0):
    # tau0 dV/dt = z - g V, g = 1 + kz*z^2  (conductance driven instantaneously by input)
    V = 0.0; out = np.zeros_like(t)
    for i in range(len(t)):
        g = 1 + kz * z[i]**2; V += dt * (z[i] - g * V) / tau0; out[i] = max(V, 0)**2
    return out
def SHUNT_FB(z, tau0=40.0, kz=30.0, tau_g=5.0):
    # conductance driven by output (feedback), low-passed with tau_g
    V = 0.0; G = 0.0; out = np.zeros_like(t)
    for i in range(len(t)):
        g = 1 + G; V += dt * (z[i] - g * V) / tau0; r = max(V, 0)**2
        G += dt * (-G + kz * r * 50.0) / tau_g; out[i] = r
    return out
def ORG(z, ty=2.0, ta=2.0, sig=0.1, b=1.0):
    y = 0.0; a = sig**2; out = np.zeros_like(t)
    for i in range(len(t)):
        sa = np.sqrt(max(a, 0)); dy = (b*z[i] - sa * y) / ty
        da = (-a + sig**2 + (max(y, 0)**2) * a) / ta
        y += dt * dy; a += dt * da; out[i] = max(y, 0)**2
    return out
def SSN(z, tauE=20.0, tauI=10.0, k=0.3, n=2.0, H=20.0):
    WEE, WEI, WIE, WII = 1.25, 0.65, 1.2, 0.5  # illustrative 2-pop values
    rE = rI = 0.0; out = np.zeros_like(t)
    for i in range(len(t)):
        h = H * z[i]
        uE = WEE*rE - WEI*rI + h; uI = WIE*rE - WII*rI + h
        rE += dt * (-rE + k*max(uE, 0)**n) / tauE
        rI += dt * (-rI + k*max(uI, 0)**n) / tauI
        out[i] = rE
    return out

def stats(r):
    on = r[t < TOFF]; tt = t[t < TOFF]; mx = on.max()
    t10 = tt[np.argmax(on >= 0.1*mx)]; t90 = tt[np.argmax(on >= 0.9*mx)]
    r_off = r[np.argmin(np.abs(t - (TOFF - dt)))]
    post = r[t >= TOFF]; tp = t[t >= TOFF]
    # half-decay measured from the time the response starts to fall (peak after offset)
    ipk = np.argmax(post); base = post[ipk]
    ih = ipk + np.argmax(post[ipk:] <= 0.5 * base)
    pk_t = t[np.argmax(on)]; ss = on[tt > 150].mean()
    return t90 - t10, tp[ih] - TOFF, pk_t, mx / max(ss, 1e-12)

for tau_in in ((5.0, 10.0) if __name__ == "__main__" else ()):
    z0 = front(tau_in)
    print(f"--- front end: two low-pass stages, tau={tau_in} ms ---")
    for nm, f in [("DN30", DN), ("SHUNT_FF", SHUNT_FF), ("SHUNT_FB", SHUNT_FB),
                  ("ORG s=.1", lambda z: ORG(z, sig=0.1)), ("ORG s=.3", lambda z: ORG(z, ty=6.0, sig=0.3)),
                  ("SSN", SSN)]:
        row = []
        for c in (0.11, 0.225, 0.45, 0.9):
            rise, hd, pk, ovs = stats(f(c * z0))
            row.append(f"c={c:.3f}: rise={rise:5.1f} halfdec={hd:5.1f} pk/ss={ovs:4.2f}")
        print(f"{nm:9s} | " + " | ".join(row))
