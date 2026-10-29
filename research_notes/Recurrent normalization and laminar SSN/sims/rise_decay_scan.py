# Derived toy scan (not a fit to data): can each model class reproduce the team's onset/offset
# diagnostics (10-90% rise 36 ms at c=.11 -> 18 ms at c=.9; half-decay ~27-32 ms, ~contrast-independent)
# when given a saturating front end f(c) = c/(c + c50in) and a two-stage low-pass input filter?
import numpy as np, itertools, sys
sys.path.insert(0, '.')
from rise_decay import t, dt, TOFF, lp, stim, DN, ORG, SHUNT_FB, SSN, stats

def front(tau_in): return lp(lp(stim, tau_in), tau_in)
CS = (0.11, 0.45, 0.9)
target = dict(r11=36.0, r90=18.0, hd=29.5)
def score(f, z0, c50):
    out = []
    for c in CS:
        drive = c / (c + c50) if c50 > 0 else c
        rise, hd, pk, ovs = stats(f(drive * z0))
        out.append((rise, hd))
    r11, r90 = out[0][0], out[-1][0]; hds = [o[1] for o in out]
    err = ((r11-36)/36)**2 + ((r90-18)/18)**2 + sum(((h-29.5)/29.5)**2 for h in hds)/3
    return err, out

fronts = {tau: front(tau) for tau in (5.0, 10.0, 15.0)}
best = {}
models = {
  'DN':    (lambda p: (lambda z: DN(z, tau2=p[0], sig=p[1])), list(itertools.product((15., 30., 60.), (0.05, 0.1, 0.3)))),
  'ORG':   (lambda p: (lambda z: ORG(z, ty=p[0], ta=2.0, sig=p[1])), list(itertools.product((2., 4., 8., 12.), (0.1, 0.2, 0.3, 0.5)))),
  'SHFB':  (lambda p: (lambda z: SHUNT_FB(z, tau0=p[0], kz=p[1])), list(itertools.product((20., 40., 60.), (3., 10., 30., 100.)))),
  'SSN':   (lambda p: (lambda z: SSN(z, tauE=p[0], tauI=p[0]/2, H=p[1])), list(itertools.product((10., 20., 40.), (10., 20., 40.)))),
}
for nm, (mk, grid) in models.items():
    res = []
    for p in grid:
        for tau_in, z0 in fronts.items():
            for c50 in (0.0, 0.1, 0.3):
                err, out = score(mk(p), z0, c50)
                res.append((err, p, tau_in, c50, out))
    res.sort(key=lambda r: r[0])
    print(f"== {nm}: best 3 (err, params, tau_in, c50_in, [(rise,halfdecay) at c=.11,.45,.9]) ==")
    for r in res[:3]:
        print(f"  err={r[0]:.3f} p={r[1]} tau_in={r[2]} c50in={r[3]} " + ", ".join(f"({a:.1f},{b:.1f})" for a, b in r[4]))
