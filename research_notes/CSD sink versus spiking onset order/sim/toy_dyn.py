# Toy comparison of dynamic normalization models for binocular summation (derived, not published).
import numpy as np
dt=0.05; T=300; t=np.arange(0,T,dt)
def lp(x,tau):
    y=np.zeros_like(x); a=dt/tau
    for i in range(1,len(x)): y[i]=y[i-1]+a*(x[i-1]-y[i-1])
    return y
step=(t>=0).astype(float)
h_in=lp(lp(step,5.0),5.0)            # common input filter (2-stage low-pass, tau 5 ms) -> LGN-like drive
sig=0.1
def DN(cL,cR,tau2=30.0,n=2.0):
    L=(cL+cR)*h_in; P=lp(np.abs(L),tau2)
    return np.abs(L)**n/(sig**n+P**n)
def ORG(cL,cR,ty=2.0,ta=2.0):
    z=(cL+cR)*h_in; y=0.0; a=sig**2; out=np.zeros_like(t)
    for i in range(len(t)):
        sa=np.sqrt(max(a,0))
        dy=(-y+z[i]+(1-sa)*y)/ty
        da=(-a+sig**2+(max(y,0)**2)*a)/ta
        y+=dt*dy; a+=dt*da; out[i]=max(y,0)**2
    return out
def SHUNT(cL,cR,tau0=25.0,tau1=5.0):
    z=(cL+cR)*h_in; k=(tau0/tau1-1)  # g/g0 = 1 + k*z^2, so tau=tau1 at z=1
    V=0.0; out=np.zeros_like(t)
    for i in range(len(t)):
        g=1+k*z[i]**2; V+=dt*(z[i]-g*V)/tau0; out[i]=max(V,0)**2
    return out
def RATE(cL,cR,tau=10.0,n=2.0):
    z=(cL+cR)*h_in; tgt=z**n/(sig**n+z**n); return lp(tgt,tau)
def summarize(name,f):
    for c in [0.11,0.45,0.9]:
        b=f(c,c); m=f(c,0.0); bsr=np.where(2*m>1e-9,b/(2*m+1e-12),np.nan)
        # onset = first time monocular response exceeds 10% of its max
        on=t[np.argmax(m>0.1*m.max())]
        ss=np.nanmean(bsr[(t>250)&(t<300)])
        pk=t[np.argmax(m)]
        vals=[bsr[np.argmin(np.abs(t-(on+d)))] for d in (2,5,10,20,35,60)]
        early=vals[0]; half=early-(early-ss)/2
        idx=np.where((t>on)&(bsr<=half))[0]; t50=(t[idx[0]]-on) if len(idx) else np.nan
        print(f"{name:6s} c={c:4.2f} onset={on:5.1f}ms monoPeak@{pk:5.1f} BSR(+2,+5,+10,+20,+35,+60ms)="+",".join(f"{v:4.2f}" for v in vals)+f" BSRss={ss:4.2f} t_half={t50:5.1f}ms")
for nm,f in [("DN30",DN),("ORG",ORG),("SHUNT",SHUNT),("RATE",RATE)]:
    summarize(nm,f)
