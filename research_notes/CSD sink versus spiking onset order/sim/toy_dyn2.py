# Derived toy sims: (i) does a suppressor in the other eye speed up the suppressed unit? (ii) offset kinetics vs contrast.
import numpy as np
dt=0.05; T=400; t=np.arange(0,T,dt); sig=0.1
def lp(x,tau):
    y=np.zeros_like(x); a=dt/tau
    for i in range(1,len(x)): y[i]=y[i-1]+a*(x[i-1]-y[i-1])
    return y
stim=((t>=0)&(t<200)).astype(float)
h=lp(lp(stim,5.0),5.0)
def DN2(cL,cR,tau2=30.0,n=2.0):     # monocular-L unit, binocular pool (L+R), fixed low-pass pool
    L=cL*h; R=cR*h; P=lp(np.abs(L)+np.abs(R),tau2)
    return np.abs(L)**n/(sig**n+P**n)
def ORG2(cL,cR,ty=2.0,ta=2.0):      # two ORGaNICs units (L,R), shared normalization (W = all ones), Rawat et al. 2024 form
    zL=cL*h; zR=cR*h; y=np.zeros(2); a=np.full(2,sig**2); out=np.zeros_like(t)
    for i in range(len(t)):
        z=np.array([zL[i],zR[i]]); yp=np.maximum(y,0)**2; sa=np.sqrt(np.maximum(a,0))
        dy=(-y+z+(1-sa)*y)/ty; da=(-a+sig**2+np.sum(yp*a))/ta
        y=y+dt*dy; a=a+dt*da; out[i]=max(y[0],0)**2
    return out
def SH2(cL,cR,tau0=25.0,tau1=5.0):  # shunting: conductance driven (instantaneously) by squared input of both eyes
    zL=cL*h; zR=cR*h; k=tau0/tau1-1; V=0.0; out=np.zeros_like(t)
    for i in range(len(t)):
        g=1+k*(zL[i]**2+zR[i]**2); V+=dt*(zL[i]-g*V)/tau0; out[i]=max(V,0)**2
    return out
def stats(r):
    on=t[np.argmax(r>0.1*r[t<200].max())]; pk=t[np.argmax(r*(t<200))]
    r50=t[np.argmax(r>0.5*r[t<200].max())]
    ss=r[(t>150)&(t<200)].mean(); post=r[t>=200]; tp=t[t>=200]
    i37=np.argmax(post<0.37*r[np.argmin(np.abs(t-199.9))]); toff=tp[i37]-200
    return on,r50,pk,ss,toff
for nm,f in [("DN30",DN2),("ORG",ORG2),("SHUNT",SH2)]:
    for cL in [0.11,0.45]:
        for cR in [0.0,0.45,0.9]:
            on,r50,pk,ss,toff=stats(f(cL,cR))
            print(f"{nm:5s} cL={cL:4.2f} cR={cR:4.2f}: t10={on:5.1f} t50={r50:5.1f} tpeak={pk:5.1f} sustained={ss:7.4f} offset_tau~{toff:5.1f}ms")
