import numpy as np
from scipy import signal
from scipy.ndimage import gaussian_filter1d
rng=np.random.default_rng(1)
fs=30000; tt=np.arange(-0.1,0.25,1/fs)*1000; T=len(tt)
def cross(t,y,frac=0.5,t1=150):
    m=(t>=0)&(t<=t1); i=np.argmax(np.where(m,y,-np.inf)); thr=frac*y[i]; j=i
    while j>0 and y[j]>thr: j-=1
    return t[j]+(thr-y[j])/(y[j+1]-y[j])*(t[j+1]-t[j]), t[i]
def ramp(rise): y=np.clip(tt/rise,0,1); return 0.5-0.5*np.cos(np.pi*y)
sos200=signal.butter(4,200,fs=fs,output='sos'); sos_bp=signal.butter(2,[500,5000],btype='band',fs=fs,output='sos')
x=(tt>=0).astype(float); y=signal.sosfiltfilt(sos200,x); pre=y[tt<0]
print("zp 200Hz LP step: pre-onset max %+.4f at %.2f ms, min %+.4f at %.2f ms; 5%% crossing %.2f ms (last upward crossing before 0.5)"%(pre.max(),tt[tt<0][pre.argmax()],pre.min(),tt[tt<0][pre.argmin()],cross(tt,y,0.05)[0]))
yc=signal.sosfilt(sos200,x); print("causal 200Hz LP step: 50%% at %+.2f ms, 10%% at %+.2f ms"%(tt[np.argmax(yc>0.5)],tt[np.argmax(yc>0.1)]))
k=np.arange(int(0.004*fs))/fs*1000; sp=-np.exp(-((k-0.5)/0.15)**2)+0.4*np.exp(-((k-0.9)/0.3)**2); k0=int(0.5e-3*fs)
def rate(rise,base,peak):
    r=ramp(rise)*np.exp(-np.clip(tt-rise,0,None)/25); return base+(peak-base)*r/r.max()
def amua(rt,noise,ntr=400):
    acc=np.zeros(T)
    for _ in range(ntr):
        cnt=rng.poisson(rt/fs); idx=np.repeat(np.arange(T),cnt); raw=np.zeros(T)
        np.add.at(raw,idx,1+0.3*rng.standard_normal(len(idx)))
        raw=np.convolve(raw,sp,'full')[k0:k0+T]+noise*rng.standard_normal(T)
        acc+=signal.sosfiltfilt(sos200,np.abs(signal.sosfiltfilt(sos_bp,raw)))
    return acc/ntr
for label,base,peak,noise in [("sparse, low noise",500,5000,0.02),("sparse, noise-dominated",500,5000,1.0),("dense, low noise",20000,200000,0.02),("dense, near-zero base",200,200000,0.02)]:
    for rise in [4,8]:
        r=rate(rise,base,peak); m=amua(r,noise); m=m-m[(tt<-10)&(tt>-90)].mean()
        mg=gaussian_filter1d(m,2*fs/1000); rg=gaussian_filter1d(r-base,2*fs/1000)
        print("%-24s rise %d: aMUA-rate 50%% %+.2f, 10%% %+.2f, peak %+.2f ms; aMUA peak/base ratio %.2f"%(label,rise,cross(tt,mg)[0]-cross(tt,rg)[0],cross(tt,mg,.1)[0]-cross(tt,rg,.1)[0],cross(tt,mg)[1]-cross(tt,rg)[1], (m.max()+0)/1))
# analytic compressive check: envelope = sqrt(base + r) - sqrt(base)
for b in [0.0,0.1,1.0]:
    for rise in [4,8]:
        r=ramp(rise)*np.exp(-np.clip(tt-rise,0,None)/25); r/=r.max(); e=np.sqrt(b+r)-np.sqrt(b)
        print("sqrt model base/peak %.1f rise %d: half-max shift %+.2f ms"%(b,rise,cross(tt,e)[0]-cross(tt,r)[0]))
