import numpy as np
from scipy import signal
from scipy.ndimage import gaussian_filter1d
rng=np.random.default_rng(0)
def halfmax(t,y,frac=0.5,t0=0,t1=150):
    m=(t>=t0)&(t<=t1); i=np.argmax(y*m); thr=frac*y[i]
    j=i
    while j>0 and y[j]>thr: j-=1
    return t[j]+(thr-y[j])/(y[j+1]-y[j])*(t[j+1]-t[j]), t[i]
# ---- 1. analog causal 4th-order Butterworth 500 Hz (approximated digitally at 30 kHz) and 0.3 Hz 1st-order HP
fs=30000; tt=np.arange(-0.1,0.3,1/fs)*1000
b,a=signal.butter(4,2*np.pi*500,analog=True); w,gd=None,None
sys=signal.lti(b,a); w=np.array([1e-3]); 
# DC group delay of analog Butterworth
den=np.poly1d(a); print("4th-order 500 Hz Butterworth DC group delay (ms): %.3f"%(1000*a[-2]/a[-1]))
sos_lp=signal.butter(4,500,fs=fs,output='sos')
def ramp(rise):  # linear-ish smooth rise: raised cosine from 0 to 1 over 'rise' ms (10-90% ~0.6*rise), starting at 0
    y=np.clip(tt/rise,0,1); return 0.5-0.5*np.cos(np.pi*y)
for rise in [4,8,16]:
    x=ramp(rise)*np.exp(-np.clip(tt-rise,0,None)/40)
    yc=signal.sosfilt(sos_lp,x); yz=signal.sosfiltfilt(sos_lp,x)
    print("rise %2d ms: causal500 half-max shift %+.2f ms, peak shift %+.2f; zero-phase500 half-max %+.2f"%(rise,
        halfmax(tt,yc)[0]-halfmax(tt,x)[0], halfmax(tt,yc)[1]-halfmax(tt,x)[1], halfmax(tt,yz)[0]-halfmax(tt,x)[0]))
# ---- 2. zero-phase 200 Hz 4th-order envelope LP (sosfiltfilt) and Gaussian sigma 2 ms on half-max/10% onset for different rises
sos200=signal.butter(4,200,fs=fs,output='sos')
for rise in [2,4,8]:
    x=ramp(rise)*np.exp(-np.clip(tt-rise,0,None)/40)
    y=signal.sosfiltfilt(sos200,x)
    g=gaussian_filter1d(x,2*fs/1000)
    on=lambda yy,f: halfmax(tt,yy,f)[0]-halfmax(tt,x,f)[0]
    print("rise %d: LP200 zp shift 10%% %+.2f 50%% %+.2f | Gauss2ms shift 10%% %+.2f 50%% %+.2f peak %+.2f"%(rise,on(y,.1),on(y,.5),on(g,.1),on(g,.5),halfmax(tt,g)[1]-halfmax(tt,x)[1]))
# step: how early does a zero-phase 200 Hz envelope / sigma-2ms Gaussian first exceed 1%,5% of a step at t=0
x=(tt>=0).astype(float); y=signal.sosfiltfilt(sos200,x); g=gaussian_filter1d(x,2*fs/1000)
for f in [0.01,0.05,0.1]:
    print("step at 0: LP200zp crosses %d%% at %+.2f ms; Gauss2ms at %+.2f ms"%(100*f,tt[np.argmax(y>f)],tt[np.argmax(g>f)]))
# ---- 3. rectified-envelope nonlinearity: aMUA vs underlying rate, half-max of each (no filter differences)
sos_bp=signal.butter(2,[500,5000],btype='band',fs=fs,output='sos')  # 2nd-order design -> 4th-order bandpass
n=int(0.004*fs); sp=np.zeros(n); k=np.arange(n)/fs*1000
sp=-np.exp(-((k-0.5)/0.15)**2)+0.4*np.exp(-((k-0.9)/0.3)**2)  # biphasic spike ~1 ms
T=len(tt); ntr=150
def rate(rise, base, peak):
    r=ramp(rise)*np.exp(-np.clip(tt-rise,0,None)/25); r=r/r.max()
    return base+(peak-base)*np.clip(r,0,None)
def amua(rt_hz, noise_sd, amp_sd=0.3):
    acc=np.zeros(T)
    for _ in range(ntr):
        cnt=rng.poisson(rt_hz/fs)
        raw=np.zeros(T)
        idx=np.repeat(np.arange(T),cnt)
        np.add.at(raw,idx,1+amp_sd*rng.standard_normal(len(idx)))
        raw=np.convolve(raw,sp,mode='same')+noise_sd*rng.standard_normal(T)
        e=np.abs(signal.sosfiltfilt(sos_bp,raw))
        acc+=signal.sosfiltfilt(sos200,e)
    return acc/ntr
def bsub(t,y): return y-y[(t<-10)&(t>-90)].mean()
for label,base,peak,noise in [("sparse,low noise (pop 2k Hz base, 20k peak)",2000,20000,0.05),
                              ("dense,low noise (20k base,200k peak)",20000,200000,0.05),
                              ("noise-dominated (2k,20k, noise sd 1.0)",2000,20000,1.0),
                              ("dense, zero base (200 Hz base, 200k peak)",200,200000,0.05)]:
    for rise in [4,8]:
        r=rate(rise,base,peak); m=amua(r,noise)
        mg=gaussian_filter1d(bsub(tt,m),2*fs/1000); rg=gaussian_filter1d(r-base,2*fs/1000)
        hm=halfmax(tt,mg)[0]-halfmax(tt,rg)[0]; q=halfmax(tt,mg,0.1)[0]-halfmax(tt,rg,0.1)[0]
        print("%-46s rise %d: aMUA-rate half-max %+.2f ms, 10%% %+.2f ms, peak %+.2f ms"%(label,rise,hm,q,halfmax(tt,mg)[1]-halfmax(tt,rg)[1]))
# ---- 4. same onset, different rise times: half-max difference
for r1,r2 in [(4,8),(4,12),(6,12)]:
    a1=gaussian_filter1d(ramp(r1)*np.exp(-np.clip(tt-r1,0,None)/40),2*fs/1000)
    a2=gaussian_filter1d(ramp(r2)*np.exp(-np.clip(tt-r2,0,None)/40),2*fs/1000)
    print("same true onset, rise %d vs %d ms: half-max diff %+.2f ms, 10%% diff %+.2f ms, peak diff %+.2f"%(r1,r2,halfmax(tt,a1)[0]-halfmax(tt,a2)[0],halfmax(tt,a1,.1)[0]-halfmax(tt,a2,.1)[0],halfmax(tt,a1)[1]-halfmax(tt,a2)[1]))
# ---- 5. 3-contact average: centre current onset 0, neighbours onset +d ms, same shape
for d in [2,4,6]:
    c=ramp(8)*np.exp(-np.clip(tt-8,0,None)/40); nb=np.interp(tt-d,tt,c)
    for wn in [0.5,1.0]:
        avg=(c+wn*2*nb)/(1+2*wn)  # wn = neighbour amplitude relative to centre
        print("neighbours +%d ms, rel amp %.1f: 3-ch average half-max shift %+.2f ms"%(d,wn,halfmax(tt,avg)[0]-halfmax(tt,c)[0]))
