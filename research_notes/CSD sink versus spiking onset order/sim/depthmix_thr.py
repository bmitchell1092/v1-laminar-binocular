import numpy as np
from math import comb
from scipy.ndimage import gaussian_filter1d
fs=30000; tt=np.arange(-0.1,0.25,1/fs)*1000
def cross(y,frac):
    m=(tt>=0)&(tt<=150); i=np.argmax(np.where(m,y,-np.inf)); thr=frac*y[i]; j=i
    while j>0 and y[j]>thr: j-=1
    return tt[j]+(thr-y[j])/(y[j+1]-y[j])*(tt[j+1]-tt[j]), tt[i]
def resp(rise,shift=0,tau=30):
    t=tt-shift; y=np.clip(t/rise,0,1); r=0.5-0.5*np.cos(np.pi*y)
    return r*np.exp(-np.clip(t-rise,0,None)/tau)
ham=np.array([0.23,0.54,0.23]); avg=np.ones(3)/3
w5=np.convolve(ham,avg); print("effective CSD-contact weights (Hamming 0.23/0.54/0.23 * 3-contact mean):",np.round(w5,3))
fr=[0.1,0.25,0.5,0.75]
print("\nDepth mixing, current only: footprint-averaged current minus centre contact (ms), Gaussian 2 ms on both")
for rise in [6,10]:
  for d in [1,2,4]:
    for a in [1.0,0.5]:
      cen=resp(rise); mix=sum(w5[k+2]*(a**abs(k))*resp(rise,d*abs(k)) for k in range(-2,3))
      cen=gaussian_filter1d(cen,2*fs/1000); mix=gaussian_filter1d(mix,2*fs/1000)
      out=[cross(mix,f)[0]-cross(cen,f)[0] for f in fr]; pk=cross(mix,.5)[1]-cross(cen,.5)[1]
      print("rise %2d d %d ms/100um amp %.1f: "%(rise,d,a)+"  ".join("%d%% %+.2f"%(f*100,o) for f,o in zip(fr,out))+"  peak %+.2f"%pk)
print("\nSign test, two-sided binomial p (n=26):")
for k in [22,17,16,13]:
    p=sum(comb(26,i) for i in range(max(k,26-k),27))/2**26*2
    print(k,"/26 ->",round(min(p,1),4))
