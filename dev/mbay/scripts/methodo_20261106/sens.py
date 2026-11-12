import importlib.util as u,numpy as np
sp=u.spec_from_file_location("h","C:/Users/masse/AppData/Local/Temp/claude/dd1494da25364fcc/09ecfbf3-1f60-4a99-9729-1c85cc9d368d/scratchpad/het_bf.py");h=u.module_from_spec(sp);sp.loader.exec_module(h)
rng=np.random.default_rng(7)
def bf(G,a,r=.5,npts=801):
    n,yb,ss=h.stats(G); U,dV=h.grid_u(n,ss,npts,width=14.)
    s2=np.exp(2*np.c_[U,np.zeros(len(U))]); s=(s2.mean(1))
    b=h.logm_given(s2,n,yb,ss,np.array([0.]))[:,0]; lw=-.5*((n-1)*np.log(s2)).sum(1)+b
    rho=1/(1+np.exp(2*U[:,0])); lw=lw+a*(np.log(rho)+np.log1p(-rho))
    w=np.exp(lw-lw.max()); g,wg=h.g_weights(r); s2r=s2/s[:,None]
    b0=h.logm_given(s2r,n,yb,ss,np.array([0.]))[:,0]; b1=h.logm_given(s2r,n,yb,ss,g)
    return (w*(np.exp(b1-b0[:,None])*wg).sum(1)).sum()/w.sum()
for lab,ns,m,sd in [("n=5,5 sd=1,1",(5,5),(0,1),(1,1)),("n=5,5 sd=1,3",(5,5),(0,1),(1,3)),("n=15,25 sd=1,2",(15,25),(0,.8),(1,2)),("n=40,40 sd=1,3",(40,40),(0,.8),(1,3)),("n=5,40 sd=4,1",(5,40),(0,1),(4,1))]:
    G=[rng.normal(a,b,c) for c,a,b in zip(ns,m,sd)]
    print(lab," ".join(f"Beta({a},{a})={bf(G,a):.3f}" for a in (0,.5,1.5,5)),flush=True)
