import sys,importlib.util as u,numpy as np
sp=u.spec_from_file_location("h","C:/Users/masse/AppData/Local/Temp/claude/dd1494da25364fcc/09ecfbf3-1f60-4a99-9729-1c85cc9d368d/scratchpad/het_bf.py");h=u.module_from_spec(sp);sp.loader.exec_module(h)
for k in (3,2):
    d=np.loadtxt(f"C:/Users/masse/AppData/Local/Temp/claude/dd1494da25364fcc/09ecfbf3-1f60-4a99-9729-1c85cc9d368d/scratchpad/d{k}.csv",delimiter=",",skiprows=1)
    G=[d[d[:,1]==i+1,0] for i in range(k)]
    n,yb,ss=h.stats(G);g,w=h.g_weights(.5)
    s2=np.ones((1,k));b0=h.logm_given(s2,n,yb,ss,np.array([0.]))[:,0];b1=h.logm_given(s2,n,yb,ss,g)
    print("py",k,"logBF(rho=1)",np.log((np.exp(b1-b0[:,None])*w).sum()))
