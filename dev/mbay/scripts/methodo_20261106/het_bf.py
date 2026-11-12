"""BF heteroscedastique : y_ij ~ N(mu+a_i, sigma_i^2), a|sigma,g ~ N(0,g tau^2 P),
g~IG(1/2,r^2/2), prior plat sur (mu, log sigma_i). Conditionnellement a (rapports rho, g)
tout est conjugue ; on integre u=log(sigma_i/sigma_k) (k-1 dim) et log g par quadrature.
Controle : estimateur independant E_{sigma|y,H0}[BF(sigma)] par tirages exacts."""
import numpy as np, sys, time

def taus2(s2, n, kind):
    # s2 : (...,k) variances relatives ; renvoie tau^2 homogene de degre 1 en s2
    if kind == "arith": return s2.mean(-1)
    if kind == "geom": return np.exp(np.log(s2).mean(-1))
    if kind == "harm": return 1.0 / (1.0 / s2).mean(-1)
    if kind == "nweight": return (s2 * n).sum(-1) / n.sum()
    if kind == "n1weight": return (s2 * (n - 1)).sum(-1) / (n - 1).sum()
    if kind == "welch":  # echelle unite d'information : nbar * mean(s2/n)
        return n.mean() * (s2 / n).mean(-1)
    raise ValueError(kind)

_gl_t = np.linspace(-14, 26, 121)  # log g
_gl_dt = _gl_t[1] - _gl_t[0]
def g_weights(r):
    g = np.exp(_gl_t)
    w = r / np.sqrt(2 * np.pi) * g ** (-1.5) * np.exp(-r * r / (2 * g)) * g * _gl_dt
    return g, w

def stats(groups):
    n = np.array([len(x) for x in groups], float)
    yb = np.array([np.mean(x) for x in groups])
    ss = np.array([np.sum((x - np.mean(x)) ** 2) for x in groups])
    return n, yb, ss

def logm_given(s2, n, yb, ss, g):
    """s2 (m,k) variances relatives, g (q,) ; log[ |V|^-1/2 (1'V^-1 1)^-1/2 Q^-(N-1)/2 ]
    renvoie (m,q). g=0 donne H0."""
    N = n.sum(); k = len(n)
    D = s2 / n  # (m,k)
    ssw = (ss / s2).sum(-1)  # (m,)
    out = np.empty((s2.shape[0], len(g)))
    one = np.ones(k)
    P = np.eye(k) - np.ones((k, k)) / k
    for j, gj in enumerate(g):
        V = np.einsum('mi,ij->mij', D, np.eye(k)) + gj * P
        Vi = np.linalg.inv(V)
        _, ld = np.linalg.slogdet(V)
        a = Vi.sum((1, 2))
        Viy = Vi @ yb
        Viy1 = Vi.sum(1)
        qb = (yb * Viy).sum(-1) - (yb * Viy1).sum(-1) ** 2 / a
        out[:, j] = -0.5 * ld - 0.5 * np.log(a) - 0.5 * (N - 1) * np.log(ssw + qb)
    return out

def grid_u(n, ss, npts, width=9.0):
    k = len(n)
    c = 0.5 * np.log(ss / (n - 1))
    sd = np.sqrt(1 / (2 * (n - 1)))
    axes = []
    for i in range(k - 1):
        h = width * np.sqrt(sd[i] ** 2 + sd[-1] ** 2)
        axes.append(np.linspace(c[i] - c[-1] - h, c[i] - c[-1] + h, npts))
    mesh = np.meshgrid(*axes, indexing='ij')
    U = np.stack([m.ravel() for m in mesh], -1)
    dV = np.prod([ax[1] - ax[0] for ax in axes])
    return U, dV

def bf_quad(groups, r=0.5, kinds=("arith",), npts=41):
    n, yb, ss = stats(groups)
    k = len(n)
    U, dV = grid_u(n, ss, npts)
    s2 = np.exp(2 * np.c_[U, np.zeros(len(U))])  # sigma_k = 1 ; invariance d'echelle
    # log m0(u) : prod rho^-(n_i-1) |D|^-1/2 W^-1/2 Q0^-(N-1)/2 (invariant par c*rho)
    g0 = logm_given(s2, n, yb, ss, np.array([0.0]))[:, 0]
    lw = -0.5 * ((n - 1) * np.log(s2)).sum(-1) + g0
    M = lw.max()
    w0 = np.exp(lw - M)
    g, wg = g_weights(r)
    res = {}
    for kind in kinds:
        t2 = taus2(s2, n, kind)
        s2r = s2 / t2[:, None]
        # meme m0 quelle que soit la normalisation : on calcule BF(u,g) a rho normalise
        b0 = logm_given(s2r, n, yb, ss, np.array([0.0]))[:, 0]
        b1 = logm_given(s2r, n, yb, ss, g)
        bfu = (np.exp(b1 - b0[:, None]) * wg).sum(1)
        res[kind] = (w0 * bfu).sum() / w0.sum()
    return res

def bf_mc(groups, r=0.5, kind="arith", M=20000, rng=None):
    """E_{sigma|y,H0}[BF(sigma)], tirages exacts : mu sur grille, sigma_i^2|mu inverse-chi2."""
    rng = rng or np.random.default_rng(1)
    n, yb, ss = stats(groups)
    k = len(n)
    lo, hi = yb.min() - 10 * np.sqrt(ss.sum() / n.sum()), yb.max() + 10 * np.sqrt(ss.sum() / n.sum())
    mu = np.linspace(lo, hi, 200001)
    lp = sum(-(n[i] / 2) * np.log(ss[i] + n[i] * (yb[i] - mu) ** 2) for i in range(k))
    p = np.exp(lp - lp.max()); p /= p.sum()
    ms = rng.choice(mu, M, p=p)
    s2 = np.empty((M, k))
    for i in range(k):
        s2[:, i] = (ss[i] + n[i] * (yb[i] - ms) ** 2) / rng.chisquare(n[i], M)
    t2 = taus2(s2, n, kind)
    s2r = s2 / t2[:, None]
    g, wg = g_weights(r)
    g = g[::2]; wg = wg[::2] * 2
    # BF sachant sigma (tau connu) : forme gaussienne en tau^2=1 apres normalisation
    D = s2r / n
    P = np.eye(k) - np.ones((k, k)) / k
    Di = 1 / D
    a0 = Di.sum(1)
    Viy0 = Di * yb
    q0 = (yb * Viy0).sum(1) - Viy0.sum(1) ** 2 / a0
    ld0 = np.log(D).sum(1)
    tot = np.zeros(M)
    for gj, wj in zip(g, wg):
        V = np.einsum('mi,ij->mij', D, np.eye(k)) + gj * P
        Vi = np.linalg.inv(V)
        _, ld = np.linalg.slogdet(V)
        a = Vi.sum((1, 2))
        q = (yb * (Vi @ yb)).sum(1) - (yb * Vi.sum(1)).sum(1) ** 2 / a
        # expo en tau reel : q/tau^2 avec tau^2 = t2
        tot += wj * np.exp(-0.5 * (ld - ld0) - 0.5 * (np.log(a) - np.log(a0)) - 0.5 * (q - q0) / t2)
    return tot.mean(), tot.std() / np.sqrt(M), (tot.sum() ** 2) / (tot ** 2).sum()

if __name__ == "__main__":
    rng = np.random.default_rng(20261106)
    t0 = time.time()
    print("== Controle 1 : quadrature contre MC independant")
    for lab, ns, mus, sds in [("k2 modere", (15, 25), (0, .6), (1, 2)),
                               ("k3 H0 het", (5, 5, 40), (0, 0, 0), (4, 4, 1)),
                               ("k3 effet", (10, 20, 30), (0, .5, 1), (1, 2, 1))]:
        G = [rng.normal(m, s, int(n)) for n, m, s in zip(ns, mus, sds)]
        for kind in ("arith", "geom"):
            pts = 121 if len(ns) == 2 else 61
            q = bf_quad(G, kinds=(kind,), npts=pts)[kind]
            q2 = bf_quad(G, kinds=(kind,), npts=pts // 2 + 1)[kind]
            m, se, ess = bf_mc(G, kind=kind, M=20000, rng=rng)
            print(f"{lab:10s} {kind:6s} quad={q:.5f} (grille/2 {q2:.5f}) MC={m:.5f} se={se:.5f} ESS={ess:.0f}", flush=True)

    print("== Controle 2 : k=2, delta a priori ~ Cauchy(sqrt2/2), tau arith")
    g = 1 / rng.chisquare(1, 400000) * 0.25  # IG(1/2,r^2/2) = r^2/chi2_1
    a = rng.normal(0, np.sqrt(g * 2))
    print("mediane |delta|", np.median(np.abs(a)), "attendu", np.sqrt(2) / 2)

    print("== Scenario H0 n=(5,5,40) sd=(4,4,1), N=300")
    kinds = ("arith", "geom", "harm", "nweight", "welch")
    N = 300
    Ps = {k: [] for k in kinds}
    for it in range(N):
        G = [rng.normal(0, s, n) for n, s in zip((5, 5, 40), (4, 4, 1))]
        r = bf_quad(G, kinds=kinds, npts=31)
        for k in kinds: Ps[k].append(r[k] / (1 + r[k]))
        if it == 19: print("  temps 20 jeux", round(time.time() - t0), "s", flush=True)
    for k in kinds:
        p = np.array(Ps[k])
        print(f"{k:8s} med={np.median(p):.3f} P>.5={np.mean(p>.5):.3f} P>.95={np.mean(p>.95):.3f} P<.25={np.mean(p<.25):.3f}", flush=True)
    print("== Scenario H0 inverse n=(5,5,40) sd=(1,1,4), N=150")
    Ps = {k: [] for k in kinds}
    for it in range(150):
        G = [rng.normal(0, s, n) for n, s in zip((5, 5, 40), (1, 1, 4))]
        r = bf_quad(G, kinds=kinds, npts=31)
        for k in kinds: Ps[k].append(r[k] / (1 + r[k]))
    for k in kinds:
        p = np.array(Ps[k])
        print(f"{k:8s} med={np.median(p):.3f} P>.5={np.mean(p>.5):.3f} P>.95={np.mean(p>.95):.3f} P<.25={np.mean(p<.25):.3f}", flush=True)

    print("== Groupe tres dispersé : ecart 1 sigma entre groupes 1 et 2 (n=20), sd3 = S ; un jeu fixe puis S varie")
    base = [rng.normal(0, 1, 20), rng.normal(1, 1, 20)]
    z = rng.normal(0, 1, 20)
    for S in (1, 3, 10, 30, 100):
        G = base + [S * z]
        r = bf_quad(G, kinds=kinds, npts=41)
        print("S=%3d " % S + " ".join(f"{k}={r[k]:.3g}" for k in kinds), flush=True)
    print("fin", round(time.time() - t0), "s")
