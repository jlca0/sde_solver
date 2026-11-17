import numpy as np
import matplotlib.pyplot as plt
import sde

X0, mu0, sig0, T = 1.0, 0.5, 0.4, 1.0   # ajusta a los parámetros del libro
kmax = 14
rng = np.random.default_rng(1)

mu      = lambda t, x: mu0 * x
sigma   = lambda t, x: sig0 * x
sigma_x = lambda t, x: sig0


def refinar(camino, T, rng):
    W = np.asarray(camino)
    h = T / (len(W) - 1)
    nuevo = np.empty(2 * len(W) - 1)
    nuevo[::2] = W
    nuevo[1::2] = 0.5 * (W[:-1] + W[1:]) + 0.5 * np.sqrt(h) * rng.standard_normal(len(W) - 1)
    return nuevo


# camino inicial con 2 pasos (log2 N = 1)
W = sde.BrownianPath(T, 2)
XX = sde.ExactGeomPath(X0, T, mu0, sig0, W).get_trayectoria()[-1]  # solo depende de W_T

ks, XE, XM = [], [], []
for k in range(1, kmax + 1):
    n = len(W.camino) - 1                    # n = 2^k
    em = sde.EulerMaruyama(n);  em.solve(X0, T, mu, sigma, W)
    ms = sde.Milstein(n);       ms.solve(X0, T, mu, sigma, sigma_x, W)
    ks.append(k)
    XE.append(em.get_trayectoria()[-1])
    XM.append(ms.get_trayectoria()[-1])
    if k < kmax:
        W.camino = refinar(W.camino, T, rng).tolist()

plt.figure(figsize=(7, 5))
plt.plot(ks, XE, "-",  color="k", label="Euler")
plt.plot(ks, XM, "--", color="k", label="Milstein")
plt.axhline(XX, ls=":", color="k", label="valor exacto")
plt.xlabel(r"$\log_2(N)$"); plt.ylabel("S(T)")
plt.title("Milstein vs Euler"); plt.legend()

plt.show()