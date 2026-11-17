import numpy as np
import matplotlib.pyplot as plt
import sde

# Parámetros del modelo
X0 = 100.0
mu0, sig0 = 0.05, 0.2
mu = lambda t, x: mu0*x
sigma = lambda t, x: sig0*x
sigma_x = lambda t, x: sig0
T = 10.0
mesh = [100, 200, 400, 800, 1600]
incs = [16, 8, 4, 2, 1]

# ---------- PARTE 1: ilustración con una trayectoria ----------
fig = plt.figure(figsize=(12, 6))
bp = sde.BrownianPath(T, 1600)
exact_path = sde.ExactGeomPath(X0, T, mu0, sig0, bp)
exact = np.array(exact_path.get_trayectoria())

fig1 = fig.add_subplot(121)
fig1.plot(np.linspace(0, T, 1601), exact, label="Exacta", color="black")
for N in mesh:
    em = sde.EulerMaruyama(N)
    em.solve(X0, T, mu, sigma, bp)
    fig1.plot(np.linspace(0, T, N + 1), em.get_trayectoria(),
              label=f"Euler-Maruyama N={N}", alpha=0.5)
fig1.set(xlabel="t", ylabel="X(t)", title="Euler-Maruyama")
fig1.legend(fontsize=8); fig1.grid(True)

fig2 = fig.add_subplot(122)
fig2.plot(np.linspace(0, T, 1601), exact, label="Exacta", color="black")
for N in mesh:
    mst = sde.Milstein(N)
    mst.solve(X0, T, mu, sigma, sigma_x, bp)
    fig2.plot(np.linspace(0, T, N + 1), mst.get_trayectoria(),
              label=f"Milstein N={N}", alpha=0.5)
fig2.set(xlabel="t", ylabel="X(t)", title="Milstein")
fig2.legend(fontsize=8); fig2.grid(True)
plt.tight_layout()

# ---------- PARTE 2: orden de convergencia fuerte (M realizaciones) ----------
M = 500
errE = np.zeros(len(mesh))
errM = np.zeros(len(mesh))

for _ in range(M):
    bp = sde.BrownianPath(T, 1600)
    XX = sde.ExactGeomPath(X0, T, mu0, sig0, bp).get_trayectoria()[-1]
    for j, N in enumerate(mesh):
        em = sde.EulerMaruyama(N);  em.solve(X0, T, mu, sigma, bp)
        ms = sde.Milstein(N);       ms.solve(X0, T, mu, sigma, sigma_x, bp)
        errE[j] += abs(em.get_trayectoria()[-1] - XX)
        errM[j] += abs(ms.get_trayectoria()[-1] - XX)

errE /= M
errM /= M

print("Orden Euler:   ", np.log2(errE[:-1] / errE[1:]))
print("Orden Milstein:", np.log2(errM[:-1] / errM[1:]))

plt.show()