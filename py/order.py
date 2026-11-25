import numpy as np
import matplotlib.pyplot as plt

import sde

X0 = 100
T = 1
n_realizaciones = 10000
n_nodos = 1000

muestra_bp = [sde.BrownianPath(T, n_nodos) for _ in range(n_realizaciones)]

muestra_gbm = np.zeros((n_realizaciones,n_nodos+1))
for k in range(n_realizaciones):
    muestra_gbm[k,:] = sde.ExactGeomPath(X0, T, 1, 0.5, muestra_bp[k]).get_trayectoria()


print('---Orden de Euler-Maruyama--- \n')
malla = [125,250,500,1000]
mu      = lambda t, x:  x
sigma   = lambda t, x: 0.5 * x
error0 = 1
for N in malla:  # Refinamiento
    error1 = 0
    for k in range(n_realizaciones):  # Montecarlo
        aprox = sde.EulerMaruyama(N)
        aprox.solve(X0, T, mu, sigma, muestra_bp[k])
        aprox_array = np.array(aprox.get_trayectoria())
        error1 += abs(aprox_array[-1]-muestra_gbm[k,-1])
    error1 = error1/n_realizaciones
    if N != malla[0]:
        print(f'El error es: {np.log(error0/error1)/np.log(2):.4f}')
    else:
        print('El error es: ---')
    error0 = error1
print('\n')   

print('---Orden de Milstein--- \n')
malla = [125,250,500,1000]
mu      = lambda t, x:  x
sigma   = lambda t, x: 0.5 * x
sigmax   = lambda t, x: 0.5
error0 = 1
for N in malla:  # Refinamiento 
    error1 = 0
    for k in range(n_realizaciones):  # Montecarlo
        aprox = sde.Milstein(N)
        aprox.solve(X0, T, mu, sigma, sigmax, muestra_bp[k])
        aprox_array = np.array(aprox.get_trayectoria())
        error1 += abs(aprox_array[-1]-muestra_gbm[k,-1])
    error1 = error1/n_realizaciones
    if N != malla[0]:
        print(f'El error es: {np.log(error0/error1)/np.log(2):.4f}')
    else:
        print('El error es: ---')
    error0 = error1
print('\n')    