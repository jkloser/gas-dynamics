import numpy as np
import sympy as sp
delta = 8*np.pi/180
M1 = 3.0
gamma = 1.4
#Beta
error_beta = 10
beta_max = 2.5
beta_min = 0.0

while abs(error_beta) > 0.0001:
    beta = (beta_max + beta_min)/2
    error_beta = ((M1**2 * (np.sin(beta))**2 - 1) * sp.cot(beta)) / (((gamma+1)/2) * M1**2 - M1**2 * (np.sin(beta))**2 + 1) - np.tan(delta)

    if error_beta > 0:
        beta_max = beta
    else:
        beta_min = beta

print(beta*180/np.pi)