import numpy as np
import sympy as sp
delta = 8*np.pi/180
M2 = 3.0
gamma = 1.4
gamma_ratio = (gamma+1)/(gamma-1)
nu2 = np.sqrt(gamma_ratio) * np.arctan(np.sqrt((1/gamma_ratio) * (M2**2 - 1))) - np.arctan(np.sqrt(M2**2 - 1))
nu3 = delta + nu2

#Hall Inverse Approximation
A = 1.3604
B = 0.0962
C = -0.5127
D = -0.6722
E = -0.3278
nu_inf = np.pi/2 * (np.sqrt(6) - 1)
y = (nu3/nu_inf)**(2/3)
M3 = (1 + A*y + B*y**2 + C*y**3)/(1 + D*y + E*y**2)

#Newton Iterations
f = 10
while abs(f) > 0.001:
    f = np.sqrt(gamma_ratio) * np.arctan(np.sqrt((1/gamma_ratio) * (M3**2 - 1))) - np.arctan(np.sqrt(M3**2 - 1)) - nu3
    #df = beta*(1-gamma**2)/(M3*(1+gamma**2*beta**2))
    df = gamma_ratio * np.sqrt(gamma_ratio)*M3 / (M3**2 - 1 + gamma_ratio * np.sqrt(M3**2 - 1)) - (1/(M3*np.sqrt(M3**2 -1)))
    M3 = M3 - f/df
print(delta)
#P3_ratio = (1+(gamma-1)*M3**2/2)**(gamma/(gamma-1)) #P03/P3
#P3_guess = (1/P3_ratio)*P02 #P02 = P03
#if abs(P3-P3_guess) < 0.01:
#    print(M3)
#    print(delta)