## Gas Dynamics Project 1\\John Kloser & Maxim Strehle\\03/06/2020
##Prandtl-Meyer Flow
##
import numpy as np
import sympy as sp
import matplotlib.pyplot as plot
import sys

M1 = 3.0
P1 = 8.7023 #psia
P3 = 10.1527 #psia
gamma = 1.4

P01 = (1/0.02722) * P1
max_del = 0.9*(M1 - 1)**(3/2)
del_0 = 29.7*np.pi/180
for delta in np.arange(del_0, 45*np.pi/180, 0.000001):

    #Beta
    error_beta = 10
    beta_max = np.pi/2
    beta_min = 0.0
    while abs(error_beta) > 0.01:
        beta = (beta_max + beta_min)/2
        error_beta = ((M1**2 * (np.sin(beta))**2 - 1) * sp.cot(beta)) / (((gamma+1)/2) * M1**2 - M1**2 * (np.sin(beta))**2 + 1) - np.tan(delta)

        if error_beta > 0:
            beta_max = beta
        else:
            beta_min = beta
    print('delta = ' + str(delta*180/np.pi))

    #Normal Mach
    Mn1 = M1 * np.sin(beta)
    Mn2 = np.sqrt((2+(gamma-1)*Mn1**2) / (2*gamma*Mn1**2 - (gamma - 1)))
    M2 = Mn2/np.sin(beta - delta)
    P02_ratio = (Mn1/Mn2) * (1+(gamma-1)*Mn2**2/2)**((gamma+1)/(2*gamma-2)) / ((1+(gamma-1)*Mn1**2/2)**((gamma+1)/(2*gamma-2)))
    P02 = P02_ratio*P01

    #Prandtl-Meyer Expansion
    #delta(theta) = delta or turn angle
    gamma_ratio = (gamma+1)/(gamma-1)
    nu2 = np.sqrt(gamma_ratio) * np.arctan(np.sqrt((1/gamma_ratio) * (M2**2 - 1))) - np.arctan(M2**2 - 1)
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
        df = gamma_ratio * np.sqrt(gamma_ratio)*M3 / (M3**2 - 1 + gamma_ratio * np.sqrt(M3**2 - 1)) - (1/(M3*np.sqrt(M3**2 -1)))
        M3 = M3 - f/df
    P3_ratio = (1+(gamma-1)*M3**2/2)**(gamma/(gamma-1)) #P03/P3
    P3_guess = (1/P3_ratio)*P02 #P02 = P03
    if abs(P3-P3_guess) < 0.00001:
        print('Pressures are equal motherfucker')
        break
    else:
        print('Pressures not equal')
print('M2 = ' + str(M2))
print ('M3 = ' + str(M3))
print('delta = ' + str(delta*180/np.pi))
print('beta = ' + str(beta*180/np.pi))
#p2/p1
#p3/p2
print('P3 = ' + str(P3_guess))
print(P3_ratio)