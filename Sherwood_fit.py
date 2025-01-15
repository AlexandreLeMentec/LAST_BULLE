import numpy  as np
from scipy import stats as st
import statsmodels.api as sm
import matplotlib.pyplot as plt

data = ['data_alpha_0_69.csv', 'data_alpha_1_46.csv', 'data_alpha_2_55.csv', 'data_alpha_5_13.csv',
        'data_alpha_9_42.csv', 'data_alpha_15_16.csv', 'data_alpha_19_44.csv', 'data_alpha_24_16.csv',
        'data_alpha_29_58.csv']
Sherwood = []
alphaG = np.array([0.69,1.46,2.55,5.13,9.42,15.16,19.44,24.16,29.58])/100
alphaL = 1 - alphaG

############################# ETUDE SH #####
for i in range(len(data)):
    data = np.loadtxt(data[i], delimiter=',')
    t = data[:, 0] # t (s)
    Y = data[:, 1] # C (mol/vol)

    dt = t[1]-t[0]

    # on enlève le lag au début 
    n = 0
    while Y[n] < 1e-3:
        n += 1

    t = t[n-1:] # on tronque jusqu'à une variation satisfaisante
    Cl = Y[n-1:]
    t = t - t[0]

    ############# values #############
    # Cte paramètres phy
    He = 4.34*10**9 # constante de Henry de l'ox dans l'eau
    Cg = 0.2095 # Concentration Ox dans l'air 
    Dl = 2.1*10**-9 # coeff diff

    # Cte paramètres exp
    d0 = 2.1*10**-3 # min size of bulles
    aI0 = 3011 
    vz0 = 0.32

    # var exp 
    alpha_exp = alphaL[i]
    alphaG_exp = alphaG[i]

    if alpha_exp <= 0.023 :
        d = 15*alpha_exp*d0 + d0
        alphaI = 0.402*alpha_exp**0.85
    else:
        d = 2.3*alpha_exp**0.5*d0 + d0
        alphaI = 0.336*alpha_exp**0.8
    
    def renormalize(Cl_l):
        for i in range(np.shape(Cl_l)[0]):
            if Cl_l[i] > 1 :
                Cl_l[i] = 0.99
        return Cl_l
        
    Y = -np.log(1-renormalize(Cl))
    print( 'Y calculé')
    n = 0
    # while Y[n] != np.nan and Y[n] != np.inf:
    #     n += 1
    
    def tau_fit(y,x):
        model_y1 = sm.RLM(y,x) 
        # on utilise un algorithme de régression linéaire "Robust Linear Model"
        result_y1 = model_y1.fit()
        return result_y1.params

    tau = 1 / tau_fit(Y[0:n],t[0:n])[0]

    Sherwood.append(alpha_exp/alphaG_exp*d/Dl*1/tau)

print(Sherwood)
input("Cliquez pour couper")