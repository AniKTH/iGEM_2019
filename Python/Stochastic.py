# Importing packages #

import matplotlib.pyplot as plt
import numpy as np
import Function as f

# EXPLICIT EULER SCHEME #

# Parameters
alpha = 1   #P_pTet = alpha, strength of the pTet promoter
beta = 1     #I_TetR = beta, Inhibition by TetR 
gamma = 1    #M_(C_2) = gamma, Dimerization of C
delta = 1    #P_Pe = delta, strength of the Pe promoter
epsilon = 0 #P_pBAD = epsilon, strength of the pBAD promoter
zeta = 1   #I_(C_2) = zeta, Inhibition by the C2 dimer
eta = 1       #T_(Cox_4) = eta, Tetramization of Cox
theta = 1     #I_Cox = theta, Inhibition by the Cox4 tetramer
iota = 0      #I_Ara = iota, Inducment from Arabinose
DC = 1     #D_C = DC, Degredation rate of C
DCox = 1   #D_Cox = DCox, Degredation rate of Cox
DTet = 0   #D_TetR = DTet, Degredation of TetR
Ara = 0      #Arabinose concentration

dt = 0.01
C = 0.5
t1 = 0.00
t2 = 10
Cox = 0.5
Tet = 0
tau = 0.01


# Creates containers for the different values the concentrations and time will take
# Updated by the tau-leap method #

time = []
C_list = []
Cox_list = []
Tet_list = []

#For loopfor the tau-loop method
for i in f.time_step(t1 ,t2,dt):
    #Tau-leap method for the C protien
    #f = alpha - beta*Tet-DC*C-gamma*C
    fC= f.functionforC(alpha,beta,Tet,DC, C,gamma)
    n1C = np.random.poisson(abs(alpha)*tau,1)
    n2C = np.random.poisson(beta*Tet*tau,1)
    n3C = np.random.poisson(abs(DC*C)*tau,1)
    n4C = np.random.poisson(abs(gamma*C)*tau,1)
    fC = fC + n1C - n2C - n3C -n4C #Add the positive impact and remove the negative
    
    #Tau-leap method for the Cox protein
    #f = delta - zeta*gamma*C + epsilon-DCox*Cox-eta*Cox-theta*Cox
    fCox = f.functionforCox(delta, zeta, gamma, C, epsilon,DCox,Cox,eta,theta)
    n1Cox = np.random.poisson(delta*tau,1)
    n2Cox = np.random.poisson(abs(zeta*gamma*C*Tet)*tau,1)
    n3Cox = np.random.poisson(epsilon*tau,1)
    n4Cox = np.random.poisson(abs(DCox*Cox*C)*tau,1)
    n5Cox = np.random.poisson(eta*Cox*tau, 1)
    n6Cox = np.random.poisson(theta*Cox*tau, 1)
    fCox = fCox + n1Cox - n2Cox + n3Cox - n4Cox - n5Cox - n6Cox
    
    #Tau-leap method for the TetR protein
    #f = epsilon + iota*Ara-DTet*Tet    
    fTet = f.functionforTetR(epsilon,iota,Ara,DTet,Tet)
    #f = epsilon + iota*Ara-DTet*Tet
    n1Tet = np.random.poisson(epsilon*tau,1)
    n2Tet = np.random.poisson(iota*Ara*tau,1)
    n3Tet = np.random.poisson(DTet*Tet*tau,1)
    fTet = fTet + n1Tet + n2Tet - n3Tet

    #Updating and storing values
    C = C + dt*fC
    Cox = Cox + dt*fCox
    Tet = Tet + dt*fTet
    t = i + dt
    time.append(t)
    C_list.append(C)
    Cox_list.append(Cox)
    Tet_list.append(Tet)


# Plot c, cox and tetR protein dynamics over time #
fig = plt.figure()
#plt.axis([-1,10.05,-1,6.5])
plt.scatter (time,C_list, s=1)
plt.scatter (time,Cox_list, s=1)
plt.scatter (time,Tet_list, s=1)
hand = ['[C]','[Cox]',"[TetR]"]
plt.legend(hand, title='Protein')
plt.title ("C, Cox and TetR protein dynamics over time (AB)")
plt.xlabel ("Time")
plt.ylabel ("Protein concentration")
plt.show()