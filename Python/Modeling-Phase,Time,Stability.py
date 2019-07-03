# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 13:18:18 2019

@author: anali
"""

import matplotlib.pyplot as plt
import Function as f
    
# Parameters
alpha = 10   #P_pTet = alpha, strength of the pTet promoter
beta = 1     #I_TetR = beta, Inhibition by TetR 
gamma = 2    #M_(C_2) = gamma, Dimerization of C
delta = 4    #P_Pe = delta, strength of the Pe promoter
epsilon = 3  #P_pBAD = epsilon, strength of the pBAD promoter
zeta = 1.5   #I_(C_2) = zeta, Inhibition by the C2 dimer
eta =1       #T_(Cox_4) = eta, Tetramization of Cox
theta =2     #I_Cox = theta, Inhibition by the Cox4 tetramer
iota =6      #I_Ara = iota, Inducment from Arabinose
DC = 0.2     #D_C = DC, Degredation rate of C
DCox = 0.2   #D_Cox = DCox, Degredation rate of Cox
DTet = 0.2   #D_TetR = DTet, Degredation of TetR
Ara = 0      #Arabinose concentration

# Initial concentrations
C = [[0.5],[1.5],[5],[10]]
Cox = [[0.5],[1.5],[5],[10]]
Tet = [[0.5],[1.5],[5],[10]]
#[[0],[0],[0],[0]]

# For the loops
dt = float(0.01) ; t1 = float(0.00) ; t2 = float(40.00) ; time = []

# Calculating f(x) and storing it in a container
fC =    [[f.functionforC(alpha,beta,Tet[0][0],DC, C[0][0],gamma)],  [f.functionforC(alpha,beta,Tet[1][0],DC, C[1][0],gamma)],[f.functionforC(alpha,beta,Tet[2][0],DC, C[2][0],gamma)],[f.functionforC(alpha,beta,Tet[3][0],DC, C[3][0],gamma)]] #Creates an empty list for values of f(C)
fCox =  [[f.functionforCox(delta, zeta, gamma, C[0][0], epsilon,DCox,Cox[0][0],eta,theta)],[f.functionforCox(delta, zeta, gamma, C[1][0], epsilon,DCox,Cox[1][0],eta,theta)],[f.functionforCox(delta, zeta, gamma, C[2][0], epsilon,DCox,Cox[2][0],eta,theta)],[f.functionforCox(delta, zeta, gamma, C[3][0], epsilon,DCox,Cox[3][0],eta,theta)]]  #Creates an empty list for values of f(Cox)
fTet =  [[f.functionforTetR(epsilon,iota,Ara,DTet,Tet[0][0])],[f.functionforTetR(epsilon,iota,Ara,DTet,Tet[1][0])],[f.functionforTetR(epsilon,iota,Ara,DTet,Tet[2][0])],[f.functionforTetR(epsilon,iota,Ara,DTet,Tet[3][0])]]  #Creates an empty list for values of f(TetR)

for i in range(0,4):
    g = 0
    time.append(0)
    for t in f.time_step(t1 ,t2,dt):
        C[i].append(C[i][g] + dt*f.functionforC(alpha,beta,Tet[i][g],DC, C[i][g],gamma))
        Cox[i].append(Cox[i][g]+dt*f.functionforCox(delta, zeta, gamma, C[i][g], epsilon,DCox,Cox[i][g],eta,theta))
        Tet[i].append(Tet[i][g]+dt*f.functionforTetR(epsilon,iota,Ara,DTet,Tet[i][g]))
        g += 1
        
        fC[i].append(f.functionforC(alpha,beta,Tet[i][g],DC, C[i][g],gamma)) #...the corresponding value of f(x) and save it with the same index by appending it
        fCox[i].append(f.functionforCox(delta, zeta, gamma, C[i][g], epsilon,DCox,Cox[i][g],eta,theta))
        fTet[i].append(f.functionforTetR(epsilon,iota,Ara,DTet,Tet[i][g]))
        
        time.append(t+dt)
        

time2 =[]
for i in time:
    time2.append(i)
    if i > t2:
        break
'''
#Time series plot
fig, axs = plt.subplots(2,2)
fig.suptitle('Time series')

axs[0,0].plot(time2,C[0],'b',time2, Cox[0], 'g',time2, Tet[0], 'm')
axs[0,0].set_title('Initial concentration: 0.5')



axs[0,1].plot(time2,C[1],'b',time2, Cox[1], 'g',time2, Tet[1], 'm')
axs[0,1].set_title('Initial concentration: 1.5')


axs[1,0].plot(time2,C[2],'b',time2, Cox[2], 'g',time2, Tet[2], 'm')
axs[1,0].set_title('Initial concentration: 5')


axs[1,1].plot(time2,C[3],'b',time2, Cox[3], 'g',time2, Tet[3], 'm')
axs[1,1].set_title('Initial concentration: 10')

hand = ['[C]','[Cox]',"[TetR]"]
for ax in axs.flat:
    ax.set(xlabel='Time', ylabel='Concentration')
    ax.legend(hand, title='Molecules')

'''

#Phase portrait plot
plt.plot(C[:][0],Cox[:][0],'k',C[:][1],Cox[:][1],'b',C[:][2],Cox[:][2],'g',C[:][3],Cox[:][3],'m')
hand = ['0.5','1.5','5','10']
plt.legend(hand, title='initial conce.')
plt.xlabel('[C]')
plt.ylabel('[Cox]')
plt.title('Phase portrait for C and Cox')
'''
#Stability plot
plt.plot(C[1],fC[1], 'b', Cox[1], fCox[1], 'g', Tet[1], fTet[1], 'm')
plt.xlabel('[Concentration]')
plt.ylabel('f(x)')
plt.title('Stability')
plt.grid()
hand = ['[C]','[Cox]',"[TetR]"]
plt.legend(hand, title='Molecules')
plt.show()
'''
