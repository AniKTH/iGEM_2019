# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 13:29:49 2019

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
#DCox = 0.2   #D_Cox = DCox, Degredation rate of Cox
DTet = 0.2   #D_TetR = DTet, Degredation of TetR
Ara = 0      #Arabinose concentration

#Altered parameter
p1 = float(0) ; p2 = float(14) ; dp = float(0.001) #Lower limit ; Upper limit ; Stepwise change
DCox = [] # Creates an empty list for the different values of the aptered parameter (DTet)

for i in f.time_step(p1,p2,dp):
    DCox.append(i) # The list of values of parameter

# Creating the list of concentration values for C
C1 = float(0) ; C2 = float(10) ; dC = float(0.001) #Lower limit ; Upper limit ; Stepwise change
C = [] # Creates an empty list for the different values of C

for i in f.time_step(C1,C2,dC):
    C.append(i) # The list of concentrations of C
    
# Creating the list of concentration values for C
Cox1 = float(0) ; Cox2 = float(10) ; dCox = float(0.001) #Lower limit ; Upper limit ; Stepwise change
Cox = [] # Creates an empty list for the different values of Cox

for i in f.time_step(Cox1,Cox2,dCox):
    Cox.append(i) # The list of concentrations of Cox

# Creating the list of concentration values for TetR
Tet1 = float(0) ; Tet2 = float(10) ; dTet = float(0.001) #Lower limit ; Upper limit ; Stepwise change
Tet = [] # Creates an empty list for the different values of TetR

for i in f.time_step(Tet1,Tet2,dTet):
    Tet.append(i) # The list of concentrations of TetR
    
# Creating containers
Cy = [] #Container for the steady states
Cx = [] #Container for the cooresponding values of parameter
Coxy = [] #Container for the steady states
Coxx = [] #Container for the cooresponding values of parameter
Tety = [] #Container for the steady states
Tetx = [] #Container for the cooresponding values of parameter

for p in DCox: # For every value of gamma...
    for i in range(0,len(C)): #... for every value of c....
        fC = f.functionforC(alpha,beta,Tet[i], DC, C[i],gamma) #... calculate f(x)...
        fCox = f.functionforCox(delta, zeta, gamma, C[i], epsilon,p,Cox[i],eta,theta)
        fTet = f.functionforTetR(epsilon,iota,Ara,DTet,Tet[i])
        if abs(fC) < 0.001: #... if f(x) is cloes enough to 0...
            Cy.append(C[i]) #... append the used value of C to the list y...
            Cx.append(p) #... and append the used value of gamma to the list x. 
            
        if abs(fCox) < 0.001: #... if f(x) is cloes enough to 0...
            Coxy.append(Cox[i]) #... append the used value of C to the list y...
            Coxx.append(p) #... and append the used value of gamma to the list x. 
        
        if abs(fTet) < 0.001: #... if f(x) is cloes enough to 0...
            Tety.append(Tet[i]) #... append the used value of C to the list y...
            Tetx.append(p) #... and append the used value of gamma to the list x. 

    print(p) #For Analis peace of mind
    
plt.plot(Cx,Cy, 'k.', Coxx, Coxy, 'b.', Tetx, Tety, 'g.') # Need to plot for every point
plt.xlabel('DCox')
plt.ylabel('[Concentration]')
plt.title('Bifurcation diagram: DCox')
hand = ['[C]','[Cox]',"[TetR]"]
plt.legend(hand, title='Molecules')
plt.grid()
plt.show()