# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 12:54:08 2019

@author: anali
"""

#Function.py
#The rate of the C protein
def functionforC(alpha,beta,Tet,DC, C,gamma):
    f = alpha - beta*Tet-DC*C-gamma*C
    return f

                #d[C]/dt=P_pTet-I_TetR*[TetR]-D_C*[C]-M_(C_2)*[C]  
                #P_pTet = alpha, strength of promoter pTet
                #I_TetR = beta, Inhibition by TetR
                #M_(C_2) = gamma, Dimerization of C
                #D_C = DC, Degradation of C
                        
#The rate of the Cox protien
def functionforCox(delta, zeta, gamma, C, epsilon,DCox,Cox,eta,theta):
    f = delta - zeta*gamma*C + epsilon-DCox*Cox-eta*Cox-theta*Cox
    return f

                #d[Cox]/dt=P_Pe+P_pBAD-D_Cox*[Cox]-I_(C_2)*M_(C_2)*[C]-T_(Cox_4)*[Cox]-I_Cox*T_(Cox_4)*[Cox]
                #P_Pe = delta, strength of the Pe promoter
                #P_pBAD = epsilon, strength of the pBAD promoter
                #I_(C_2) = zeta, Inhibition by the C2 dimer
                #M_(C_2) = gamma, Dimerization of C
                #T_(Cox_4) = eta, tetramization of Cox
                #I_Cox = theta, Inhibition by the Cox4 tetramer
                #D_Cox = DCox, Degredation of Cox

#The rate of the TetR protein
def functionforTetR(epsilon,iota,Ara,DTet,Tet):
    f = epsilon + iota*Ara-DTet*Tet
    return f

                #d[TetR]/dt=P_pBAD+I_Ara*[Arabinose]-D_TetR*[TetR]
                #P_pBAD = epsilon, strength of the pBAD promoter
                #I_Ara = iota, Inducement from Arabinose
                #D_TetR = DTet, degredation of TetR
                        
# For using small time steps
def time_step(start, stop=None, step=None):
    if stop == None:
        stop = start + 0.0
        start = 0.0
    if step == None:
        step = 1.0
    while start < stop:
        yield start
        start += step
