# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 15:05:05 2020

@author: rcard
"""

#As próximas 3 linhas são para selecionar entre plot inline ou em nova janela
#Útil para rlocus
from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'qt')
#get_ipython().run_line_magic('matplotlib', 'inline')

import numpy as np                #Biblioteca para cálculo numérico
import math                       #Funções matemáticas
import matplotlib.pyplot as plt   #Funções de plot similares as do MATLAB
import control as ctrl            #Biblioteca para controle
from control.matlab import *      #Funções para controle similares as do MATLAB
import scipy.linalg as la         #Biblioteca de álgebra 
import scipy.signal as sps        #Biblioteca de processamento de sinais


#Declaração do sistema
R=1
L=1e-3
Cap=150e-6

A=np.array([[0, 1/Cap],[-1/L, -R/L]])
B=np.array([[0], [1/L]])
C=np.array([[1, 0]])
D=np.array([[0]])
sys=ss(A,B,C,D)

#Polos do sistema em malha aberta
polos=la.eigvals(A)
print(polos)

#Matriz de controlabilidade e cálculo do posto da matriz
M=ctrb(A,B)
print(np.linalg.matrix_rank(M))

#Cálculo do vetor de ganhos K
polos_desejados=[-2582, -2582]
K=acker(A, B, polos_desejados)
print(K)

#Cálculo dos polos do sistema em malha fechada
polosMF=la.eigvals(A-B@K)
print(polosMF)

#Sistema em malha fechada
Amf=A-B@K
Bmf=np.zeros((2,1))
Cmf=C
Dmf=D
sysMF=ss(Amf,Bmf,Cmf,Dmf)

#Condições iniciais
x0=[[10], [0]]

#Tempo de simulação
t=np.linspace(0,0.02,1000)

#Resposta às condições iniciais
YMA, tplotMA=initial(sys,t,x0)
YMF, tplotMF=initial(sysMF,t,x0)

plt.figure()
plt.subplot(2,1,1)
plt.plot(t,YMA)
plt.ylabel('$V_c(t)$')
plt.xlabel('t(s)')
plt.legend(['Malha aberta'])

plt.subplot(2,1,2)
plt.plot(t,YMF)
plt.ylabel('$V_c(t)$')
plt.xlabel('t(s)')
plt.legend(['Malha fechada'])
plt.show()


















