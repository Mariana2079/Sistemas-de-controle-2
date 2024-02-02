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


#Sistemas descrito por função de transferência
num=[1, 3, 2]
den=[1, 12, 47, 60]
G=tf(num,den)
print(G)

#Descrição do sistema na forma controlável
A=np.array([[0, 1, 0],[0, 0, 1],[-60, -47, -12]])
B=np.array([[0], [0], [1]])
C=np.array([[2, 3, 1]])
D=np.array([[0]])

#Representação na forma controlável
sysC=ss(A,B,C,D)
print(sysC)

#Representação na forma observável
sysO=ss(np.transpose(A),np.transpose(C),np.transpose(B),D)
print(sysO)

#Matriz de transformação para diagonalização
P=np.array([[1, 1, 1],[-3, -4, -5],[9, 16, 25]])

#Obtenção de outra matriz P através do Matlab
#[aVet, aVal]=eig(A)
#P=aVet;

#Sistema diagonalizado
Ad1=np.linalg.inv(P)@A@P
Bd1=np.linalg.inv(P)@B
Cd1=C@P
Dd1=D;
sysD1=ss(Ad1,Bd1,Cd1,Dd1)
print(sysD1)

#Representação na forma diagonal a partir
#da função de transferência
[Res,Polos,K] = sps.residue(num,den) 
print(Res)
print(Polos)

Ad2=[[-3, 0, 0],[0, -4, 0],[0, 0, -5]]
Bd2=[[1],[1],[1]]
Cd2=[[1, -6, 6]]
Dd2=[[0]]
sysD2=ss(Ad2,Bd2,Cd2,Dd2)
print(sysD2)

#Comparação das respostas ao degrau para
#diferentes representações de um mesmo sistema
t=np.linspace(0,2,100)
YC, tC=step(sysC,t)
YO, tO=step(sysO,t)
YD1, tD1=step(sysD1,t)
YD2, tD2=step(sysD2,t)
YG, tG=step(G,t)
plt.figure()
plt.subplot(5,1,1)
plt.plot(tC,YC)
plt.legend(['Saída y - forma controlável'])
plt.ylabel('Amplitude')
plt.subplot(5,1,2)
plt.plot(tO,YO)
plt.legend(['Saída y - forma observável'])
plt.ylabel('Amplitude')
plt.subplot(5,1,3)
plt.plot(tD1,YD1)
plt.legend(['Saída y - forma diagonal 1'])
plt.ylabel('Amplitude')
plt.subplot(5,1,4)
plt.plot(tD2,YD2)
plt.legend(['Saída y - forma diagonal 2'])
plt.ylabel('Amplitude')
plt.subplot(5,1,5)
plt.plot(tG,YG)
plt.legend(['Saída y - G(s)'])
plt.xlabel('t(s)')
plt.ylabel('Amplitude')
plt.show()

#Representação na forma de Jordan a partir
#da função de transferência
numj=[1, 3, 2]
denj=[1, 11, 39, 45]
Gj=tf(numj,denj)
print(Gj)
[Rj,Pj,Kj] = sps.residue(numj,denj) 
print(Rj)
print(Pj)

#Sistema na forma de Jordan
Aj=[[-3, 1, 0],[0, -3, 0],[0, 0, -5]]
Bj=[[0],[1],[1]]
Cj=[[1, -2, 3]]
Dj=[0]
sysj=ss(Aj,Bj,Cj,Dj)
print(sysj)

#Comparação das respostas ao degrau do sistema representado
#por função de transferência e na forma de Jordan
YGj, tGj=step(Gj,t)
Yj, tj=step(sysj,t)
plt.figure()
plt.subplot(2,1,1)
plt.plot(tGj,YGj)
plt.legend(['Saída y - G(s)'])
plt.ylabel('Amplitude')
plt.subplot(2,1,2)
plt.plot(tj,Yj)
plt.legend(['Saída y - forma de Jordan'])
plt.xlabel('t(s)')
plt.ylabel('Amplitude')
plt.show()