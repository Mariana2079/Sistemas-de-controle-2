# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 20:38:39 2020

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
A=np.array([[0,1, 0],[0, 0, 1],[-6, -11, -6]])
B=np.array([[0], [0], [1]])
C=np.array([[1, 0, 0]])
D=np.array([[0]])
sys=ss(A,B,C,D)

#Polos da planta
print(la.eigvals(A))

#Sistema aumentado
Abar=np.concatenate((A, -C))
Abar=np.concatenate((Abar, np.zeros((4,1))), axis=1)
Bbar=np.concatenate((B, [[0]]))
Cbar=np.concatenate((C, [[0]]), axis=1)
Dbar=D

#Matriz de controlabilidade
M=ctrb(Abar,Bbar)
#Posto da matriz de controlabilidade
print(np.linalg.matrix_rank(M))

#Polos de MF desejados
polosMF=[-2+2*math.sqrt(3)*1j, -2-2*math.sqrt(3)*1j, -10, -50]
#Projeto do vetor de ganhos K do controlador
Kbar=acker(Abar, Bbar, polosMF)
print(Kbar)

#Sistema em MF sem controlador
Amf=A-B@C
Bmf=B
Cmf=C
Dmf=D
sysMF=ss(Amf,Bmf,Cmf,Dmf)

#Polos do sistema em MF sem controlador
print(la.eigvals(Amf))

#Sistema em MF com controlador
Amfc=Abar-Bbar@Kbar
Bmfc=np.array([[0], [0], [0], [1]])
Cmfc=Cbar
Dmfc=Dbar
sysMFc=ss(Amfc,Bmfc,Cmfc,Dmfc)

#Polos do sistema em MF com controlador
print(la.eigvals(Amfc))


#Resposta ao degrau em MF
t=np.linspace(0,20,1000)
YMF, tplotMF=step(sysMF,t)
YMFc, tplotMFc=step(sysMFc,t)

plt.figure()
plt.plot(t,YMF,t,YMFc)
plt.ylabel('$Amplitude$')
plt.xlabel('t(s)')
plt.legend(['Sistema sem controlador', 'Sistema com controlador'])
plt.show()



##########################################
# Análise das variações paramétricas
a=np.array([[0,1, 0],[0, 0, 1],[-60, -11, -6]]) #Variação em A
abar=np.concatenate((a, -C))
abar=np.concatenate((abar, np.zeros((4,1))), axis=1)

#Sistema controlado em MF com variação em A
Amfp1=abar-Bbar@Kbar
Bmfp1=np.array([[0], [0], [0], [1]])
Cmfp1=Cbar
Dmfp1=Dbar
sysMFp1=ss(Amfp1,Bmfp1,Cmfp1,Dmfp1)
G1=ss2tf(sysMFp1)
print(G1)
print(pole(G1))
print(dcgain(G1))


c=np.array([[0.5, 0, 0]]) #Variação em C
aba2r=np.concatenate((A, -c))
abar2=np.concatenate((aba2r, np.zeros((4,1))), axis=1)
cbar=np.concatenate((c, [[0]]), axis=1)

#Sistema controlado em MF com variação em C
Amfp2=abar2-Bbar@Kbar
Bmfp2=np.array([[0], [0], [0], [1]])
Cmfp2=cbar
Dmfp2=Dbar
sysMFp2=ss(Amfp2,Bmfp2,Cmfp2,Dmfp2)
G2=ss2tf(sysMFp2)
print(G2)
print(pole(G2))
print(dcgain(G2))

#Resposta ao degrau em MF para os sistemas com variações paramétricas
YMFp1, tplotMFp1=step(sysMFp1,t)
YMFp2, tplotMFp2=step(sysMFp2,t)

plt.figure()
plt.plot(t,YMF,t,YMFc,t,YMFp1,t,YMFp2)
plt.legend(['Sem controlador','Com controlador e sem var. paramétrica','Com controlador e com var. paramétrica em A','Com controlador e com var. paramétrica em C'])
plt.title('Resposta ao degrau em malha fechada')
plt.ylabel('$Amplitude$')
plt.xlabel('t(s)')
plt.show()


