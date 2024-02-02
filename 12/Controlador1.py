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
A=np.array([[0,1, 0],[0, 0, 1],[0, -2, -3]])
B=np.array([[0], [0], [1]])
C=np.array([[1, 0, 0]])
D=np.array([[0]])
sys=ss(A,B,C,D)

#Matriz de controlabilidade
M=ctrb(A,B)
#Posto da matriz de controlabilidade
print(np.linalg.matrix_rank(M))

#Polos de MF desejados
polos_desejados=[-2+2*math.sqrt(3)*1j, -2-2*math.sqrt(3)*1j, -10]
#Projeto do vetor de ganhos K do controlador
K=acker(A, B, polos_desejados)
print(K)

#Sistema em MF com kr=1
Amfc=A-B@K
Bmfc=B
Cmfc=C
Dmfc=D
sysMFc=ss(Amfc,Bmfc,Cmfc,Dmfc)

#Cálculo de kr
GanhoCC=dcgain(sysMFc)
kr=1/GanhoCC
print(kr)

#Sistema em MF com kr ajustado
Amfca=A-B@K
Bmfca=B*kr
Cmfca=C
Dmfca=D
sysMFca=ss(Amfca,Bmfca,Cmfca,Dmfca)

#Sistema em MF sem controlador
Amf=A-B@C
Bmf=B
Cmf=C
Dmf=D
sysMF=ss(Amf,Bmf,Cmf,Dmf)

#Resposta ao degrau em MF
t=np.linspace(0,20,1000)
YMF, tplotMF=step(sysMF,t)
YMFca, tplotMFca=step(sysMFca,t)
YMFc, tplotMFc=step(sysMFc,t)

plt.figure()
plt.subplot(2,1,1)
plt.plot(t,YMF,t,YMFca)
plt.ylabel('$Amplitude$')
plt.xlabel('t(s)')
plt.legend(['Sistema sem controlador', 'Sistema com controlador (kr=160)'])

plt.subplot(2,1,2)
plt.plot(t,YMFc)
plt.ylabel('$Amplitude$')
plt.xlabel('t(s)')
plt.legend(['Sistema com controlador (kr=1)'])
plt.show()


















