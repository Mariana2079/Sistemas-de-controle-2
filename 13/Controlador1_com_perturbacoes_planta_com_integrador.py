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
polosMF=[-2+2*math.sqrt(3)*1j, -2-2*math.sqrt(3)*1j, -10]
#Projeto do vetor de ganhos K do controlador
K=acker(A, B, polosMF)
print(K)

#Sistema controlado em MF com kr=1
Amfc=A-B@K
Bmfc=B
Cmfc=C
Dmfc=D
sysMFc=ss(Amfc,Bmfc,Cmfc,Dmfc)

#Cálculo de kr
GanhoCC=dcgain(sysMFc)
kr=1/GanhoCC
print(kr)

#Sistema controlado em MF com kr ajustado
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



##########################################
# Análise das variações paramétricas
a=np.array([[0,1, 0],[0, 0, 1],[0, -20, -3]]) #Variação em A
c=np.array([[0.5, 0, 0]])  #Variação em C

#Sistema controlado em MF variação em A
Amfp1=a-B@K
Bmfp1=B*kr
Cmfp1=C
Dmfp1=D
sysMFp1=ss(Amfp1,Bmfp1,Cmfp1,Dmfp1)
G1=ss2tf(sysMFp1)
print(G1)
pole(G1)

#Sistema controlado em MF variação em C
Amfp2=A-B@K
Bmfp2=B*kr
Cmfp2=c
Dmfp2=D
sysMFp2=ss(Amfp2,Bmfp2,Cmfp2,Dmfp2)
G2=ss2tf(Amfp2,Bmfp2,Cmfp2,Dmfp2)
print(G2)
pole(G2)

#Resposta ao degrau em MF para os sistemas com variações paramétricas
YMFp1, tplotMFp1=step(sysMFp1,t)
YMFp2, tplotMFp2=step(sysMFp2,t)

plt.figure()
plt.plot(t,YMF,t,YMFca,t,YMFp1,t,YMFp2)
plt.legend(['Sem controlador','Com controlador e sem var. paramétrica','Com controlador e com var. paramétrica em A','Com controlador e com var. paramétrica em C'])
plt.title('Resposta ao degrau em malha fechada')
plt.ylabel('$Amplitude$')
plt.xlabel('t(s)')
plt.show()


