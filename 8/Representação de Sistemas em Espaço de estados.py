# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 20:42:40 2020

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

#Parâmetros do sistema
R=0.5
L=5e-3
Cap=1500e-6

#Matrizes do sistema - modelo 1
A1=[[0, 1/Cap], [-1/L, -R/L]]
B1=[[0], [1/L]]
C1=[1, 0]
D1=0

#Declaração do sistema - modelo 1
sys1=ss(A1,B1,C1,D1)

#Resposta ao degrau
t=np.linspace(0,0.12,1000)
y1, t1, x1 = step(sys1,t,return_x=True)
plt.figure()
plt.subplot(3,1,1)
plt.plot(t1,y1)
plt.legend(['Saída y'])
plt.ylabel('Amplitude')
plt.subplot(3,1,2)
plt.plot(t1,x1[:,0])
plt.legend(['Estado $x_1$'])
plt.ylabel('Amplitude')
plt.subplot(3,1,3)
plt.plot(t1,x1[:,1])
plt.legend(['Estado $x_2$'])
plt.xlabel('t(s)')
plt.ylabel('Amplitude')
plt.show()

###########################################################################


#Matrizes do sistema - modelo 2
A2=[[0, 1], [-1/(L*Cap), -R/L]]
B2=[[0], [1/(L*Cap)]]
C2=[1, 0]
D2=0

#Declaração do sistema - modelo 2
sys2=ss(A2,B2,C2,D2)

#Resposta ao degrau
y2, t2, x2 = step(sys2,t,return_x=True)
plt.figure()
plt.subplot(3,1,1)
plt.plot(t2,y2)
plt.legend(['Saída y'])
plt.ylabel('Amplitude')
plt.subplot(3,1,2)
plt.plot(t2,x2[:,0])
plt.legend(['Estado $x_1$'])
plt.ylabel('Amplitude')
plt.subplot(3,1,3)
plt.plot(t2,x2[:,1])
plt.legend(['Estado $x_2$'])
plt.xlabel('t(s)')
plt.ylabel('Amplitude')
plt.show()


#############################################################################

#Declaração da função de transferência
num=1/(L*Cap)
den=[1, R/L, 1/(L*Cap)]
G=tf(num,den)
y3, t3 = step(G,t)
plt.figure()
plt.plot(t3,y3)
plt.legend(['Saída y'])
plt.xlabel('t(s)')
plt.ylabel('Amplitude')
plt.axis([0, 0.12, 0, 2])
plt.show()


#############################################################################
#Autovalores modelo 1
print(la.eigvals(A1))

#Autovalores modelo 2
print(la.eigvals(A2))

#Polos da função de transferência
print(pole(G))

#############################################################################
#Conversão de espaço de estados para função de transferência
Gconvertida=ss2tf(A1,B1,C1,D1)
print(Gconvertida)


#Conversão de função de transferência para espaço de estados
sysc=tf2ss(G)
print(sysc)