# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 16:54:58 2020

@author: rcard
"""

#As próximas 3 linhas são para selecionar entre plot inline ou em nova janela
#Útil para rlocus
from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'qt')
#get_ipython().run_line_magic('matplotlib', 'inline')

import numpy as np                #Biblioteca para cálculo numérico (matrizes)
import math                       #Funções matemáticas
import matplotlib.pyplot as plt   #Funções de plot similares as do MATLAB
import control as ctrl            #Biblioteca para controle
from control.matlab import *      #Funções para controle similares as do MATLAB
import scipy.linalg as la         #Biblioteca de álgebra 
import scipy.signal as sps        #Biblioteca de processamento de sinais
import sympy as sym               #Biblioteca para cálculo simbólico


#Cálculo da matriz de transição de estados
a=sym.Matrix([[0, 1],[-2, -3]])
s, t= sym.symbols('s t')
phi_s=(s*sym.eye(2)-a).inv()
phi_s=sym.cancel(phi_s)
phi=sym.inverse_laplace_transform(phi_s,s,t)
print(phi)

#Outra forma de calcular a matriz de transição de estados
phi2=sym.exp(a*t)
print(phi2)

#Declaração do sistema
A=[[0, 1],[-2, -3]]
B=[[0], [1]]
C=[[1, 0]]
D=[[0]]
sys=ss(A,B,C,D)

#Condições iniciais
x0=[[1], [-1]]

#Rempo de simulação
t=np.linspace(0,10,100)

#Resposta da função forçante
Yu, tplot, Xu=step(sys,t,return_x=True)
plt.figure()
plt.plot(tplot,Xu)
plt.title('Resposta à função forçante')
plt.xlabel('t(s)')
plt.ylabel('Amplitude')
plt.legend(['$x_1$','$x_2$'])
plt.show()


#Resposta da condição inicial
Y0, tplot0, X0=initial(sys,t,x0,return_x=True)
plt.figure()
plt.plot(t,X0)
plt.title('Resposta às condições iniciais')
plt.xlabel('t(s)')
plt.legend(['$x_1$','$x_2$'])
plt.show()

#Resposta completa
XT=Xu+X0
plt.figure()
plt.plot(tplot,XT)
plt.title('Resposta completa')
plt.xlabel('t(s)')
plt.legend(['$x_1$','$x_2$'])
plt.show()
