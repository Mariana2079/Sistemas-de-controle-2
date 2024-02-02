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

#Polos da planta
print(la.eigvals(A))

#########################################
#         Projeto do Controlador        #
#########################################

#Sistema aumentado
Abar=np.concatenate((A, -C))
Abar=np.concatenate((Abar, np.zeros((4,1))), axis=1)
Bbar=np.concatenate((B, [[0]]))
Cbar=np.concatenate((C, [[0]]), axis=1)
Dbar=D

#Matriz de controlabilidade do sistema aumentado
M=ctrb(Abar,Bbar)
#Posto da matriz de controlabilidade
print(np.linalg.matrix_rank(M))

#Polos de MF desejados
polosMF=[-2+2*math.sqrt(3)*1j, -2-2*math.sqrt(3)*1j, -10, -50]
#Projeto do vetor de ganhos Kbar do controlador
Kbar=acker(Abar, Bbar, polosMF)
print(Kbar)

#Verificação do projeto do controlador
print(la.eigvals(Abar-Bbar@Kbar))

#########################################
#         Projeto do Observador         #
#########################################

#Matriz de observabilidade
N=obsv(A,C)
#Posto da matriz de observabilidade
print(np.linalg.matrix_rank(N))

#Autovalores do observador
autovaloresOBSV=[-10, -10, -50]

#Projeto do vetor de ganhos K do controlador
Ke=acker(A.transpose(), C.transpose(), autovaloresOBSV).transpose()
print(Ke)

#Verificação do projeto
print(la.eigvals(A-Ke@C))