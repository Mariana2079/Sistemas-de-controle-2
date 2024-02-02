close all;
clear all;
clc;

%% INICIALIZACAO

num = 2;
den = [1 3 2];
G = tf(num, den);

%% CANONICA CONTROLAVEL
%Descri��o do sistema na forma control�vel
A = [0 1; -2 -3];
B = [0; 1];
C = [2 0];
D = 0;

%Representa��o na forma control�vel
sysC = ss(A,B,C,D)

%Representa��o na forma observ�vel
sysO=ss(A',C',B',D)


% Diagonal do EE

%Matriz de transforma��o para diagonaliza��o
P = [1 1; -1 -2];

%Obten��o de outra matriz P atrav�s do Matlab
%[aVet, aVal] = eig(A)
%P=aVet;

%Sistema diagonalizado
Ad1=inv(P)*A*P;
Bd1=inv(P)*B;
Cd1=C*P;
Dd1=D;
sysD1=ss(Ad1,Bd1,Cd1,Dd1)