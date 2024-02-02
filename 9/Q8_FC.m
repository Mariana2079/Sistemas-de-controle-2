close all;
clear all;
clc;


%% CANONICA observ�vel
%Descri��o do sistema na forma observ�vel
A = [0 0 -20; 1 0 -32; 0 1 -13];
B = [20; 0; 0];
C = [0 0 1];
D = 0;

%Representa��o na forma observ�vel
sysO=ss(A,B,C,D)

%[aVet, aVal] = eig(A)
% Diagonal do EE

%Matriz de transforma��o para diagonaliza��o
P = [[20; 12; 1] [10; 11; 1] [2; 3; 1]];

%Obten��o de outra matriz P atrav�s do Matlab
%[aVet, aVal] = eig(A)
%P=aVet;

%Sistema diagonalizado
Ad1=inv(P)*A*P;
Bd1=inv(P)*B;
Cd1=C*P;
Dd1=D;
sysD1=ss(Ad1,Bd1,Cd1,Dd1)