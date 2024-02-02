close all;
clear all;
clc;


%% CANONICA observável
%Descrição do sistema na forma observável
A = [0 0 -20; 1 0 -32; 0 1 -13];
B = [20; 0; 0];
C = [0 0 1];
D = 0;

%Representação na forma observável
sysO=ss(A,B,C,D)

%[aVet, aVal] = eig(A)
% Diagonal do EE

%Matriz de transformação para diagonalização
P = [[20; 12; 1] [10; 11; 1] [2; 3; 1]];

%Obtenção de outra matriz P através do Matlab
%[aVet, aVal] = eig(A)
%P=aVet;

%Sistema diagonalizado
Ad1=inv(P)*A*P;
Bd1=inv(P)*B;
Cd1=C*P;
Dd1=D;
sysD1=ss(Ad1,Bd1,Cd1,Dd1)