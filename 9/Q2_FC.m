close all;
clear all;
clc;

%% INICIALIZACAO

num = 2;
den = [1 3 2];
G = tf(num, den);

%% CANONICA CONTROLAVEL
%Descrição do sistema na forma controlável
A = [0 1; -2 -3];
B = [0; 1];
C = [2 0];
D = 0;

%Representação na forma controlável
sysC = ss(A,B,C,D)

%Representação na forma observável
sysO=ss(A',C',B',D)


% Diagonal do EE

%Matriz de transformação para diagonalização
P = [1 1; -1 -2];

%Obtenção de outra matriz P através do Matlab
%[aVet, aVal] = eig(A)
%P=aVet;

%Sistema diagonalizado
Ad1=inv(P)*A*P;
Bd1=inv(P)*B;
Cd1=C*P;
Dd1=D;
sysD1=ss(Ad1,Bd1,Cd1,Dd1)