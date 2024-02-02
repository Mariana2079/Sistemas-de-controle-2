clear all;
close all;
clc;

%Declaração da planta
A=[0 1 0;0 0 1;-6 -11 -6];
B=[0 0 10]';
C=[1 0 0];
D=0;

%Polos da planta
eig(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Projeto do Controlador        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sistema aumentado
Abar=[A zeros(3,1);-C 0];
Bbar=[B;0];
Cbar=[C 0];
Dbar=[0];

%Matriz de controlabilidade do sistema aumentado
M=ctrb(Abar,Bbar)
%Posto da matriz de controlabilidade do sistema aumentado
rank(M)

%Polos de MF desejados
polosMF=[-2+j*2*sqrt(3) -2-j*2*sqrt(3) -10 -50];
%Projeto do vetor de ganhos Kbar do controlador
Kbar=acker(Abar,Bbar,polosMF)
vpa(Kbar)

%Verificação do projeto do controlador
eig(Abar-Bbar*Kbar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Projeto do Observador         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Matriz de observabilidade
N=obsv(A,C)
%Posto da matriz de observabilidade
rank(N)

%Autovalores desejados para o observador
autovaloresOBSV=[-10 -10 -50];
Ke=acker(A',C',autovaloresOBSV)'

%Verificação do projeto do observador
eig(A-Ke*C)
