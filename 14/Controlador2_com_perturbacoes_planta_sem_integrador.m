clear all;
close all;
clc;

%Declaração da planta
A=[0 1 0;0 0 1;-6 -11 -6];
B=[0 0 1]';
C=[1 0 0];
D=0;
sysMA=ss(A,B,C,D);
%Polos da planta
eig(A)

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

%Sistema em MF sem controlador
Amf=A-B*C;
Bmf=B;
Cmf=C;
Dmf=D;
sysMF=ss(Amf,Bmf,Cmf,Dmf);
%Polos do sistema em MF sem controlador
eig(Amf)

%Sistema em MF com controlador
Amfc=Abar-Bbar*Kbar;
Bmfc=[0 0 0 1]';
Cmfc=Cbar;
Dmfc=Dbar;
sysMFc=ss(Amfc,Bmfc,Cmfc,Dmfc);
%Polos do sistema em MF com controlador
eig(Amfc)
[n d]=ss2tf(Amfc,Bmfc,Cmfc,Dmfc);
Gmfc=tf(n,d)
pole(Gmfc)
dcgain(Gmfc)



%Resposta ao degrau em MF
figure;
step(sysMF,20);
hold;
step(sysMFc,20);
legend('Sistema sem controlador','Sistema com controlador');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Análise das variações paramétricas
a=[0 1 0;0 0 1;-60 -11 -6]; %Variação em A
abar=[a zeros(3,1);-C 0];

%Sistema controlado em MF com variação em A
Amfp1=abar-Bbar*Kbar;
Bmfp1=[0 0 0 1]';
Cmfp1=Cbar;
Dmfp1=Dbar;
sysMFp1=ss(Amfp1,Bmfp1,Cmfp1,Dmfp1);
[n1 d1]=ss2tf(Amfp1,Bmfp1,Cmfp1,Dmfp1);
G1=tf(n1,d1)
pole(G1)
dcgain(G1)


c=[0.5 0 0]; %Variação em C
abar2=[A zeros(3,1);-c 0];
cbar=[c 0];

%Sistema controlado em MF com variação em C
Amfp2=abar2-Bbar*Kbar;
Bmfp2=[0 0 0 1]';
Cmfp2=cbar;
Dmfp2=Dbar;
sysMFp2=ss(Amfp2,Bmfp2,Cmfp2,Dmfp2);
[n2 d2]=ss2tf(Amfp2,Bmfp2,Cmfp2,Dmfp2);
G2=tf(n2,d2)
pole(G2)
dcgain(G2)

%Resposta ao degrau em MF para os sistemas com variações paramétricas
figure;
step(sysMF,20);
hold;
step(sysMFc,20);
step(sysMFp1,20)
step(sysMFp2,20)
legend('Sem controlador','Com controlador e sem var. paramétrica','Com controlador e com var. paramétrica em A','Com controlador e com var. paramétrica em C');
title('Resposta ao degrau em malha fechada');