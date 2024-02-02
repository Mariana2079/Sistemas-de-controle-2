clear all;
close all;
clc;

%Declaração da planta
A=[0 1 0;0 0 1;0 -2 -3];
B=[0 0 1]';
C=[1 0 0];
D=0;
sysMA=ss(A,B,C,D);

%Matriz de controlabilidade
M=ctrb(A,B)
%Posto da matriz de controlabilidade
rank(M)

%Polos de MF desejados
polosMF=[-2+j*2*sqrt(3) -2-j*2*sqrt(3) -10];
%Projeto do vetor de ganhos K do controlador
K=acker(A,B,polosMF)

%Sistema controlado em MF com kr=1
Amfc=A-B*K;
Bmfc=B;
Cmfc=C;
Dmfc=D;
sysMFc=ss(Amfc,Bmfc,Cmfc,Dmfc);

%Cálculo de kr
GanhoCC=dcgain(sysMFc);
kr=1/GanhoCC

%Sistema controlado em MF com kr ajustado
Amfca=A-B*K;
Bmfca=B*kr;
Cmfca=C;
Dmfca=D;
sysMFca=ss(Amfca,Bmfca,Cmfca,Dmfca);

%Sistema em MF sem controlador
Amf=A-B*C;
Bmf=B;
Cmf=C;
Dmf=D;
sysMF=ss(Amf,Bmf,Cmf,Dmf);

%Resposta ao degrau em MF
figure;
subplot(2,1,1);
step(sysMF,20);
hold;
step(sysMFca,20);
legend('Sistema sem controlador','Sistema com controlador (kr=160)');
subplot(2,1,2);
step(sysMFc,20)
legend('Sistema com controlador (kr=1)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Análise das variações paramétricas
a=[0 1 0;0 0 1;0 -20 -3]; %Variação em A
c=[0.5 0 0]; %Variação em C

%Sistema controlado em MF com variação em A
Amfp1=a-B*K;
Bmfp1=B*kr;
Cmfp1=C;
Dmfp1=D;
sysMFp1=ss(Amfp1,Bmfp1,Cmfp1,Dmfp1);
[n1 d1]=ss2tf(Amfp1,Bmfp1,Cmfp1,Dmfp1);
G1=tf(n1,d1)
pole(G1)

%Sistema controlado em MF com variação em C
Amfp2=A-B*K;
Bmfp2=B*kr;
Cmfp2=c;
Dmfp2=D;
sysMFp2=ss(Amfp2,Bmfp2,Cmfp2,Dmfp2);
[n2 d2]=ss2tf(Amfp2,Bmfp2,Cmfp2,Dmfp2);
G2=tf(n2,d2)
pole(G2)

%Resposta ao degrau em MF para os sistemas com variações paramétricas
figure;
step(sysMF,20);
hold;
step(sysMFca,20);
step(sysMFp1,20)
step(sysMFp2,20)
legend('Sem controlador','Com controlador e sem var. paramétrica','Com controlador e com var. paramétrica em A','Com controlador e com var. paramétrica em C');
title('Resposta ao degrau em malha fechada');