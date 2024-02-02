clear all;
close all;
clc;

%Declaração da planta
A=[0 1 0;0 0 1;0 -200 -3000];
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