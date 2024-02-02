close all;
clear all;
clc;

%Declaração do sistema
R=1;
L=1e-3;
Cap=150e-6;

A=[0 1/Cap;-1/L -R/L];
B=[0 1/L]';
C=[1 0];
D=0;
sys=ss(A,B,C,D)

%Polos do sistema em malha aberta
eig(A)

%Matriz de controlabilidade e cálculo do posto da matriz
M=ctrb(A,B);
rank(M)

%Cálculo do vetor de ganhos K
polos_desejados=[-2582 -2582];
K=acker(A,B,polos_desejados)

%Cálculo dos polos do sistema em malha fechada
eig(A-B*K)

%Sistema em malha fechada
Amf=A-B*K;
Bmf=zeros(2,1);
Cmf=C;
Dmf=D;
sysMF=ss(Amf,Bmf,Cmf,Dmf);

%Condições iniciais
x0=[10; 0];

%Tempo de simulação
t=linspace(0,0.02,1000);

%Resposta às condições iniciais
[YMA tplotMA]=initial(sys,x0,t);
[YMF tplotMF]=initial(sysMF,x0,t);

figure;
subplot(2,1,1);
plot(t,YMA,'linewidth',2);
xlabel('t(s)');
ylabel('V_c(t)');
legend('Malha aberta');
subplot(2,1,2);
plot(t,YMF,'linewidth',2);
xlabel('t(s)');
ylabel('V_c(t)');
legend('Malha fechada');