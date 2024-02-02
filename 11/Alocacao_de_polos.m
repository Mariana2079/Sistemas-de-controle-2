close all;
clear all;
clc;

%Declara��o do sistema
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

%Matriz de controlabilidade e c�lculo do posto da matriz
M=ctrb(A,B);
rank(M)

%C�lculo do vetor de ganhos K
polos_desejados=[-2582 -2582];
K=acker(A,B,polos_desejados)

%C�lculo dos polos do sistema em malha fechada
eig(A-B*K)

%Sistema em malha fechada
Amf=A-B*K;
Bmf=zeros(2,1);
Cmf=C;
Dmf=D;
sysMF=ss(Amf,Bmf,Cmf,Dmf);

%Condi��es iniciais
x0=[10; 0];

%Tempo de simula��o
t=linspace(0,0.02,1000);

%Resposta �s condi��es iniciais
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