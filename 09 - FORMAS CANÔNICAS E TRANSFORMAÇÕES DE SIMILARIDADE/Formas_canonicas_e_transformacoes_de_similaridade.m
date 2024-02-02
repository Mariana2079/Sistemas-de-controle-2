close all;
clear all;
clc;

%Sistema descrito por fun��o de transfer�ncia
num=[1 1];
den=[1 6 9];
G=tf(num,den)

%Descri��o do sistema na forma control�vel
A=[0 1 0;0 0 1;-60 -47 -12];
B=[0; 0; 1];
C=[2 3 1];
D=0;

%Representa��o na forma control�vel
sysC=ss(A,B,C,D)

%Representa��o na forma observ�vel
sysO=ss(A',C',B',D)

%Matriz de transforma��o para diagonaliza��o
P=[1 1 1;-3 -4 -5;9 16 25];

%Obten��o de outra matriz P atrav�s do Matlab
[aVet, aVal]=eig(A)
P=aVet;

%Sistema diagonalizado
Ad1=inv(P)*A*P;
Bd1=inv(P)*B;
Cd1=C*P;
Dd1=D;
sysD1=ss(Ad1,Bd1,Cd1,Dd1)

%Representa��o na forma diagonal a partir
%da fun��o de transfer�ncia
[Res,Polos,K] = residue(num,den) 

Ad2=[-3 0 0;0 -4 0;0 0 -5];
Bd2=[1; 1; 1];
Cd2=[1 -6 6];
Dd2=0;
sysD2=ss(Ad2,Bd2,Cd2,Dd2)

%Compara��o das respostas ao degrau para
%diferentes representa��es de um mesmo sistema
[YC tC]=step(sysC,2);
[YO tO]=step(sysO,2);
[YD1 tD1]=step(sysD1,2);
[YD2 tD2]=step(sysD2,2);
[YG tG]=step(G,2);
subplot(5,1,1);
plot(tC,YC,'linewidth',2);
legend('Sa�da y - forma control�vel');
ylabel('Amplitude');
subplot(5,1,2);
plot(tO,YO,'linewidth',2);
legend('Sa�da y - forma observ�vel');
ylabel('Amplitude');
subplot(5,1,3);
plot(tD1,YD1,'linewidth',2);
legend('Sa�da y - forma diagonal 1');
ylabel('Amplitude');
subplot(5,1,4);
plot(tD2,YD2,'linewidth',2);
legend('Sa�da y - forma diagonal 2');
ylabel('Amplitude');
subplot(5,1,5);
plot(tG,YG,'linewidth',2);
legend('Sa�da y - G(s)');
xlabel('t(s)');
ylabel('Amplitude');

%Representa��o na forma de Jordan a partir
%da fun��o de transfer�ncia
numj=[1 3 2];
denj=[1 11 39 45];
Gj=tf(numj,denj)
[Rj,Pj,Kj] = residue(numj,denj) 

%Sistema na forma de Jordan
Aj=[-3 1 0;0 -3 0;0 0 -5];
Bj=[0;1;1];
Cj=[1 -2 3];
Dj=0;
sysj=ss(Aj,Bj,Cj,Dj);

%Compara��o das respostas ao degrau do sistema representado
%por fun��o de transfer�ncia e na forma de Jordan
[YGj tGj]=step(Gj,2);
[Yj tj]=step(sysj,2);
figure;
subplot(2,1,1);
plot(tGj,YGj,'linewidth',2);
legend('Sa�da y - G(s)');
ylabel('Amplitude');
subplot(2,1,2);
plot(tj,Yj,'linewidth',2);
legend('Sa�da y - forma de Jordan');
xlabel('t(s)');
ylabel('Amplitude');