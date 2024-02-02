close all;
clear all;
clc;

%Par�metros do sistema
R=0.5;
L=5e-3;
Cap=1500e-6;

%Matrizes do sistema - modelo 1
A1=[0 1/Cap;-1/L -R/L];
B1=[0; 1/L];
C1=[1 0];
D1=0;

%Declara��o do sistema - modelo 1
sys1=ss(A1,B1,C1,D1);

%Resposta ao degrau
[Y1 T1 X1]=step(sys1,0.12);

%Plota sa�da e estados
subplot(3,1,1);
plot(T1,Y1,'linewidth',2);
legend('Sa�da y');
ylabel('Amplitude');
subplot(3,1,2);
plot(T1,X1(:,1),'linewidth',2);
legend('Estado x_1');
ylabel('Amplitude');
subplot(3,1,3);
plot(T1,X1(:,2),'linewidth',2);
legend('Estado x_2');
xlabel('t(s)');
ylabel('Amplitude');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Matrizes do sistema - modelo 2
A2=[0 1;-1/(L*Cap) -R/L];
B2=[0; 1/(L*Cap)];
C2=[1 0];
D2=0;

%Declara��o do sistema - modelo 2
sys2=ss(A2,B2,C2,D2);

%Resposta ao degrau
[Y2 T2 X2]=step(sys2,0.12);

%Plota sa�da e estados
figure;
subplot(3,1,1);
plot(T2,Y2,'linewidth',2);
legend('Sa�da y');
ylabel('Amplitude');
subplot(3,1,2);
plot(T2,X2(:,1),'linewidth',2);
legend('Estado x_1');
ylabel('Amplitude');
subplot(3,1,3);
plot(T2,X2(:,2),'linewidth',2);
legend('Estado x_2');
xlabel('t(s)');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Representa��o por fun��o de transfer�ncia
num=1/(L*Cap);
den=[1 R/L 1/(L*Cap)];
G=tf(num,den);
[y t]=step(G,0.12);
figure;
plot(t,y,'linewidth',2);
axis([0 0.12 0 2]);
legend('Sa�da y');
xlabel('t(s)');
ylabel('Amplitude');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Autovalores modelo 1
eig(A1)
%Autovalores modelo 2
eig(A2)
%Polos da fun��o de transfer�ncia
pole(G)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convers�o de espa�o de estados para fun��o de transfer�ncia
[N D]=ss2tf(A1,B1,C1,D1);
Gconvertida=tf(N,D)

%Convers�o de fun��o de transfer�ncia para espa�o de estados
[Ac Bc Cc Dc]=tf2ss(num,den);
sysc=ss(Ac, Bc, Cc, Dc)