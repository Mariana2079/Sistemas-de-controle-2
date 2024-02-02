clear all;
close all;
clc;

%Declaração da planta
A=[0 1 0;0 0 1;-6 -11 -6];
B=[0 0 1]';
C=[1 0 0];
D=0;

%Polos da planta
eig(A)

%Matriz de observabilidade
N=obsv(A,C)
%Posto da matriz de observabilidade
rank(N)

%Autovalores desejados para o observador
autovaloresOBSV=[-30 -30 -30];
Ke=acker(A',C',autovaloresOBSV)'

%Verificação do projeto
eig(A-Ke*C)


%Plot dos gráficos gerados no Simulink
figure;
subplot(2,1,1);
plot(Planta.time,Planta.signals.values,'linewidth',2)
xlabel('t(s)');
ylabel('Amplitude');
legend('$$x_1$$', '$$x_2$$', '$$x_3$$', 'Interpreter','latex');
subplot(2,1,2);
plot(Estimados.time,Estimados.signals.values,'linewidth',2)
xlabel('t(s)');
ylabel('Amplitude');
legend('$$\tilde{x}_1$$', '$$\tilde{x}_2$$', '$$\tilde{x}_3$$', 'Interpreter','latex');

figure;
plot(Erros.time,Erros.signals.values,'linewidth',2)
xlabel('t(s)');
ylabel('Amplitude');
legend('$$e_{x1}$$', '$$e_{x2}$$', '$$e_{x3}$$', 'Interpreter','latex');

figure;
subplot(2,1,1);
plot(Planta.time,Planta.signals.values,'linewidth',2)
xlabel('t(s)');
ylabel('Amplitude');
legend('$$x_1$$', '$$x_2$$', '$$x_3$$', 'Interpreter','latex');
subplot(2,1,2);
plot(Erros.time,Erros.signals.values,'linewidth',2)
xlabel('t(s)');
ylabel('Amplitude');
legend('$$e_{x1}$$', '$$e_{x2}$$', '$$e_{x3}$$', 'Interpreter','latex');

figure;
subplot(3,1,1);
plot(Planta.time,Planta.signals.values,'linewidth',2)
xlabel('t(s)');
ylabel('Amplitude');
legend('$$x_1$$', '$$x_2$$', '$$x_3$$', 'Interpreter','latex');
subplot(3,1,2);
plot(Estimados.time,Estimados.signals.values,'linewidth',2)
xlabel('t(s)');
ylabel('Amplitude');
legend('$$\tilde{x}_1$$', '$$\tilde{x}_2$$', '$$\tilde{x}_3$$', 'Interpreter','latex');
subplot(3,1,3);
plot(Erros.time,Erros.signals.values,'linewidth',2)
xlabel('t(s)');
ylabel('Amplitude');
legend('$$e_{x1}$$', '$$e_{x2}$$', '$$e_{x3}$$', 'Interpreter','latex');



