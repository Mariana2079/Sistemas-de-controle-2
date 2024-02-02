close all;
clear all;
clc;

%Declara��o do sistema
A=[0 1;-2 -3];
B=[0; 1];
C=[1 0];
D=0;
sys=ss(A,B,C,D);

%C�lculo da matriz de transi��o de estados
syms s;
phi_s=inv(s*eye(2)-A)
phi=ilaplace(phi_s)

%Outra maneira de calcular a matriz de transi��o de estados
syms t;
phi2=expm(A*t)

%Condi��es iniciais
x0=[1; -1];
%Tempo de simula��o
t=linspace(0,10,100);

%Resposta � fun��o for�ante
[Y_u t_plot X_u]=step(sys,t);
figure;
plot(t_plot,X_u);
title('Resposta � fun��o for�ante');
xlabel('t(s)');
ylabel('Amplitude');
legend('x_1','x_2');


%Resposta �s condi��es iniciais
[Y_0 t_plot_0 X_0]=initial(sys,x0,t);
figure;
plot(t,X_0)
title('Resposta �s condi��es iniciais');
xlabel('t(s)');
legend('x_1','x_2');

%Resposta completa
X_T=X_u+X_0;
figure;
plot(t_plot,X_T);
title('Resposta completa');
xlabel('t(s)');
legend('x_1','x_2');

%Resposta da solu��o anal�tica
k=1;
for tt=0:0.1:10
    eAt(:,:,k)=[(2*exp(-tt)-exp(-2*tt)) (exp(-tt)-exp(-2*tt));
               (-2*exp(-tt)+2*exp(-2*tt)) (-exp(-tt)+2*exp(-2*tt))];
    
    X_0analitica(:,k)=eAt(:,:,k)*x0;
    tplot(k)=tt;
    k=k+1;
end
figure;
plot(tplot,X_0analitica);
title('Resposta anal�tica �s condi��es iniciais');
xlabel('t(s)');
legend('x_1','x_2');

X_uanalitica=[((1/2)-exp(-tplot)+(1/2)*exp(-2*tplot));
               (exp(-tplot)-exp(-2*tplot))];
               
figure;
plot(tplot,X_uanalitica);
title('Resposta anal�tica � fun��o for�ante');
xlabel('t(s)');
legend('x_1','x_2');

X_Tanalitica=X_uanalitica+X_0analitica;
figure;
plot(tplot,X_Tanalitica);
title('Resposta anal�tica completa');
xlabel('t(s)');
legend('x_1','x_2');

