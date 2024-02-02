close all;
clear all;
clc;

%Declaração do sistema
A=[0 1;-2 -3];
B=[0; 1];
C=[1 0];
D=0;
sys=ss(A,B,C,D);

%Cálculo da matriz de transição de estados
syms s;
phi_s=inv(s*eye(2)-A)
phi=ilaplace(phi_s)

%Outra maneira de calcular a matriz de transição de estados
syms t;
phi2=expm(A*t)

%Condições iniciais
x0=[1; -1];
%Tempo de simulação
t=linspace(0,10,100);

%Resposta à função forçante
[Y_u t_plot X_u]=step(sys,t);
figure;
plot(t_plot,X_u);
title('Resposta à função forçante');
xlabel('t(s)');
ylabel('Amplitude');
legend('x_1','x_2');


%Resposta às condições iniciais
[Y_0 t_plot_0 X_0]=initial(sys,x0,t);
figure;
plot(t,X_0)
title('Resposta às condições iniciais');
xlabel('t(s)');
legend('x_1','x_2');

%Resposta completa
X_T=X_u+X_0;
figure;
plot(t_plot,X_T);
title('Resposta completa');
xlabel('t(s)');
legend('x_1','x_2');

%Resposta da solução analítica
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
title('Resposta analítica às condições iniciais');
xlabel('t(s)');
legend('x_1','x_2');

X_uanalitica=[((1/2)-exp(-tplot)+(1/2)*exp(-2*tplot));
               (exp(-tplot)-exp(-2*tplot))];
               
figure;
plot(tplot,X_uanalitica);
title('Resposta analítica à função forçante');
xlabel('t(s)');
legend('x_1','x_2');

X_Tanalitica=X_uanalitica+X_0analitica;
figure;
plot(tplot,X_Tanalitica);
title('Resposta analítica completa');
xlabel('t(s)');
legend('x_1','x_2');

