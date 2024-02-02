%%
close all;
clear all;
clc;

%% Declaracaoo do sistema
% x' = Ax + Bu
% Y  = Cx + Du
A=[0 1 0;0 0 1;-8 -8 -4];
B=[0; 0; 1];
C=[4 0 0];
D=0;
sysMA=ss(A,B,C,D);
sysMA

%% Polos
% tem qeu mudar a dimensao da matriz eye de acordo com A
syms s;
I = eye(3);
eqn = det((s * I - A));
polos = vpasolve(eqn, s);
disp(polos);

% outra maneira de fazer é achar os auto valores de A
disp(eig(A));

%% Verificar controlabilidade do sistema
% encontrar n, onde n é o numero de estados 
% -> n = dim(x)
% montar a matriz de controlabilidade
% M = [ B : A.B : ... : A^(n-1).B]

% Se
% posto(M) = n -> Sistema controlável
% det(M) ≠ 0  -> Sistema controlável

%matriz de controlabilidade
M=ctrb(A,B);
disp(M);
% Calcular posto
disp(rank(M));

%% Fórmula de Ackermann
% M = [0 0 .. 1][ B : A.B : ... : A^(n-1).B]⁻¹ ∅(s)
% onde ∅(s) é o polinômio característico formado
% a partir dos polos desejados
% de malha fechada (s1 , s2 , ⋯ , sn ),
% onde  ∅(s) s foi substituido por A,ficando ∅(A)
% M = [0 0 .. 1][ B : A.B : ... : A^(n-1).B]⁻¹ ∅(A)

% ∅(s) = (s − s1)(s − s2) ⋯ (s − sn)
polos_desejados=[-2 -2 -20];
K=acker(A,B,polos_desejados)

phi_s=det((s * I - (A-B*K)));
disp(phi_s);
 
% ∅(A) = (A − s1) (A − s2) ⋯ (A − sn)
phi_A = A^3 + 24*A^2 + 84*A + 80 *eye(3);
disp(phi_A);

%% Representacao em malha fechada
% x' = Ax + Bu
% Y  = Cx + Du
% u=-Kx
% substituindo u temos 
% x' = Ax + B(-Kx) 
% malha fechada
% x' = (A-(B*K))x
% Y  = Cx
Amf=(A-B*K);
Bmf=zeros(3,1);
Cmf=C;
Dmf=D;
sysMF=ss(Amf,Bmf,Cmf,Dmf);

%polos de malha fechada
disp(eig(Amf));

[NUM,DEN] = ss2tf(A,B,C,D) 
Gma =tf(NUM,DEN);
[NUMmf,DENmf] = ss2tf(Amf,Bmf,Cmf,Dmf) 
Gmf =tf(NUMmf,DENmf);

figure
step(Gma);
hold;
step(Gmf);
legend('Malha aberta','Malha fechada');

%% 
%Condições iniciais
x0=[10; 0 ;0];

%Tempo de simulação
t=linspace(0,0.02,10000);

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