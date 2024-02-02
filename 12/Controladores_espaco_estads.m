%%
close all;
clear all;
clc;

%% Declaracaoo do sistema
% x' = Ax + Bu
% Y  = Cx + Du
A=[0 1;1 0];
B=[0;1];
C=[1 0];
D=0;
sysMA=ss(A,B,C,D);
sysMA

%% Determinar os polos dominantes de malha fechada desejados
% algarismos significativos
qtd_sig = 3;
% Dado Mp=5% e Ts=2s
Mp = 16.3;
Ts = 2;

% Calcular zeta e omegaN
zeta = -log(Mp/100) / sqrt(pi^2 + log(Mp/100)^2);
wn =  4/(Ts * zeta);
disp(zeta);
disp(wn);

% Polos desejados
zeta =1;
wn = 4;
real =-zeta*wn;
real = round(real, qtd_sig);
img = wn*sqrt(1 - zeta^2) * 1j;
img = round(img, qtd_sig);
s = real + img;
disp(s);


%% Polos
% tem qeu mudar a dimensao da matriz eye de acordo com A
syms s;
I = eye(2);
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
polos_desejados=[-4 -4];
K=acker(A,B,polos_desejados);
disp(K);
phi_s=det((s * eye(3) - (A-B*K)));
disp(phi_s);
 
%% ∅(A) = (A − s1) (A − s2) ⋯ (A − sn)
phi_A = A^3 + 9*A^2 + 36*A + 80*eye(3);
disp(phi_A);

%% Representacao em malha fechada
%com compensador
% x' = Ax + Bu
% Y  = Cx + Du  
% u  = −Kx + kr.r
% substituindo u temos 
% x' = (A − Bk)x + kr.B.r
% Y  = Cx

%Sistema controlado em MF com kr=1
Amfc=(A-B*K);
Bmfc=B;
Cmfc=C;
Dmfc=D;
sysMFc=ss(Amfc,Bmfc,Cmfc,Dmfc);
%polos de malha fechada
disp(eig(Amfc));

%% Cálculo de kr 
GanhoCC=dcgain(sysMFc);
disp(GanhoCC);
kr=1/GanhoCC;
disp(kr);
%% ou obter a funcao transferencia em malha fechada com kr=1, entao aplicar
% s=0 para obter o ganhoCC   kr=1/GanhoCC
[NUM,DEN] = ss2tf(Amfc,Bmfc,Cmfc,Dmfc) 
Gmf =tf(NUM,DEN);
s = 0;
disp(evalfr(Gmf,s));
disp(1/evalfr(Gmf,s));

%% Sistema controlado em MF com kr ajustado
%com compensador
% x' = Ax + Bu
% Y  = Cx + Du  
% u  = −Kx + kr.r
% substituindo u temos 
% x' = (A − Bk)x + kr.B.r
% Y  = Cx
Amfca=A-B*K;
Bmfca=B*kr;
Cmfca=C;
Dmfca=D;
sysMFca=ss(Amfca,Bmfca,Cmfca,Dmfca);

%% Sistema em MF sem compensador
% sem compensador
% x' = Ax + Bu
% Y  = Cx + Du  
% u = r − y = r − Cx
% substituindo u temos 
% x' = (A-BC)x + B(r − Cx) 
% Y  = Cx

Amf=A-B*C;
Bmf=B;
Cmf=C;
Dmf=D;
sysMF=ss(Amf,Bmf,Cmf,Dmf);

%% Resposta ao degrau em malha fechada
figure;
subplot(2,1,1);
step(sysMF,20);
hold;
step(sysMFca,20);
legend('Sistema sem controlador','Sistema com controlador (kr=160)');
subplot(2,1,2);
step(sysMFc,20)
legend('Sistema com controlador (kr=1)');