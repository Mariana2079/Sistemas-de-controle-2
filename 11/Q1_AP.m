close all;
clear all;
clc;

%% Declaração do sistema
A = [0 1; 4 0];
B = [0; 1];
C = [2 0];
D = 0;
sysMA = ss(A, B, C, D);
sysMA

%% Polos
% tem qeu mudar a dimensao da matriz eye de acordo com A
syms s;
I = eye(2);
eqn = det((s * I - A));
polos = vpasolve(eqn, s);
disp(polos);

% outra maneira de fazer Ã© achar os auto valores de A
disp(eig(A));


%% Verificar contralabilidade do sistema
% encontrar n, onde n e o numero de estados 
% -> n = dim(x)
% montar a matriz de controlabilidade
% M = [B | AB | ... | A^(n-1)B] 

% Se
% posto(M) = n -> Sistema controlavel
% det(M) != 0  -> Sistema controlavel
% Nao a necessidade de verificar os dois, apenas 1 ja e valido

n = 2;

M = [B A*B];
det(M)      

M = ctrb(A,B);
rank(M)            


%% Aplicar a Formula de Ackermann K = [0 0 0 ... 1][B | AB | ... | A^(n-1)B] phi(A)
% M = [B | AB | ... | A^(n-1)B]^1
% onde phi(s) o polinomio caracteristico formado a partir dos polos desejados
% de malha fechada (s1 , s2 , ... , sn )
% onde  phi(s) s foi substituido por A, ficando phi(A)

s1 = -2;
s2 = -2;
% s3 = -20;

polos_desejados = [s1 s2];
K = acker(A, B, polos_desejados);

% Obter phi(s) e phi(A)
phi_s=det((s * I - (A-B*K)));
disp(phi_s);
% copiar saida e phi_a
% n = qtd variaves de estados
syms s
I = eye(n);
s = A;
phi_A  = s^2 + 4*s + 4*I;
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


Amf = (A - B*K);
Bmf = zeros(n, 1);
Cmf = C;
Dmf = D;
sysMF = ss(Amf, Bmf, Cmf, Dmf);

%polos de malha fechada
disp(eig(Amf));

% [NUM,DEN] = ss2tf(A, B, C, D) 
% Gma =tf(NUM,DEN);
% [NUMmf,DENmf] = ss2tf(Amf,Bmf,Cmf,Dmf) 
% Gmf =tf(NUMmf,DENmf);




