close all;
clear all;
clc;

%% Determinar os polos dominantes de malha fechada desejados
% quantidade de significativos
qtd_sig = 3;

% Dado Mp=n% e Ts=Ns
Mp = 0;
% Ts = 2;

% Calcular zeta e omegaN
zeta = -log(Mp/100) / sqrt(pi^2 + log(Mp/100)^2);
wn =  4/(Ts * zeta);
disp(zeta);
disp(wn);

zeta = 1;
wn =  4;
% Polos desejados

real = -zeta*wn;
real = round(real, qtd_sig);
img = wn*sqrt(1 - zeta^2) * 1j;
img = round(img, qtd_sig);
s = real + img;
disp(s);

% arredondamento requerido 3 algarismos significativo
s_real = -2;
s_img =  3.464;
s = s_real + s_img*1j;


%% Declaração do sistema Malha Aberta
A = [0 1; 1 0];
B = [0; 1];
C = [1 0];
D = 0;
sysMA = ss(A, B, C, D);
sysMA

%% Representacao em malha fechada s/ controlador
% x' = Ax + Bu
% Y  = Cx + Du
% u = r - y = r - Cx
% substituindo u temos 
% x' = Ax + B(r - Cx) 
% malha fechada s/ compesador
% x' = (A-(B*C))x + B*r
% Y  = Cx

Amf = (A - B*C);
Bmf = B;
Cmf = C;
Dmf = D;
sysMF = ss(Amf, Bmf, Cmf, Dmf);

eig(sysMF)

%% Polos
% tem qeu mudar a dimensao da matriz eye de acordo com A
syms s;
I = eye(3);
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

s1 = -4;
% s2 = -2 + 3.464*1j;
% s3 = -50;

polos_desejados = [s1 s1];
K = acker(A, B, polos_desejados);
disp(vpa(K));
K = place(A, B, polos_desejados);

% Obter phi(s) e phi(A)
syms s
phi_s=vpa(det((s * I - (A-B*K))));
disp(phi_s);

% n = qtd variaves de estados

I = eye(n);
s = A;
phi_A  = s^3 + 9*s^2 + 36*s + 80*I;
disp(phi_A);



%% Representacao em malha fechada com Compesador e achar ganho CC e Kr
% x' = Ax + Bu
% Y  = Cx + Du
% u=-Kx +  Kr*r
% substituindo u temos 
% x' = Ax + B(-Kx + Kr*r) 
% malha fechada
% x' = (A-(B*K))x + B*Kr*r
% Y  = Cx

%Sistema controlado em MF com kr=1
Amfc = A-B*K;
Bmfc = B;
Cmfc = C;
Dmfc = D;
sysMFc = ss(Amfc, Bmfc, Cmfc, Dmfc);
eig(sysMFc)

%% Cálculo de kr
GanhoCC = dcgain(sysMFc);
disp(double(GanhoCC));
kr = 1/GanhoCC;
disp(vpa(kr))

syms s
[numMF, denMF] = ss2tf(Amf, Bmf, Cmf, Dmf);
G = tf(numMF, denMF);
s = 0;
Gcc = evalfr(G, s);
kr2 = 1/Gcc


%Sistema controlado em MF com kr ajustado
Amfca = A - B*K;
Bmfca = B*kr;
Cmfca = C;
Dmfca = D;
sysMFca = ss(Amfca, Bmfca, Cmfca, Dmfca);


%% Resposta ao degrau em MF
figure;
subplot(2,1,1);
step(sysMF,20);
hold;
step(sysMFca,20);
legend('Sistema sem controlador','Sistema com controlador (kr=16)');
subplot(2,1,2);
step(sysMFc,20)
legend('Sistema com controlador (kr=1)');



