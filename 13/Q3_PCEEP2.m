close all;
clear all;
clc;

%% Determinar os polos dominantes de malha fechada desejados
% quantidade de significativos
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
A = [0 1; -8 -4];
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
n = 2;
I = eye(n);
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

polos_desejados = [s1 s2];
K = acker(A, B, polos_desejados);
disp(vpa(K));
%K = place(A, B, polos_desejados);
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
sysMFc
eig(sysMFc)

%% Cálculo de kr e Planta de malha fechada com compensador MFca
GanhoCC = dcgain(sysMFc);
disp(double(GanhoCC));
kr = 1/GanhoCC;
disp(vpa(kr))


%Sistema controlado em MF com kr ajustado
Amfca = A - B*K;
Bmfca = B*kr;
Cmfca = C;
Dmfca = D;
sysMFca = ss(Amfca, Bmfca, Cmfca, Dmfca);
sysMFca
GanhoCC = dcgain(sysMFca);
disp(GanhoCC);



%% Erro em regime permanente

% Teorema do valor final
% E = R - Y = R - R*Gmf
% E = R(1 - Gmf)
% e(inf) = s*E = s*R*(1 - Gmf)
syms s
R = 1/s;
I = eye(n);
E = R * (1 - (Cmfca*(inv(s*I - Amfca))*Bmfca));
erro = limit(s * E, s, 0);
disp(erro)


%% Saida em regime permanente
% Teorema do valor final
% Y = e(inf) = lim s * R * E, s->0

syms s
I = eye(n);
E = vpa((Cmfca*(inv(s*I - Amfca))*Bmfca));

% R = Entrada;
% degrau = 1/s; rampa = 1/s^2 ; parabola = 1/s^3
R = 1/s;
Y = limit(s * R * E , s, 0);
disp(double(Y));

%% Resposta ao degrau em MF
figure;
subplot(2,1,1);
step(sysMF,20);
hold;
step(sysMFca,20);
legend('Sistema sem controlador','Sistema com controlador (kr=2600)');
subplot(2,1,2);
step(sysMFc,20)
legend('Sistema com controlador (kr=1)');


%% Análise das variações parametricas
%a=[0 1 0;0 0 1;0 -20 -3]; %Variacao em A
b = [0; 1.5]; %Variacao em B
% c = [0.5 0 0]; %Variacao em C

%Sistema controlado em MF com variação em A
% Amfp1=a-B*K;
% Bmfp1=B*kr;
% Cmfp1=C;
% Dmfp1=D;
% sysMFp1=ss(Amfp1,Bmfp1,Cmfp1,Dmfp1);
% [n1 d1]=ss2tf(Amfp1,Bmfp1,Cmfp1,Dmfp1);
% G1=tf(n1,d1)
% pole(G1)

%Sistema controlado em MF com variação em C
Amfp2 = A-b*K;
Bmfp2 = b*kr;
Cmfp2 = C;
Dmfp2 = D;
sysMFp2 = ss(Amfp2,Bmfp2,Cmfp2,Dmfp2);
[n2 d2] = ss2tf(Amfp2,Bmfp2,Cmfp2,Dmfp2);
G2 = tf(n2,d2)
pole(G2)


%% Erro em regime permanente das variações paramétricas

% Teorema do valor final
% E = R - Y = R - R*Gmf
% E = R(1 - Gmf)
% e(inf) = s*E = s*R*(1 - Gmf)
syms s
R = 1/s;
I = eye(n);
E = R * (1 - (Cmfp2*(inv(s*I - Amfp2))*Bmfp2));
erro = limit(s * E, s, 0);
disp(double(erro))

%% Saida em regime permanente das variações paramétricas

syms s
I = eye(n);
E = vpa((Cmfp2*(inv(s*I - Amfp2))*Bmfp2));

% R = Entrada;
% degrau = 1/s; rampa = 1/s^2 ; parabola = 1/s^3
R = 1/s;
erro = limit(s * R * E , s, 0);
disp(vpa(erro));

