clear all;
close all;
clc;

%% Declaração do sistema
A = [0 1; -3 -4];
B = [0 1]';
C = [3 0];
D = 0;
sysMA = ss(A, B, C, D);
sysMA

n = 2;


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

%% Polos,  autovalores da matriz
eig(A)
eig(Amf)

%% ############## PROJETO DO CONTROLADOR ##############
%% Ampliacao do sistema para incluir integrador
mzero = zeros(n,1);
Abarra = [A mzero; -C 0];
Bbarra = [B; 0];
Cbarra = [C 0];
Dbarra = 0;

% sistema aumentado
sysMAaum = ss(Abarra, Bbarra, Cbarra, Dbarra);
sysMAaum

% quantidade de estados da planta aumentada
n_aum = 3
eig(sysMAaum)


%% Polos da Planta original
% tem qeu mudar a dimensao da matriz eye de acordo com A
syms s;
I = eye(n);
eqn = det((s * I - A));
polos = vpasolve(eqn, s);
disp(polos);

% outra maneira de fazer para achar os auto valores de A
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

M = [Bbarra Abarra*Bbarra (Abarra^2)*Bbarra];
det(M)       

M = ctrb(Abarra,Bbarra);
rank(M)             


%% Aplicar a Formula de Ackermann K = [0 0 0 ... 1][B | AB | ... | A^(n-1)B] phi(A)
% M = [B | AB | ... | A^(n-1)B]^1
% onde phi(s) o polinomio caracteristico formado a partir dos polos desejados
% de malha fechada (s1 , s2 , ... , sn )
% onde  phi(s) s foi substituido por A, ficando phi(A)
s1 = -2 + 2*1i;
s2 = -2 - 2*1i;
s3 = -10;
polos_desejados = [s1 s2 s3];
K = acker(Abarra, Bbarra, polos_desejados);
disp(vpa(K));
%K = place(A, B, polos_desejados);

% Obter phi(s) e phi(A)
syms s
I = eye(n_aum);
phi_s=vpa(det((s * I - (Abarra-Bbarra*K))));
disp(phi_s);

% n = qtd variaves de estados

I = eye(n_aum);
s = Abarra;
phi_A  = s^3 + 14.0*s^2 + 48.0*s + 80.0*I;
disp(phi_A);



%% Representacao em malha fechada do sistema aumentado controlado em MF

%Sistema aumentado controlado em MF com integrador
Amfca = Abarra-Bbarra*K;
Bmfca = [mzero; 1]; % Rbarra = [zero(n, 1), 1]
Cmfca = Cbarra;
Dmfca = Dbarra;
sysMFca = ss(Amfca, Bmfca, Cmfca, Dmfca);

GanhoCC = dcgain(sysMFca);

fprintf('O ganho do sistema aumentado em malha fechada\n')
disp(double(GanhoCC));

fprintf('Polos sistema aumentado em malha fechada\n')
eig(sysMFca)


%% ############## ANALISE DAS VARIACOES PARAMETRICAS ##############
%a=[0 1 0;0 0 1;0 -20 -3]; %Variacao em A
b = [0; 1.5]; %Variacao em B
% c = [0.5 0 0]; %Variação em C

%Sistema controlado em MF com variação em A
% Amfp1=a-B*K;
% Bmfp1=B*kr;
% Cmfp1=C;
% Dmfp1=D;
% sysMFp1=ss(Amfp1,Bmfp1,Cmfp1,Dmfp1);
% [n1 d1]=ss2tf(Amfp1,Bmfp1,Cmfp1,Dmfp1);
% G1=tf(n1,d1)
% pole(G1)

% %Sistema controlado em MF com variação em B
mzero = zeros(n,1);

% Valores atualizado do sistema aumentado 
Abarra2 = [A mzero; -C 0];
Bbarra2 = [b; 0];
Cbarra2 = [C 0];
Dbarra2 = 0;


% %Sistema controlado em MF com variação em C
% mzero = zeros(n,1);
% 
% % Valores atualizado do sistema aumentado 
% Abarra2 = [A mzero; -c 0];
% Bbarra2 = [B; 0];
% Cbarra2 = [c 0];
% Dbarra2 = 0;

% malha fechada com variacao parametrica
Amfp2 = Abarra2-Bbarra2*K;
Bmfp2 = [mzero; 1]; % Rbarra = [zero(n, 1), 1];
Cmfp2 = Cbarra2;
Dmfp2 = Dbarra2;
sysMFp2 = ss(Amfp2,Bmfp2,Cmfp2,Dmfp2);
[n2 d2] = ss2tf(Amfp2,Bmfp2,Cmfp2,Dmfp2);
G2 = tf(n2, d2);
pole(G2)


%% ############## PROJETO DO OBSERVADOR ##############
%%  Verificar observabilidade do sistema, com base no projeto do controlador onde for LETRA = LETRAmfca
% encontrar n, onde n e o numero de estados 
% -> n = dim(x)
% montar a matriz de observabilidade N
% N = [C | CA | ... | CA^(n-1)] 

% Se
% posto(N) = n -> Sistema observavel
% det(N) != 0  -> Sistema observavel
% Nao a necessidade de verificar os dois, apenas 1 ja e valido

n2 = n_aum;

%Matriz de observabilidade
% N = [C; C*A];
N = obsv(Amfca,Cmfca)
rank(N)     
det(N)


%% Definir autovalores do observador
eig(Amfca)

% autovalor 5 vezes autovalores do sistema em malha fechado compensado (5 * autovalorMF)
u1 = 5 * (-2);
u2 = 5 * (-2);
u3 = 5 * (-10);


%% Calcular polinomio caracteristico desejado
syms s
phi_s = conv([1 -u1], [1 -u2])
phi_s = conv(phi_s, u3)
poly2sym(phi_s, s)

% lembrar de matriz identidade multiplicar constante
I = eye(n2)
s = A;
phi_a = - 50*s^2 - 1000*s - 5000*I

%% Aplicar a Formula de Ackermann Ke = phi(A)N [0; 0; ...0; 1]
% N = [C | CA | ... | CA^(n-1)] ^-1
% onde phi(s) o polinomio caracteristico formado a partir dos polos desejados
% de malha fechada (s1 , s2 , ... , sn )
% onde  phi(s) s foi substituido por A, ficando phi(A)


Ke = phi_a * (inv(N) * [0; 0; 1]);
%Ke = phi_a * (N\[0; 0; 1]);
vpa(Ke)

% todos os passo anteriores se resumem no acker
autovaloresOBSV = [u1 u2 u3];
Ke = acker(A', C', autovaloresOBSV)';
disp(Ke);

%% Representacao do observador
% xtil' = (A - Ke C)xtil + Bu + Ke y
% xtil' = (A - Ke C)xtil + [B Ke] [u; y]
% ytil = I xtil

I = eye(n2)

Aob = (A - Ke*C);
Bob = [B Ke];
Cob = I;
Dob = 0;
sysOB = ss(Aob, Bob, Cob, Dob);
sysOB

% Vale ressalta que é utililizado uma tecnica de calculo numerico para
% resolver a equacao de estados, para isso e necessario discretiza-la.

%% Verificacao
%Autovalores de Aob = (A - Ke*C);
eig(Aob)
vpa(eig(Aob))
% como verificado os autovalores deram iguais a -50, validando o resultado
% esperado
