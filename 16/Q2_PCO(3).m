clear all;
close all;
clc;

%% Declaração do sistema
A = [0 1 0; 0 0 1; -12 -19 -8];
B = [0 0 1]';
C = [2 1 0];
D = 0;
sysMA = ss(A, B, C, D);
sysMA

n = 3;


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
sysMF = ss(Amf, Bmf, Cmf, Dmf)




%% Polos,  autovalores da matriz
eig(A)
eig(Amf)

%% ############## PROJETO DO CONTROLADOR ##############
%% Ampliacao do sistema para incluir integrador

% sistema aumentado
mzero = zeros(n,1);
Abarra = [A mzero; -C 0];
Bbarra = [B; 0];
Cbarra = [C 0];
Dbarra = 0;
sysMAaum = ss(Abarra, Bbarra, Cbarra, Dbarra);
sysMAaum

% quantidade de estados da planta aumentada
n_aum = 4

% polos do sistema aumentado
eig(sysMAaum)


%% Polos da Planta original

% tem qeu mudar a dimensao da matriz eye de acordo com A da planta
syms s;
I = eye(n);
eqn = det((s * I - A));
polos = vpasolve(eqn, s);
disp(polos);

% outra maneira de fazer para achar os autovalores de A
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

M = [Bbarra Abarra*Bbarra (Abarra^2)*Bbarra (Abarra^3)*Bbarra];
det(M)       

M = ctrb(Abarra,Bbarra);
rank(M)             

% soma dos elementos da matriz de contrabilidade
% sum(sum(M))
sum(M, 'all')

%% Aplicar a Formula de Ackermann K = [0 0 0 ... 1][B | AB | ... | A^(n-1)B] phi(A)
% M = [B | AB | ... | A^(n-1)B]^1
% onde phi(s) o polinomio caracteristico formado a partir dos polos desejados
% de malha fechada (s1 , s2 , ... , sn )
% onde  phi(s) s foi substituido por A, ficando phi(A)
s1 = -2;
s2 = -2;
s3 = -20;
s4 = -20;
polos_desejados = [s1 s2 s3 s4];

%Kbarra = [K | -ki], onde o K sao referentes aos polos e ki do integrador
Kbarra = acker(Abarra, Bbarra, polos_desejados);
disp(vpa(Kbarra));
% Lembrar que devido a insersao do polo adicional, ele vem com valor
% negativo mas colocar positivo para responder no moodle.

%% Representacao em malha fechada do sistema aumentado controlado em MF

%Sistema aumentado controlado em MF com integrador
Amfca = Abarra-Bbarra*Kbarra;
Bmfca = [mzero; 1]; % Rbarra = [zero(n, 1), 1]
Cmfca = Cbarra;
Dmfca = Dbarra;
sysMFca = ss(Amfca, Bmfca, Cmfca, Dmfca);

GanhoCC = dcgain(sysMFca);

fprintf('O ganho do sistema aumentado em malha fechada\n')
disp(double(GanhoCC));

fprintf('Polos sistema aumentado em malha fechada\n')
eig(sysMFca)




%% ############## PROJETO DO OBSERVADOR ##############
%%  Verificar observabilidade do sistema
% encontrar n, onde n e o numero de estados 
% -> n = dim(x)
% montar a matriz de observabilidade N
% N = [C | CA | ... | CA^(n-1)] 

% Se
% posto(N) = n -> Sistema observavel
% det(N) != 0  -> Sistema observavel
% Nao a necessidade de verificar os dois, apenas 1 ja e valido

n2 = n;

%Matriz de observabilidade
% N = [C; C*A, C*A²];
N = obsv(A, C)
rank(N)     
det(N)

% soma dos elementos da matriz de observabilidade
% sum(sum(N))
sum(N, 'all')


%% Definir autovalores do observador
eig(Amfca)

% autovalores 10 vezes maior, estipulado pelo exercicio
u1 = 10 * (-2);
u2 = 10 * (-2);
u3 = 10 * (-20);


%% Calcular polinomio caracteristico desejado
syms s
phi_s = conv([1 -u1], [1 -u2])
phi_s = conv(phi_s, -u3)
poly2sym(phi_s, s)

% lembrar de matriz identidade multiplicar constante
I = eye(n2)
s = A;
phi_a = 200*s^2 + 8000*s + 80000*I

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

%% ############## REPRESENTACAO MF COM CONTROLADOR E OBSERVADOR ##############
% K vem do Kbarra = [K... | -ki], onde n representa a quantidade de estados
% da planta
K =  Kbarra(1:n)

Abarra - Bbarra*Kbarra
%deducao x = x aumentado, xo = x observado
%x' = Abarra*x + Bbarra * u + Rbarra * r    video 14 minuto 7:16
%sendo u = -Kbarra*xo
%x' = Abarra*x + Bbarra * (-Kbarra*x) + Rbarra * r 
%x' = Abarra*x - Bbarra*Kbarra*xo
%adicionando e subtraindo  Bbarra*Kbarra*x video 15 minuto 31:20
%x' = Abarra*x - Bbarra*Kbarra*xo + Bbarra*Kbarra*x - Bbarra*Kbarra*x
%x' = (Abarra - Bbarra*Kbarra)x + Bbarra*Kbarra(x - xo)

%sendo xo' = (A - KeC)xo + B*-Kbarra*xo + Ke*y
%xo' = (A - KeC -KbarraB) + Ke*y


%Aco = Amf-Bbarra*K; %Amfc e o mesmo
%Bco = Bmfp2



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



%% Erro em regime permanente das variações paramétricas

% Teorema do valor final
% E = R - Y = R - R*Gmf
% E = R(1 - Gmf)
% e(inf) = s*E = s*R*(1 - Gmf)
syms s
R = 1/s;
I = eye(n_aum);
E = R * (1 - (Cmfp2*(inv(s*I - Amfp2))*Bmfp2));
erro = limit(s * E, s, 0);
fprintf('Saida em regime permanente do sistema aumentado com variacao\n')
disp(double(erro))

%% Saida em regime permanente das variações paramétricas

syms s
I = eye(n_aum);
E = vpa((Cmfp2*(inv(s*I - Amfp2))*Bmfp2));

% R = Entrada;
% degrau = 1/s; rampa = 1/s^2 ; parabola = 1/s^3
R = 1/s;
erro = limit(s * R * E , s, 0);
fprintf('Erro em regime permanente do sistema aumentado com variacao\n')
disp(vpa(erro));
