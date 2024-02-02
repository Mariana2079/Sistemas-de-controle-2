close all;
clear all;
clc;

%% Declaração do sistema
A = [0 1 0; 0 0 1; -6 -11 -6];
B = [0; 0; 1];
C = [1 0 0];
D = 0;
sysMA = ss(A, B, C, D);
sysMA

%% Polos ,  autovalores doa martiz
disp(eig(A));


%% Verificar observabilidade do sistema
% encontrar n, onde n e o numero de estados 
% -> n = dim(x) 
% montar a matriz de observabilidade N
% N = [C | CA | ... | CA^(n-1)] 

% Se
% posto(N) = n -> Sistema observavel
% det(N) != 0  -> Sistema observavel
% Nao a necessidade de verificar os dois, apenas 1 ja e valido

n = 3;

%Matriz de observabilidade
N=obsv(A,C)

N = [C; C*A; C*A^2];
det(N)      

N = ctrb(A,B);
rank(N)           


%% Definir autovalores do observador
eig(A)

% autovalor mais rapido
aut = -3;
% O -3 e o autovalor mais rapido da planta, por isso 5 * autovalor
u1 = 5 * aut;
u2 = u1
u3 = u2


%% Calcular polinomio caracteristico desejado


syms s
phi_s = conv(conv([1 -u1], [1 -u2] ), [1 -u3])
poly2sym(phi_s, s)

% lembrar de matriz identidade multiplicar constante
I = eye(3)
s = A;
phi_a = s^3 + 45*s^2 + 675*s + 3375*I

%% Aplicar a Formula de Ackermann Ke = phi(A)N [0; 0; ...0; 1]
% N = [C | CA | ... | CA^(n-1)] ^-1
% onde phi(s) o polinomio caracteristico formado a partir dos polos desejados
% de malha fechada (s1 , s2 , ... , sn )
% onde  phi(s) s foi substituido por A, ficando phi(A)


Ke = phi_a * inv(N) * [0; 0; 1]


% todos os passo anteriores se resumem no acker
u1 = 5 * -3;
u2 = u1
u3 = u2

autovaloresOBSV = [u1 u2 u3];
Ke = acker(A', C', autovaloresOBSV)';


%% Representacao do observador
% xtil' = (A - Ke C)xtil + Bu + Ke y
% xtil' = (A - Ke C)xtil + [B Ke] [u; y]
% ytil = I xtil

I = eye(3)

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

% como verificado os autovalores deram iguais a 15, validando o resultado
% esperado
