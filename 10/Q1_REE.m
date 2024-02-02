close all;
clear all;
clc;

%% Declara��o do sistema
A=[0 1;-6 -5];
B=[0; 1];
C=[1 0];
D=0;
sys=ss(A,B,C,D);


%% Polos
syms s;

I = eye(2);
eqn = det((s * I - A));
polos = vpasolve(eqn, s);
disp(polos);

% retorna os autovalores que � o mesmo que os polos
% disp(eig(A))

%% C�lculo da matriz de transi��o de estados
syms s;

% matriz identidade
I = eye(2);
phi_s = inv(s * I - A);

% isolar phi
disp(phi_s * eqn);

% inversa de laplace
eat = ilaplace(phi_s);

% outra forma de obter eat 
% syms t
% expm(A*t)

% Calculo do residuo se necessario
% num = [1 5]
% den = [1 5 6]
% residue(num, den)

%% ----------------------- Equa��o de estado -----------------------------

%% Integral
% x(t) = e^At * x(0) + integral e^A(t-tal) * B * entrada(tal) dtal

% Integral
syms tal t
eatal = eat * B;

% pegar saidas
% substituir t por (t - tal)
% resolver integral

% separada as integrais int = int(func1) + int(func2)
func1 = exp(-2*(t - tal)) 
func2 = - exp(-3*(t - tal))

func3 = 3*exp(-3*(t - tal)) 
func4 = -2*exp(-2*(t - tal))

% int(func, variavel, inicio, fim)
I1 = int(func1, tal, 0, t)
I2 = int(func2, tal, 0, t)

I3 = int(func3, tal, 0, t)
I4 = int(func4, tal, 0, t)

INT1 = I1 + I2
INT2 = I3 + I4


%% Estruturando saida

f1 = INT1
vpa(f1)
f2 = INT2
vpa(f2)
f = [f1; f2]

syms x1_0 x2_0

x0 = [x1_0; x2_0];
% x = [x1; x2]
x = eat * x0 + f

x1(x1_0, x2_0) = x(1, 1);
x2(x1_0, x2_0) = x(2, 1);

%% Estruturando saida y = C*x em funcao do tempo
syms t
% Y e obtido a EE inicial
% y(t) = C * Estados
syms e1 e2
C * [e1; e2]

% e1 = x1 e e2 = x2

% estruturando y com base em x1 com valores aplicados
x1_0 = 0
x2_0 = 0

y1(t) = x1(x1_0, x2_0)
y2(t) = x2(x1_0, x2_0)

double(y(1))
double(y(inf))

double(limit(y1,t,1))
double(limit(y1,t,inf))

double(limit(y2,t,1))
double(limit(y2,t,inf))
