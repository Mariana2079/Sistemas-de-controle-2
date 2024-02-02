clear;
clc;

//Matrizes da planta
A=[0 1 0;0 0 1;-6 -11 -6];
B=[0 0 1]';
C=[1 0 0];

//Polos da planta
spec(A)

//Matriz de observabilidade
N=obsv_mat(A,C)
//Posto da matriz de observabilidade
rank(N)

//Vetor de ganhos do observador para autovalores -15, -15, -15
autovaloresOBSV=[-15 -15 -15];
Ke=ppol(A',C',autovaloresOBSV)'

//Verificação dos autovalores do observador
spec(A-Ke*C)

//Vetor de ganhos do observador para autovalores -30, -30, -30
autovaloresOBSV=[-30 -30 -30];
Ke2=ppol(A',C',autovaloresOBSV)'

//Verificação dos autovalores do observador
spec(A-Ke2*C)
