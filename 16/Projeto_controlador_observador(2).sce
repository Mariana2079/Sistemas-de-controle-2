//Matrizes da planta
A=[0 1 0;0 0 1;-6 -11 -6];
B=[0 0 1]';
C=[1 0 0];
D=0;

//Polos da planta
spec(A)

/****************************************
*         Projeto do Controlador        *
*****************************************/

//Sistema aumentado
Abar=[A zeros(3,1);-C 0];
Bbar=[B;0];
Cbar=[C 0];
Dbar=[0];

//Matriz de controlabilidade do sistema aumentado
M=cont_mat(Abar,Bbar)
//Posto da matriz de controlabilidade do sistema aumentado
rank(M)

//Polos de MF desejados
polosMF=[-2+2*sqrt(3)*%i -2-2*sqrt(3)*%i -10 -50];
//Projeto do vetor de ganhos Kbar do controlador
Kbar=ppol(Abar,Bbar,polosMF)

//Verificação do projeto do controlador
spec(Abar-Bbar*Kbar)


/****************************************
*         Projeto do Observador         *
*****************************************/

//Matriz de observabilidade
N=obsv_mat(A,C)
//Posto da matriz de observabilidade
rank(N)

//Autovalores desejados para o observador
autovaloresOBSV=[-10 -10 -50];
Ke=ppol(A',C',autovaloresOBSV)'

//Verificação do projeto do observador
spec(A-Ke*C)
