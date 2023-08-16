#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

int main(){
//ENTRADA
int m=32;
double x0=0;
double xf=1;

//VALORES AUXILIARES
int i;

//GERADOR DE MALHA
double h=(xf-x0)/m;
double x[m+1];

//cout << "Malha:\n";
for(i=0;i<m+1;i++){
    x[i]=x0+i*h;
//    cout << scientific << setprecision(4) << x[i] << '\n';
}

//PREPARA SISTEMA LINEAR
double sys_diag[m-1], sys_diag_sub[m-1], sys_diag_up[m-1], sys_x[m-1], sys_b[m-1];

//Escreve primeiro valor nos vetores
sys_diag[0]=-2;
sys_diag_sub[0]=0;
sys_diag_up[0]=1;
sys_b[0]=pow(h,2)*(-pow(M_PI,2)*sin(M_PI*x[1])) - (-pow(M_PI,2)*sin(M_PI*x0));

//Escreve ultimo valor nos vetores
sys_diag[m-2]=-2;
sys_diag_sub[m-2]=1;
sys_diag_up[m-2]=0;
sys_b[m-2]=pow(h,2)*(-pow(M_PI,2)*sin(M_PI*x[m-1])) - (-pow(M_PI,2)*sin(M_PI*xf));

//Escreve valores intermediarios nos vetores
for(i=1;i<m-2;i++){
    sys_diag[i]=-2;
    sys_diag_sub[i]=1;
    sys_diag_up[i]=1;
    sys_b[i]=pow(h,2)*(-pow(M_PI,2)*sin(M_PI*x[i+1]));
}

//Debug: Mostra na tela os valores
//cout << "Sistema:\n";
//for(i=0;i<m-1;i++){
//    cout << i << "  ";
//    cout << scientific << setprecision(4) << sys_diag[i] << "  " << sys_diag_sub[i] << "  " << sys_diag_up[i] << "  " << sys_b[i] << '\n';
//}

//RESOLVE O SISTEMA USANDO METODO DE THOMAS
double u[m+1];
double w;

//Ajusta coeficientes (perdendo os originais)
for(i=1;i<m-1;i++){
    w=sys_diag_sub[i]/sys_diag[i-1];
    sys_diag[i]=sys_diag[i]-w*sys_diag_up[i-1];
    sys_b[i]=sys_b[i]-w*sys_b[i-1];
}

//Resolve o sistema
u[m-1]=sys_b[m-2]/sys_diag[m-2];
for(i=m-3;i>-1;i--){
    u[i+1]=(sys_b[i]-sys_diag_up[i]*u[i+2])/sys_diag[i];
}

//Coloca os valores nas bordas
u[0]=sin(M_PI*x[0]);
u[m]=sin(M_PI*x[m]);


//SOLUCAO EXATA
double u_hat[m+1];

for(i=0;i<m+1;i++){
    u_hat[i]=sin(M_PI*x[i]);
}


//SALVA RESULTADO
ofstream saida;
saida.open ("saida.csv");

//cout << "RESULTADOS:\n";
//cout << "_____x_____  _____U_____  _____Û_____  ____|E|____\n";

for(i=0;i<m+1;i++){
    saida << scientific << setprecision(5) << x[i] << "," << u[i] << "," << u_hat[i] << "," << abs(u[i]-u_hat[i]) << '\n';
}

saida.close();

return 0;
}