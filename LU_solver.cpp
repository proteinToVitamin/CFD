#include <iostream>
#include <string>
using namespace std;


int const Nx = 2;
int const Ny = 2;
int const Ncells = Nx*Ny; // number of cells

// 線型方程式Ax=bを解くための関数
// AをLU分解して，そのままAに入れる.ピボット選択はしていない.
void LU_decomp(double A[][Ncells]) {
    for (int k = 0; k < Ncells - 1; ++k) {
        double w = 1./A[k][k];
        for (int i = k+1; i < Ncells; ++i) {
            A[i][k] *= w;
            for (int j = k+1; j < Ncells; ++j){
                A[i][j] -= A[i][k]*A[k][j]; 
            }
        }
    }
}

//線型方程式LU*x=bを解いて，根をbに格納する.
void LU_solve(double const LU[][Ncells], double b[]) {
    // LUはLU分解済みの係数行列
    for (int k = 1; k < Ncells; ++k){
        for (int i = 0; i < k; ++i){
            b[k] -=LU[k][i]*b[i];
        }
    }
    for (int k = Ncells-1; k >= 0; --k){
        for (int i = k+1; i < Ncells; ++i){
            b[k] -=LU[k][i]*b[i];
        }
        b[k] /= LU[k][k];
    }
}

int main(){
    double A[Ncells][Ncells]={
        {1, 1, 1, 1},
        {3, 2, 1, 0},
        {27, 9, 3, 1},
        {27, 6, 1, 0}
    };
    double b[Ncells]={5, 0, 1, 0};

    LU_decomp(A);
    LU_solve(A,b);

    for(int i = 0; i<Ncells; i++){
        cout<<b[i]<<endl;
    }
    

}