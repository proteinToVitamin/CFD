#include <iostream>
#include <string>
#include <cmath>
using namespace std;


int const Nx = 2;
int const Ny = 2;
int const Ncells = Nx*Ny; // number of cells


void printMatrix(double A[Ncells]) {
    for (int i = 0; i < Ncells; ++i) {
            cout << A[i]<<endl;
        }
    cout << endl;
}

// 線型方程式Ax=bを解くための関数
void jacobi(double A[][Ncells], double b[]){
    //更新前後のxを保存する。
    //xの初期値は全要素0とする
    double x_new[Ncells];
    double x_old[Ncells];
    for (int i = 0; i < Ncells; ++i){
        x_old[i]=b[i];
        x_new[i]=0;
    }

    //収束条件
    const double delta = 0.0000001;
    int iteration=0;
    while (true){
        iteration+=1;
#ifdef DEBUG
        cout<<"_______________"<<endl;
        cout<<"iteration : "<<iteration<<endl;
        cout<<endl;
        cout<<"updated x is"<<endl;
        cout<<endl;
        printMatrix(x_new);
#endif
        //xを更新する
        for (int i = 0; i < Ncells; ++i){
            double sum=0.0;
            for (int j = 0; j < Ncells; ++j){
                if(j != i){
                    sum+=A[i][j]*x_old[j];
                }
            } 
            x_new[i]=(b[i]-sum)/A[i][i];
        }

        //ステップ進展による差を計算
        double error=0;
        for (int k = 0; k < Ncells; ++k){
            error+=(x_new[k]-x_old[k])*(x_new[k]-x_old[k]);
        }
#ifdef DEBUG
        cout<<"error : "<<error<<endl;
        cout<<"delta : "<<delta<<endl;
        cout<<endl;
#endif

        //収束判定
        if (error<delta){
            cout<<"error is small enough"<<endl;
            break;
        }
        else if(isnan(error)){
            cout<<"error is not defined"<<endl;
            break;
        }
        else if(isinf(error)){
            cout<<"error is infinite"<<endl;
            break;
        }
        else{
            for(int h = 0; h < Ncells; ++h){
                x_old[h]=x_new[h];
            }
        }
    }
    //収束した後
    for(int h = 0; h < Ncells; ++h){
                b[h]=x_new[h];
            }
}



int main(){
    double A[Ncells][Ncells]={
        {10, 1, 1, 1},
        {3, 20, 1, 0},
        {3, 9, 30, 1},
        {27, 6, 1, 40}
    };
    double b[Ncells]={5, 0, 1, 0};

    jacobi(A,b);
    cout<<"answer is"<<endl;
    cout<<endl;
    for(int i = 0; i<Ncells; i++){
        cout<<b[i]<<endl;
    }
    cout<<endl;
}