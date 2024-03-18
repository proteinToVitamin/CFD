#include <iostream>
#include <string>
#include <cmath> 


using namespace std;

// 計算空間のサイズと天井の速度 
double const Lx = 1.;
double const Ly = 1.;
double const U = 1.;
// Reynolds 数
double const Re = 1;
// 計算空間の分割数
int const Nx = 16;
int const Ny = 16;
int const Ncells = Nx*Ny; // number of cells

// 流れの物理量 
double u[Nx+1][Ny];
double v[Nx][Ny+1]; 
double p[Nx][Ny];
// SMAC法の予測値 
double up[Nx+1][Ny]; 
double vp[Nx][Ny+1];
// 圧力の修正量 
double phi[Nx][Ny];
// 解析時間と時間刻み 
double t_end = 10.; 
double dt = 0.0001;

// 係数行列
double M[ Ncells ][ Ncells ];
// 右辺ベクトルを格納する配列.線型方程式を解き終わった時には，根が格納される. 
double b[Ncells];


// 線型方程式を解くための関数

//LU分解
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

//--------------------//
//Jacobi法
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
#ifdef DEBUG_JACOBI
        cout<<"_______________"<<endl;
        cout<<"iteration : "<<iteration<<endl;
        cout<<endl;
        cout<<"updated x is"<<endl;
        cout<<endl;
        //printMatrix(x_new);
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
#ifdef DEBUG_JACOBI
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

#ifdef DEBUG
// 行列・ベクトル積 b = A*x
void Matvec(double const A[][Ncells], double const x[], double b[]) {
    for (int i = 0; i < Ncells; ++i) { 
        b[i] = 0.;
        for (int j = 0; j < Ncells; ++j) 
            b[i] += A[i][j]*x[j];
    } 
    cout<<"DEBUG was conducted"<<endl;
}
#endif // DEBUG

int main(){
    // セルの幅
    double dx = Lx/Nx;
    double dy = Ly/Ny;

    // 動粘性係数をレイノルズ数から計算
    double nu = U*Lx/Re;

    // 現在の時間tと時間発展
    double t = dt;

    // Poisson方程式の係数行列を作る 
    // Mは0で初期化されている前提
    for (int i = 0; i <= Nx - 1; ++i)
        for (int j = 0; j <= Ny - 1; ++j) {
            int P=i+j*Nx; 
            int S = i +(j-1)*Nx; 
            int W=(i-1)+j *Nx; 
            int E=(i+1)+j *Nx; 
            int N=i +( j +1)*Nx; 
            if (j >= 1) {
                M[P][S] = 1./(dy*dy);
                M[P][P] -= 1./(dy*dy); 
                }
            if (i >= 1) {
                M[P][W] = 1./(dx*dx); 
                M[P][P] -= 1./(dx*dx);
                }
            if (i <= Nx - 2) {
                M[P][E] = 1./(dx*dx);
                M[P][P] -= 1./(dx*dx); 
            }
            if (j <= Ny - 2) { 
                M[P][N] = 1./(dy*dy); 
                M[P][P] -= 1./(dy*dy);
            } 
        }
    
#ifdef DEBUG
    double OrigM[Ncells][Ncells];
    for (int i = 0; i <= Nx*Ny - 1; ++i)
        for (int j = 0; j <= Nx*Ny - 1; ++j) 
            OrigM[i][j] =M[i][j];

#endif // DEBUG

#if 1
    LU_decomp(M) ;
#endif

    while(t <= t_end){
        std::cout << "Time: " << t << std::endl; 
        std::cerr << "Time: " << t << std::endl;

        //step1：流速の予測値の計算
        for(int i = 1; i <= Nx - 1; ++i){// upの計算
            for(int j = 0; j <= Ny - 1; ++j){
                
                double uw = (u[i-1][j]+u[i][j])/2.; 
                double ue = (u[i][j]+u[i][j+1])/2.;

                double uux = ((ue*ue)-(uw*uw))/dx;

                double us, vs, un, vn;
                if (0 == j) { 
                    us = vs = 0.;
                } 
                else {
                us = 0.5*(u[i][j-1] + u[i][j]); vs = 0.5*(v[i-1][j] + v[i][j]);
                }
                if (Ny-1 == j) {
                un = U;
                vn = 0.; 
                } 
                else {
                un = 0.5*(u[i][j] + u[i][j+1]);
                vn = 0.5*(v[i-1][j+1] + v[i][j+1]); 
                }
                double uvy = (un*vn-us*vs)/dy;

                // 圧力勾配項
                double dpdx = (p[i][j] - p[i-1][j])/dx;

                //　粘性項
                double uxx = (u[i-1][j]-2*u[i][j] + u[i+1][j])/(dx*dx); 
                double uyy ;
                if (0 == j){
                    double us = 0.;
                    uyy = (8/3*us - (8/3+4/3)*u[i][j] + 4/3*u[i][j +1])/(dy*dy);
                }else if (Ny - 1 == j) {
                    double un = U;
                    uyy = (4./3.*u[ i ][ j -1] - (4./3.+8./3.)*u[ i ][ j ] + 8./3.*un)/(dy*dy);
                }else{
                    uyy = (u[i][j-1] - 2*u[i][j] + u[i][j+1])/(dy*dy);
                }
                
                up[i][j] = u[i][j] + dt*(-(uux+uvy) -dpdx + nu*(uxx+uyy));
            }
        }

        for(int i = 0; i <= Nx - 1; ++i){// vpの計算
            for(int j = 1; j <= Ny - 1; ++j){
                //対流項x
                double uw, vw, ue, ve;
                if (0 == i) {
                    uw = vw = 0;
                }
                else {
                    uw = 0.5*(u[i][j-1] + u[i][j]); vw = 0.5*(v[i-1][j] + v[i][j]);
                }
                if (Nx - 1 == i) 
                { ue = ve = 0;
                } 
                else {
                    ue = 0.5*(u[i+1][j-1] + u[i+1][j]); ve = 0.5*(v[i][j] + v[i+1][j]);
                    }
                    double vux = (-(vw*uw)+(ve*ue))/dx;
                //対流項y
                double vs = 0.5*(v[i][j-1] + v[i][j]); 
                double vn = 0.5*(v[i][j] + v[i][j+1]);
                double vvy = (-vs*vs + vn*vn)/dy; 

                // 圧力勾配項 y
                double dpdy = (-p[i][j-1] + p[i][j])/dy;
                //粘性項y 
                double vxx ;
                if (0 == i) {
                double vw = 0.;
                vxx = (8./3.*vw - (8./3.+4./3)*v[i][j] + 4./3.*v[i+1][j])/(dx*
                dx);
                } 
                else if (Nx - 1 == i) {
                double ve = 0;
                vxx = (4./3.*v[i-1][j] - (8./3.+4./3)*v[i][j] + 8./3.*ve)/(dx*
                dx); } 
                else{
                     vxx = (v[i-1][j] - 2*v[i][j] + v[i+1][j])/(dx*dx);
                }

                double vyy = (v[i][j-1] - 2*v[i][j] + v[i][j+1])/(dy*dy);

                vp[i][j] = v[i][j] + dt*(-(vux+vvy) -dpdy + nu*(vxx+vyy));
            }
        }
        //step2：圧力の修正量の計算
        for (int i = 0; i <= Nx - 1; ++i)
            for (int j = 0; j <= Ny - 1; ++j)
                b[i+j*Nx] = ((-up[i][j]+up[i+1][j])/dx + (-vp[i][j]+vp[i][j+1])/dy ) / dt ;
        
#ifdef DEBUG
        double RHS[Ncells];
        for (int i = 0; i < Ncells; ++i)
            RHS[i] = b[i];
#endif // DEBUG

#if 1
        LU_solve(M, b);
#else
        for (int i = 0; i <= Ncells - 1; ++i){
            double x = b[i];
            b[i] = -x;
            for (int j = 0; j <= Ncells - 1; ++j) {
                double m = M[i][j];
                M[i][j] = -m;
            }
        }
        jacobi(M,b);
#endif

#ifdef DEBUG
        // 根のテスト
        double prod [ Ncells ] ; 
        Matvec(OrigM, b, prod);
        double sum_diff = 0;
        for (int i = 0; i < Ncells; ++i)
            sum_diff += (RHS[i]-prod[i])*(RHS[i]-prod[i]);
        std::cerr << "RMS of Diff == " << sqrt(sum_diff / Ncells) << std::endl ;
#endif // DEBUG

        for (int i = 0; i <= Nx - 1; ++i) 
            for (int j = 0; j <= Ny - 1; ++j)
                phi[i][j] = b[i+j*Nx];
        //step3：圧力と流速の修正
        for (int i = 0; i <= Nx - 1; ++i) 
            for (int j = 0; j <= Ny - 1; ++j)
                p[i][j] += phi[i][j];
        for (int i = 1; i <= Nx - 1; ++i) 
            for (int j = 0; j <= Ny - 1; ++j)
                u[i][j] = up[i][j] - dt*(-phi[i-1][j] + phi[i][j])/dx;
        for (int i = 0; i <= Nx - 1; ++i) 
            for (int j = 1; j <= Ny - 1; ++j)
                v[i][j] = vp[i][j] - dt*(-phi[i][j-1] + phi[i][j])/dy;

        //計算結果の出力
        //P
        cout << "P" << endl; 
        for (int i = 0; i <= Nx - 1; ++i)
            for (int j = 0; j <= Ny - 1; ++j)
                cout << p[i][j] << endl;
        //U
        cout << "U" << endl; 
        for (int i = 1; i <= Nx - 1; ++i)
            for (int j = 0; j <= Ny - 1; ++j) 
                cout << u[i][j] << endl;
        //V
        cout << "V" << endl; 
        for (int i = 0; i <= Nx - 1; ++i)
            for (int j = 1; j <= Ny - 1; ++j) 
                cout << v[i][j] << endl;

#ifdef DEBUG
        // 連続の式に対する誤差のチェック
        double sum_voldiff = 0;
        for (int i = 0; i <= Nx - 1; ++i)
            for (int j = 0; j <= Ny - 1; ++j) {
                double voldiff = (u[i][j] - u[i+1][j])*dy + (v[i][j] - v[i][j +1])*dx;
                sum_voldiff += voldiff * voldiff ; 
                }
        cerr<<"RMS error of Continuity: "<<sqrt (sum_voldiff/Ncells)<<endl;
#endif

        //時間の更新
        t+=dt;
    }
}
