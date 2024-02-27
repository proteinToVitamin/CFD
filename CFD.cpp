#include <iostream>
#include <string>
using namespace std;

// 計算空間のサイズと天井の速度 
double const Lx = 1.;
double const Ly = 1.;
double const U = 1.;
// Reynolds 数
double const Re = 1.;
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

int main(){
    // セルの幅
    double dx = Lx/Nx;
    double dy = Ly/Ny;

    // 動粘性係数をレイノルズ数から計算
    double nu = U*Lx/Re;

    // 現在の時間tと時間発展
    double t = dt;
    while(t <= t_end){
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
                vn = 0.; } 
                else {
                un = 0.5*(u[i][j] + u[i][j+1]);
                vn = 0.5*(v[i-1][j+1] + v[i][j+1]); 
                }
                double uvy = (un*vn-us*vs)/dy;

                double dpdx = (p[i][j] - p[i-1][j])/dx;

                double uxx = (u[i-1][j]-2*u[i][j] + u[i+1][j])/(dx*dx); 
                double uyy ;
                if (0 == j){
                    double us = 0.;
                    uyy = (8/3*us - (8/3+4/3)*u[i][j] + 4/3*u[i][j +1])/(dy*dy);
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
        //step3：圧力と流速の修正
        //計算結果の出力
        t+=dt;
    }
}
