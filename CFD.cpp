#include <iostream>

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

    // 同年生係数をレイノルズ数から計算
    double nu = U*Lx/Re;

    // 現在の時間tと時間発展
    double t = dt;
    while(t <= t_end){
        //step1：流速の予測値の計算
        for(int i = 1; i <= Nx − 1; ++i){// upの計算
            for(int j = 0; j <= Ny - 1; ++j){

            }
        }
        for(int i = 0; i <= Nx - 1; ++i){// vpの計算
            for(int j = 1; j <= Ny - 1; ++j)
        }
        //step2：圧力の修正量の計算
        //step3：圧力と流速の修正
        //計算結果の出力
        t+=dt;
    }
}