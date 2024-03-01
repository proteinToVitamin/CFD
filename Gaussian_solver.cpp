#include <iostream>
#include <cmath>

int main() {
    const int N = 4;
    double A[N][N] = {{1, 1, 1, 1},
                       {3, 2, 1, 0},
                       {27, 9, 3, 1},
                       {27, 6, 1, 0}
                       };
    double B[N] = {5, 0, 1, 0};
    double X[N]; // 解の配列

    // ガウスの消去法で連立方程式を解く
    for (int i = 0; i < N - 1; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double ratio = A[j][i] / A[i][i];
            for (int k = i; k < N; ++k) {
                A[j][k] -= ratio * A[i][k];
            }
            B[j] -= ratio * B[i];
        }
    }

    // 後退代入で解を求める
    for (int i = N - 1; i >= 0; --i) {
        X[i] = B[i] / A[i][i];
        for (int j = i - 1; j >= 0; --j) {
            B[j] -= A[j][i] * X[i];
        }
    }

    // 解の出力
    for (int i = 0; i < N; ++i) {
        std::cout << "x[" << i << "] = " << X[i] << std::endl;
    }

    return 0;
}
