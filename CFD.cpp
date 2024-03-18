    for (int i = 0; i <= Nx - 1; ++i){
        for (int j = 0; j <= Ny - 1; ++j) {
            double m = M[i][j];
            M[i][j] = -m;
        }
    }
