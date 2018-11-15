//
//  POLYFIT.c
//  Created by Baoqiang Cao on 11/13/18.
//  Copyright © 2018 Baoqiang Cao. All rights reserved.
//

#include "POLYFIT.h"

/*
 *polyfit: a + b x + c x^2 + d x^3 + e x^4
 Reference is https://en.wikipedia.org/wiki/Polynomial_interpolation
 I use https://en.wikipedia.org/wiki/Polynomial_regression to solve it.
 */



int inverse(float **a, int N) {
    //Assume a is a N by N matrix
    float d = determinant(a, N);
    if (d == 0) {
        printf("\nERROR: Inverse of Entered Matrix is not possible\n");
        return -1;
    } else {
        printf("Determinant is %.4f\n", d);
        return(cofactor(a, N, d));
    }
}

/*For calculating Determinant of the Matrix */
float determinant(float **a, int k) {
    float det = 0.0;
    float **b = NULL;
    if (k < 0) {
        printf("Erro, k should never be negative!\n");
    } else if (k == 1) {
        det = a[0][0];
    } else if (k == 2) {
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    } else {
        int j, m, n, jj;
        for (j = 0; j < k; j++) {
            b = matrix_zero_init(k-1, k-1);
            if (b != NULL) {
                for (m = 1; m < k; m++) {
                    for(n = 0, jj=0; n < k; n++) {
                        if (n != j) {
                            b[m-1][jj] = a[m][n];
                            jj++;
                        }
                    }
                }
                det = det + pow(-1, j) * a[0][j] * determinant(b, k - 1);
                free_p2p(b, k-1);
            } else {
                printf("Error in matrix_zero_init for b, k=%d\n", k-1);
            }
        }
    }
    return (det);
}

int cofactor(float **num, int k, float Det) {
    float **b = matrix_zero_init(k-1, k-1); //calloc(k*k, sizeof(float));
    float **fac = matrix_zero_init(k, k); //calloc(k*k, sizeof(float));
    int i, j, m, n, bi, bj;
    int status = -1;
    if (fac != NULL && b != NULL) {
        for(i = 0; i < k; i++) {
            for (j = 0; j < k; j++) {
                fac[i][j] = 0.0;
                for (m = 0, bi=0; m < k; m++) {
                    if (m != i) {
                        for (n = 0, bj = 0; n < k; n++) {
                            if (n != j) {
                                b[bi][bj] = num[m][n];
                                bj++;
                            }
                        }
                        bi++;
                    }
                }
                fac[i][j] = pow(-1, i+j) * determinant(b, k-1);
            }
        }
        //transpose_inverse(num, fac, k);
        for (i = 0; i < k; i++) {
            for (j = 0; j < k; j++) {
                num[i][j] = fac[i][j] / Det;
            }
        }
        free_p2p(b, k-1);
        free_p2p(fac, k);
        status = 1;
    }
    return(status);
}

int transpose(float **a, int N, int M, float **b) {
    int i, j, success = 1;
    if (a != NULL && b != NULL) {
        for(i = 0; i < N; i++) {
            for(j = 0; j< M; j++)
                b[j][i] = a[i][j];
        }
    } else {
        success = -1;
    }
    return (success);
}

float** matrix_multi(float **a, int ax, int ay, float **b, int bx, int by) {
    int i, j, k;
    float **prod = matrix_zero_init(ax, by);
    if (ay == bx && prod != NULL) {
        
            for (i = 0; i < ax; i++) {
                for (j = 0; j < by; j++) {
                    for (k = 0; k < ay; k++) {
                        prod[i][j] += a[i][k] * b[k][j];
                    }
                }
            }
        
    } else {
        if (prod != NULL) {
            free_p2p(prod, ax);
        }
    }
    return (prod);
}

void free_p2p(float **p, int row) {
    for(int i=0; i < row; i++) free(p[i]);
    free(p);
}

void vandermode(float **vand, float *data, int num, int degree) {
    /*
     *data is num * int allocated and assigned already.
     degree > 0;
     */
    int i, j;
    for (i = 0; i < num; i++) {
        //vand[i] = (float *) malloc((degree+1) * sizeof(float));
        vand[i][0] = 1.0;
        for (j = 1; j < (degree+1); j++) {
            vand[i][j] = pow(*(data+i), j);
        }
    }
}

float **matrix_zero_init(int row, int colm) {
    float **mat = NULL;
    if ( (mat = calloc(row, sizeof(float *))) != NULL) {
        for (int i = 0; i < row; i++) {
            if ( (mat[i] = calloc(colm, sizeof(float))) == NULL ) {
                free_p2p(mat, i-1);
                break;
            }
        }
    }
    return(mat);
}

float* norm_input(int *dat, int num) {
    float* ret = calloc(num+2, sizeof(float));
    float mean, std;
    int i;
    for(i = 0, mean=0.0; i < num; i++) {
        mean += *(dat+i);
    }
    mean /= num;
    if (num == 1) {
        //std is undefined, so just give the minimum float, in my definintion the minimum.
        std = STD_MIN_FLOAT;
    } else {
        for(i = 0, std=0.0; i < num; i++) {
            std += pow(*(dat+i)-mean, 2);
        }
        std = sqrt(std /(num-1));
    }
    for (i =0; i < num; i++) {
        *(ret+i) = (*(dat+i) - mean) / std;
    }
    *(ret + num) = mean;
    *(ret + num + 1) = std;
    return(ret);
}

int polyfit2(int newx, int *x, int *y, int num, int degree) {
    /*Will not do check here, so before call this method make sure:
     1. x and y are not null.
     2. num > degree
     change to the new strategy: https://en.wikipedia.org/wiki/Polynomial_regression
     Y = B X + E
     B = (X_t X)_inv X_t Y
     */
    double ret = 0.0f;
    
    float *nx = norm_input(x, num);
    float *ny = norm_input(y, num);
    if ( *(ny+num+1) < STD_PARALELL ) {
        //this means, the line is horizontal, we treat this as non-exist.
        return(POLYFIT_FAILURE);
    }
    
    float **vand = matrix_zero_init(num, degree+1);
    float **t_vand = matrix_zero_init(degree+1, num);
    
    if (vand == NULL || t_vand == NULL) {
        return(POLYFIT_FAILURE);
    }
    vandermode(vand, nx, num, degree);
    transpose(vand, num, degree+1, t_vand);
    float **pny = matrix_zero_init(num, 1);

    for (int i = 0; i < num; i++) {
        pny[i][0] = *(ny+i);
    }
  
    float **mprod = matrix_multi(t_vand, degree+1, num, vand, num, degree+1);
    
    int flag_inv = inverse(mprod, (degree+1));
    if (flag_inv == 1) {
        float **prod_x = matrix_multi(mprod, degree+1, degree+1, t_vand, degree+1, num);
        float **B = matrix_multi(prod_x, degree+1, num, pny, num, 1);
        if (B!= NULL) {
            for (int i = 0; i < (degree+1); i++) {
                ret += B[i][0] * pow( (newx- *(nx+num)) / *(nx+num+1), i);
                //printf("B:i=%d, Bi=%.4f\n", i, B[i][0]);
            }
            ret = *(ny+num+1) * ret + *(ny + num);
            //printf("newx=%d, ret=%d, avgy=%.1f, stdy=%.1f, avgx=%.1f, stdx=%.1f\n", newx, (int) ret, *(ny+num), *(ny+num+1), *(nx+num), *(nx+num+1));
        } else {
            ret = POLYFIT_FAILURE;
        }
        free_p2p(B, degree+1);
        free_p2p(prod_x, degree+1);
    } else {
        ret = POLYFIT_FAILURE;
    }
    free(nx);
    free(ny);
    free_p2p(mprod, degree+1);
    free_p2p(vand, num);
    free_p2p(t_vand, degree+1);
    free_p2p(pny, num);
    return (int)ret;
}
