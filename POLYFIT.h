//
//  POLYFIT.h
//  Created by Baoqiang Cao on 11/13/18.
//  Copyright © 2018 Baoqiang Cao. All rights reserved.
//

#ifndef POLYFIT_h
#define POLYFIT_h
#define STD_MIN_FLOAT 0.000000001
#define STD_PARALELL 2
#define POLYFIT_FAILURE -999999
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct datax {
    double px;
    int trace;
};

typedef struct datax datax;

int polyfit(int newx, int *x, int *y, int num, int degree);


void sample_equal_X(int *x, int *selectedID, int num, int degree) ;
//polyfit2:
float* norm_input(int *dat, int num);
int polyfit2(int newx, int *x, int *y, int num, int degree);
float determinant(float **a, int N);
int cofactor(float **a, int N, float);
void transpose_inverse(float *num, float *fac, int N);
int inverse(float **a, int N);
float** matrix_multi(float **a, int ax, int ay, float **b, int bx, int by);
int transpose(float **a, int N, int M, float **b);
void vandermode(float **vand, float *data, int num, int degree);
void free_p2p(float **p, int row);
float **matrix_zero_init(int row, int colm);

#endif /* POLYFIT_h */
