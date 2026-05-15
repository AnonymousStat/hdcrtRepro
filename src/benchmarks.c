#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

void _Gram(double *x, int *param, double *xbyx, double **xTx)
{
    // input: 
    // x in R^{n*p}, saved by row.
    // param = c(n, p).
    // xbyx is the upper triangle of x^T x without diagonal, saved by row.
    // xTx is the matrix x^T x, 
    // only saves the upper triangle without diagonal, saved by row.
    unsigned int n, p, i, j, k, c;
    n = param[0];
    p = param[1];
    c = 1;
    double tmp, *xi, *xj;

    for (i = 0; i < n-1; ++i) {
        xi = x + i*p; // point to i-th sample
        for (j = i+1; j < n; ++j) {
            xj = x + j*p; // point to j-th sample
            for (tmp = 0.0, k = 0; k < p; ++k) 
                tmp += xi[k]*xj[k];
            xbyx[c] = tmp;
            c++;
        }
    }

    c = 0;
    for (i = 0; i < n-1; ++i) {
        xTx[i] = xbyx + c;
        c += n-2-i;
    }
}

double _tr_sigma2_hat(double **xTx, int n)
{
    int i, j, k, l;
    double trsigma2, tmp, y1n, y2n, y3n;
    y1n = y2n = y3n = 0.0;

    for (i = 0; i < n-1; ++i) {
        for (j = i+1; j < n; ++j) {
            y1n += xTx[i][j] * xTx[i][j];
            for (k = j+1; k < n; ++k) {
                tmp = (xTx[i][j] + xTx[i][k])*xTx[j][k] + xTx[i][j]*xTx[i][k];
                y2n += tmp/(n-2);
                for (l = k+1; l < n; ++l) {
                    tmp = xTx[i][j]*xTx[k][l] + xTx[i][k]*xTx[j][l] + xTx[i][l]*xTx[j][k];
                    y3n += tmp*4/(n-2)/(n-3);
                }
            }
        }
    }

    trsigma2 = (y1n - y2n*2 + y3n)/(n*(n-1)/2);
    return trsigma2;
}

void centralization(double *x, int n, int p)
{
    int i, j;
    double *meanx, *xi;
    meanx =  (double*)malloc(sizeof(double) *p);
    for (i = 0; i < p; i++) meanx[i] = 0.0;

    for (i = 0; i < n; i++) {
        xi = x + i*p;
        for (j = 0; j < p; j++) meanx[j] += xi[j];
    }
    for (i = 0; i < p; i++) meanx[i] /= n;
    for (i = 0; i < n; i++) {
        xi = x + i*p;
        for (j = 0; j < p; j++) xi[j] -= meanx[j];
    }
    free(meanx);
}

double _ZC_F_test(double **xTx, double *Y, int n)
{
    // avoid memory overflow
    // Zhong and Chen 2011 JASA
    // input: 
    // xTx in R^{n*n}.
    // Y in R^{n}.
    int i, j, k, l;
    double tmp, T0, T = 0.0;
    double c1 = (n-2)*(n-3);

    for (i = 0; i < n-3; ++i) {
        for (j = i+1; j < n-2; ++j) {
            T0 = 0.0;
            for (k = j+1; k < n-1; ++k) {
                for (l = k+1; l < n; ++l) {
                    tmp = xTx[i][k] - xTx[i][l] - xTx[j][k] + xTx[j][l];
                    T0 += tmp * (Y[i]-Y[j]) * (Y[k] - Y[l]);
                    tmp = xTx[i][j] - xTx[i][l] - xTx[j][k] + xTx[k][l];
                    T0 += tmp * (Y[i]-Y[k]) * (Y[j] - Y[l]);
                    tmp = xTx[i][j] - xTx[i][k] - xTx[j][l] + xTx[k][l];
                    T0 += tmp * (Y[i]-Y[l]) * (Y[j] - Y[k]);
                }
            }
            T += T0/c1;
        }
    }

    T /= n*(n-1)/2.0;

    return T;
}

double _FZW_RBS_test(double **xTx, int n)
{
    // Feng, Zou and Wang 2013 EJS
    // input: 
    // xTx in R^{n*n}.
    // Y has ascending ordered, Y_i <= Y_j
    int i, j;
    double *R, T = 0.0;

    R = (double*)malloc(sizeof(double)*n);
    for (i = 0; i < n; ++i) R[i] = (i+1.0)/(n+1.0) - 0.5;

    for (i = 0; i < n-1; ++i) {
        for (j = i+1; j < n; ++j) {
            T += R[i]*R[j]*xTx[i][j];
        }
    }

    T /= n*(n-1.0)/24.0;
    free(R);
    return T;
}

double _CGZ_F_test(double **xTx, double *X, double *Y, int n, int p)
{
    // Cui, Guo and Zhong 2018 AOS
    // input: 
    // xTx in R^{n*n}.
    // x in R^{n*p}, saved by row.
    // Y in R^{n}.
    int i, j;
    double T = 0.0;

    double *Xbar, *XbarTXi, *Xnorm, *xi;
    double Xbarnorm = 0.0, Ybar = 0.0;
    double deltax, deltay;
    Xbar    = (double*)malloc(sizeof(double)*p);
    XbarTXi = (double*)malloc(sizeof(double)*n);
    Xnorm   = (double*)malloc(sizeof(double)*n);
    for (j = 0; j < p; ++j) Xbar[j] = 0.0;

    for (i = 0; i < n; ++i) {
        XbarTXi[i] = 0.0;
        Xnorm[i] = 0.0;
        Ybar += Y[i];
        xi = X + i*p; // point to i-th sample
        for (j = 0; j < p; ++j) {
            Xbar[j] += xi[j];
            Xnorm[i] += xi[j] * xi[j];
        }
        Xnorm[i] /= 2.0*n;
    }
    Ybar /= n;
    for (j = 0; j < p; ++j) {
        Xbar[j] /= n;
        Xbarnorm += Xbar[j] * Xbar[j];
    }
    for (i = 0; i < n; ++i) {
        xi = X + i*p; // point to i-th sample
        for (j = 0; j < p; ++j) XbarTXi[i] += xi[j] * Xbar[j];
    }

    for (i = 0; i < n-1; ++i) {
        for (j = i+1; j < n; ++j) {
            deltax = Xbarnorm + Xnorm[i] + Xnorm[j] - XbarTXi[i] - XbarTXi[j];
            deltax += xTx[i][j] * (n-1.0) / n;
            deltay = (Y[i]-Ybar)*(Y[j]-Ybar) + (Y[i]-Y[j])*(Y[i]-Y[j])/n/2.0;
            T += deltax * deltay;
        }
    }

    T /= (1.0-1.0/n)*(n-2)*(n-2)/2.0;

    return T;
}

double _sigma2_hat_h0(double *Y, int n) {
    double i, Ey = 0.0, Ey2 = 0.0;
    for ( i = 0; i < n; i++ )
    {
        Ey += *Y;
        Ey2 += (*Y)*(*Y);
        Y++;
    }
    Ey /= n;
    double sigma2_hat = Ey2/n - Ey*Ey;
    sigma2_hat /= 1.0-1.0/n;
    return(sigma2_hat);
}

SEXP _ZC_Test(SEXP _X, SEXP _Y, SEXP _Param)
{
    int n;
    n = INTEGER(_Param)[0];
    double *y;
    y = REAL(_Y);

    double *xbyx;
    double **xTx;
    xbyx =  (double*)malloc(sizeof(double) * (n*n-n+2)/2);
    xTx  = (double**)malloc(sizeof(double*)*(n-1));

    SEXP _Test, _Tr_sigma2, _sigma2;
    PROTECT(_Test = allocVector(REALSXP, 1));
    PROTECT(_Tr_sigma2 = allocVector(REALSXP, 1));
    PROTECT(_sigma2 = allocVector(REALSXP, 1));

    _Gram(REAL(_X), INTEGER(_Param), xbyx, xTx);
    REAL(_Test)[0] = _ZC_F_test(xTx, y, n);
    REAL(_Tr_sigma2)[0] = _tr_sigma2_hat(xTx, n);
    REAL(_sigma2)[0] = _sigma2_hat_h0(y, n);

    SEXP Result, R_names;
    char *names[3] = {"test", "tr_sigma", "sigma2"};
    PROTECT(Result    = allocVector(VECSXP,  3));
    PROTECT(R_names   = allocVector(STRSXP,  3));
    SET_STRING_ELT(R_names, 0,  mkChar(names[0]));
    SET_STRING_ELT(R_names, 1,  mkChar(names[1]));
    SET_STRING_ELT(R_names, 2,  mkChar(names[2]));
    SET_VECTOR_ELT(Result, 0, _Test);
    SET_VECTOR_ELT(Result, 1, _Tr_sigma2);
    SET_VECTOR_ELT(Result, 2, _sigma2);
    setAttrib(Result, R_NamesSymbol, R_names); 

    free(xbyx);
    free(xTx);
    UNPROTECT(5);
    return Result;
}

SEXP _CGZ_Test(SEXP _X, SEXP _Y, SEXP _Param)
{
    int n, p;
    n = INTEGER(_Param)[0];
    p = INTEGER(_Param)[1];
    double *y;
    y = REAL(_Y);

    double *xbyx;
    double **xTx;
    xbyx =  (double*)malloc(sizeof(double) * (n*n-n+2)/2);
    xTx  = (double**)malloc(sizeof(double*)*(n-1));

    SEXP _Test, _Tr_sigma2;
    PROTECT(_Test = allocVector(REALSXP, 1));
    PROTECT(_Tr_sigma2 = allocVector(REALSXP, 1));

    _Gram(REAL(_X), INTEGER(_Param), xbyx, xTx);
    REAL(_Test)[0] = _CGZ_F_test(xTx, REAL(_X), y, n, p);
    REAL(_Tr_sigma2)[0] = _tr_sigma2_hat(xTx, n);

    
    // Rprintf("Tn is %f.\n", REAL(_Test)[0]);
    // Rprintf("trSigma2 is %f.\n", REAL(_Tr_sigma2)[0]);

    SEXP Result, R_names;
    char *names[2] = {"test", "tr_sigma"};
    PROTECT(Result    = allocVector(VECSXP,  2));
    PROTECT(R_names   = allocVector(STRSXP,  2));
    SET_STRING_ELT(R_names, 0,  mkChar(names[0]));
    SET_STRING_ELT(R_names, 1,  mkChar(names[1]));
    SET_VECTOR_ELT(Result, 0, _Test);
    SET_VECTOR_ELT(Result, 1, _Tr_sigma2);
    setAttrib(Result, R_NamesSymbol, R_names); 

    free(xbyx);
    free(xTx);
    UNPROTECT(4);
    return Result;
}

SEXP _FZW_Test(SEXP _X, SEXP _Param)
{
    int n, p;
    n = INTEGER(_Param)[0];
    p = INTEGER(_Param)[1];

    double *xbyx;
    double **xTx;
    xbyx =  (double*)malloc(sizeof(double) * (n*n-n+2)/2);
    xTx  = (double**)malloc(sizeof(double*)*(n-1));

    SEXP _Test, _Tr_sigma2;
    PROTECT(_Test = allocVector(REALSXP, 1));
    PROTECT(_Tr_sigma2 = allocVector(REALSXP, 1));

    // FZW's method needs to do centralization first
    centralization(REAL(_X), n, p);
    _Gram(REAL(_X), INTEGER(_Param), xbyx, xTx);
    REAL(_Tr_sigma2)[0] = _tr_sigma2_hat(xTx, n);
    REAL(_Test)[0] = _FZW_RBS_test(xTx, n);

    SEXP Result, R_names;
    char *names[2] = {"test", "tr_sigma"};
    PROTECT(Result    = allocVector(VECSXP,  2));
    PROTECT(R_names   = allocVector(STRSXP,  2));
    SET_STRING_ELT(R_names, 0,  mkChar(names[0]));
    SET_STRING_ELT(R_names, 1,  mkChar(names[1]));
    SET_VECTOR_ELT(Result, 0, _Test);
    SET_VECTOR_ELT(Result, 1, _Tr_sigma2);
    setAttrib(Result, R_NamesSymbol, R_names); 

    free(xbyx);
    free(xTx);
    UNPROTECT(4);
    return Result;
}


