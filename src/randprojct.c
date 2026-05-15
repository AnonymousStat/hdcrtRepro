#include <math.h> 		// required for sqrt(), fabs();
#include <stdio.h>  	// required for exit
#include <stdlib.h> 	// required for malloc(), free();
#include <string.h> 	// required for memcpy()
#include <R.h>
#include <Rinternals.h> // required for SEXP et.al.;
#include "testhb.h"

double _RPtest0(double *y, double *x, int p, int n){
	int i,j;
	double tmp, nume = 0.0, deno=0.0, Tn, *Q, *R;

	Q  = (double*)malloc(sizeof(double)*n*p);
	R  = (double*)malloc(sizeof(double)*p*p);


	QRDecompN(Q, R, x, n, p);

	for(i=0;i<n;i++){
		deno += y[i]*y[i];
	}
	for(j=0;j<p;j++){
		tmp = 0.0;
		for(i=0;i<n;i++){
			tmp += Q[j*n+i]*y[i];
		}
		nume += tmp*tmp;
	}
	deno -= nume;
	Tn = (n-p)*nume/p/deno;

	free(Q);
	free(R);

	return(Tn);

}

double _RPtest(double *y, double *x, double *Pk, int p, int n, int k1){
	int i,j,s;
	double tmp, nume = 0.0, deno=0.0, Tn, *Uk, *XP, *Q, *R;

	Uk = (double*)malloc(sizeof(double)*n*k1);
	XP = (double*)malloc(sizeof(double)*n*k1);
	Q  = (double*)malloc(sizeof(double)*n*k1);
	R  = (double*)malloc(sizeof(double)*k1*k1);

	for(j=0;j<k1;j++){
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(s=0;s<p;s++){
				tmp += x[s*n+i]*Pk[j*p+s];
			}
			XP[j*n+i] = tmp;
		}
	}


	for(j=0;j<k1;j++){
		tmp = 0.0;
		for(i=0;i<n;i++){
			tmp += XP[j*n+i];
		}
		tmp /= n;
		for(i=0;i<n;i++){
			Uk[j*n+i] = XP[j*n+i] - tmp;
		}
	}

	QRDecompN(Q, R, Uk, n, k1);

	tmp = 0.0;
	for(i=0;i<n;i++){
		tmp += y[i];
		deno += y[i]*y[i];
	}
	deno -= tmp*tmp/n;
	for(j=0;j<k1;j++){
		tmp = 0.0;
		for(i=0;i<n;i++){
			tmp += Q[j*n+i]*y[i];
		}
		nume += tmp*tmp;
	}
	deno -= nume;
	Tn = (n-k1-1)*nume/k1/deno;

	free(Uk);
	free(XP);
	free(Q);
	free(R);

	return(Tn);

}

SEXP RPtest0(SEXP X_, SEXP Y_, SEXP DIM_)
{
	// dimensions
	int *dims 		= INTEGER(DIM_);
	int n     		= dims[0];
	int p     		= dims[1];

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);

	// Outcome
	SEXP rMDD;
	PROTECT(rMDD 	= allocVector(REALSXP, 1));

	REAL(rMDD)[0] = _RPtest0(y, x, p, n);

	UNPROTECT(1);
	return rMDD;

}

SEXP RPtest(SEXP X_, SEXP Y_, SEXP Pk_, SEXP DIM_)
{
	// dimensions
	int *dims 		= INTEGER(DIM_);
	int n     		= dims[0];
	int p     		= dims[1];
	int q     		= dims[2];

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);
	double *Pk 		= REAL(Pk_);

	// Outcome
	SEXP rMDD;
	PROTECT(rMDD 	= allocVector(REALSXP, 1));

	REAL(rMDD)[0] = _RPtest(y, x, Pk, p, n, q);

	UNPROTECT(1);
	return rMDD;

}

