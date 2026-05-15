#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "testhb.h"

double gradient_wts[16] = {1.0,0.37267800,0.22395160,0.11938701,0.06250000,0.03125000,0.01562500,
	0.00781250,0.00390625,0.00195313,0.00097656,0.00004883,0.00000244,0.00000012,0.00000001,0.0};

int sgn(double x){
	if(x>MEPS){
		return 1;
	}
	else if(x<-MEPS){
		return -1;
	}
	else {
		return 0;
	}
}

void dist(double *x, double *x0, int n, int p){
	int i,j,k;
	double tmp, tmp1;

	for(i=0; i<n-1; i++){
		for(j=i+1; j<n; j++){
			tmp = 0.0;
			for(k=0; k<p; k++){
				tmp1 = x0[k*n+i]-x0[k*n+j];
				tmp += tmp1*tmp1;
			}
			x[i*n+j] = sqrt(tmp);
		}
	}
	for(i=0; i<n-1; i++){
		for(j=i+1; j<n; j++){
			x[j*n+i] = x[i*n+j];
		}
	}
	for(i=0; i<n; i++){
		x[i*n+i] = 0.0;
	}
}

void rkhs(double *x, double *x0, int n, int p, int type, double sigma2, int d){
	int i,j,k;
	double tmp, tmp1;

	if(type == 1){  // exp(-sigma*||x-y||_2^2)
		for(i=0; i<n-1; i++){
			for(j=i+1; j<n; j++){
				tmp = 0.0;
				for(k=0; k<p; k++){
					tmp1 = x0[k*n+i] - x0[k*n+j];
					tmp += tmp1*tmp1;
				}
				x[i*n+j] = exp(-sigma2*tmp);
			}
		}
		for(i=0; i<n; i++){
			x[i*n+i] = 1.0;
		}
	}
	else if(type == 2){ // exp(-sigma*||x-y||_1)
		for(i=0; i<n-1; i++){
			for(j=i+1; j<n; j++){
				tmp = 0.0;
				for(k=0; k<p; k++){
					tmp	+= fabs(x0[k*n+i] - x0[k*n+j]);
				}
				x[i*n+j] = exp(-sigma2*tmp);
			}
		}
		for(i=0; i<n; i++){
			x[i*n+i] = 1.0;
		}
	}
	else if(type == 3){ // sigma/(sigma+||x-y||_2^2)
		for(i=0; i<n-1; i++){
			for(j=i+1; j<n; j++){
				tmp = 0.0;
				for(k=0; k<p; k++){
					tmp1 	= x0[k*n+i] - x0[k*n+j];
					tmp		+= tmp1*tmp1;
				}
				x[i*n+j] 	= sigma2/(sigma2+tmp);
			}
		}
		for(i=0; i<n; i++){
			x[i*n+i] = 1.0;
		}
	}
	else{
		for(i=0; i<n; i++){
			for(j=i; j<n; j++){ // (sigma+ <x,y>)^d
				tmp = 0.0;
				for(k=0; k<p; k++){
					tmp += x0[k*n+i]*x0[k*n+j];
				}
				x[i*n+j] = pow(sigma2 + tmp, d);
			}
		}
	}
	for(i=0; i<n-1; i++){
		for(j=i+1; j<n; j++){
			x[j*n+i] = x[i*n+j];
		}
	}
}

void rkhs2(double *x, double *x0, double *y0, int n, int p, int type, double sigma2, int d){
	int i,j,k;
	double tmp, tmp1;

	if(type == 1){  // exp(-sigma*||x-y||_2^2)
		for(i=0; i<n-1; i++){
			for(j=i+1; j<n; j++){
				tmp = 0.0;
				for(k=0; k<p; k++){
					tmp1 = x0[k*n+i] - y0[k*n+j];
					tmp += tmp1*tmp1;
				}
				x[i*n+j] = exp(-sigma2*tmp);
			}
		}
		for(i=0; i<n; i++){
			x[i*n+i] = 1.0;
		}
	}
	else if(type == 2){ // exp(-sigma*||x-y||_1)
		for(i=0; i<n-1; i++){
			for(j=i+1; j<n; j++){
				tmp = 0.0;
				for(k=0; k<p; k++){
					tmp	+= fabs(x0[k*n+i] - y0[k*n+j]);
				}
				x[i*n+j] = exp(-sigma2*tmp);
			}
		}
		for(i=0; i<n; i++){
			x[i*n+i] = 1.0;
		}
	}
	else if(type == 3){ // sigma/(sigma+||x-y||_2^2)
		for(i=0; i<n-1; i++){
			for(j=i+1; j<n; j++){
				tmp = 0.0;
				for(k=0; k<p; k++){
					tmp1 	= x0[k*n+i] - y0[k*n+j];
					tmp		+= tmp1*tmp1;
				}
				x[i*n+j] 	= sigma2/(sigma2+tmp);
			}
		}
		for(i=0; i<n; i++){
			x[i*n+i] = 1.0;
		}
	}
	else{
		for(i=0; i<n; i++){
			for(j=i; j<n; j++){ // (sigma+ <x,y>)^d
				tmp = 0.0;
				for(k=0; k<p; k++){
					tmp += x0[k*n+i]*y0[k*n+j];
				}
				x[i*n+j] = pow(sigma2 + tmp, d);
			}
		}
	}
	for(i=0; i<n-1; i++){
		for(j=i+1; j<n; j++){
			x[j*n+i] = x[i*n+j];
		}
	}
}

void Std(double *std, double *x, int n, int p, int flag){
	int i, j;
	double s, s1;
	for(j=0;j<p;j++)
	{
		s = 0; s1 = 0;
		for(i=0;i<n;i++) {
			s += x[j*n+i];
		}
		s = s/n;
		for(i=0;i<n;i++) {
			s1  += x[j*n+i]*x[j*n+i];
		}
		s1 = s1/n - s*s;
		if(flag){
			std[j] = sqrt(s1);
		}
		else{
			std[j] = sqrt(n*s1/(n-1));
		}
	}
}

void Standarize(double *y, double *std, int n, int p, double *x, int flag)
{
	int i, j;
	double s, s1;
	for(j=0;j<p;j++)
	{
		s = 0; s1 = 0;
		for(i=0;i<n;i++) {
			s += x[j*n+i];
		}
		s = s/n;
		for(i=0;i<n;i++) {
			s1  += x[j*n+i]*x[j*n+i];
		}
		s1 = s1/n - s*s;
		if(flag)
			std[j] = sqrt(s1);
		else
			std[j] = sqrt(n*s1/(n-1));
		for(i=0;i<n;i++) {
			y[j*n+i] = (x[j*n+i]-s)/std[j];
		}
	}
}

void sortN0(int *ind, double *x, int n, int dd){
	int i, j, MaxInd, d;
	double tmp, tmp1;

	d = (dd==n?dd-1:dd);
	for(i=0; i<d; i++){
		tmp = x[0]; MaxInd = ind[0];
		for(j=1; j<n-i; j++)
		{
			tmp1 = x[j];
			if(tmp1<tmp){
				x[j-1]      = tmp1;
				x[j]        = tmp;
				ind[j-1]    = ind[j];
				ind[j]      = MaxInd;
			}
			else{
				tmp     = tmp1;
				MaxInd  = ind[j];
			}
		}
	}
}

void sortN(int *ind0, double *x0, int n, int dd){
	int i, j, MaxInd, d, *ind;
	double tmp, *x;
	ind = (int*)malloc(sizeof(int)*n);
	x = (double*)malloc(sizeof(double)*n);
	for(i=0;i<n;i++){
		ind[i] = i;
		x[i] = x0[i];
	}

	d = (dd==n?dd-1:dd);
	for(i=0;i<d;i++)
	{
		tmp = x[0]; MaxInd = ind[0];
		for(j=1;j<n-i;j++)
		{
			if(x[j]<tmp)
			{
				x[j-1] = x[j];
				x[j] = tmp;
				ind[j-1] = ind[j];
				ind[j] = MaxInd;
			}
			else
			{
				tmp = x[j];
				MaxInd = ind[j];
			}
		}
	}
	for(j=0;j<dd;j++) ind0[j] = ind[n-j-1];

	free(ind);
	free(x);
}

void sortN_LS(int *ind, double *x, int n, int dd){
	int i, j, MaxInd, d;
	double tmp, tmp1;

	d = (dd==n?dd-1:dd);
	for(i=0; i<n; i++){
		ind[i] = i;
	}
	for(i=0; i<d; i++){
		tmp = x[0]; MaxInd = ind[0];
		for(j=1; j<n-i; j++)
		{
			tmp1 = x[j];
			if(tmp1>tmp){
				x[j-1]      = tmp1;
				x[j]        = tmp;
				ind[j-1]    = ind[j];
				ind[j]      = MaxInd;
			}
			else{
				tmp     = tmp1;
				MaxInd  = ind[j];
			}
		}
	}
}

void sortN_SL(int *ind, double *x, int n, int dd){
	int i, j, MaxInd, d;
	double tmp, tmp1;

	d = (dd==n?dd-1:dd);
	for(i=0; i<d; i++){
		tmp = x[0]; MaxInd = ind[0];
		for(j=1; j<n-i; j++)
		{
			tmp1 = x[j];
			if(tmp1<tmp){
				x[j-1]      = tmp1;
				x[j]        = tmp;
				ind[j-1]    = ind[j];
				ind[j]      = MaxInd;
			}
			else{
				tmp     = tmp1;
				MaxInd  = ind[j];
			}
		}
	}
}

void SortQ(double *s, int l, int r)
{
    int i, j;
	double x;
    if (l < r)
    {
        i = l;
        j = r;
        x = s[i];
        while (i < j)
        {
            while(i < j && s[j] > x) j--;
			if(i < j) s[i++] = s[j];
            while(i < j && s[i] < x) i++;
			if(i < j) s[j--] = s[i];
        }
        s[i] = x;
        SortQ(s, l, i-1);
        SortQ(s, i+1, r);
    }
}

double minx(double *x, int n){
	int i;
	double mx;
	mx = x[0];
	for(i=1; i<n; i++){
		if(x[i]<mx){
			mx = x[i];
		}
	}
	return mx;
}

double maxx(double *x, int n)
{
	int i;
	double x0=x[0];
	for(i=1;i<n;i++){
		if(x0<x[i]) x0=x[i];
	}
	return(x0);
}

double cumsum(double *x, int n){
	int i;
	double mx;
	mx = x[0];
	for(i=1; i<n; i++){
		mx += x[i];
	}
	return mx;
}

double meanx(double *x, int n){
	int i;
	double mx=0.0;
	for(i=0; i<n; i++){
		mx += x[i];
	}
	return mx/n;
}

void AbyB(double *outVector, double *A, double *v, int n, int p, int q){
	int i,j,k;
	double tmp;
	for (i=0;i<n;i++){
		for(k=0;k<q;k++){
			tmp = 0;
			for(j=0;j<p;j++){
				tmp += A[j*n + i]*v[k*p + j];
			}
			outVector[k*n+i] = tmp;
		}
	}
}

void tAbyB(double *outMatrix, const double *A, const double *B, int n, int p, int q){
    int i,j,k;
    double temp;
	for (i = 0; i<p; i++){
		for (k = 0; k<q; k++){
			temp = 0.0;
			for (j = 0; j < n; j++)
				temp += A[i + j*p] * B[k + j*q];
			outMatrix[i*q + k] = temp;
		}
	}
}

int LowTriangularInv(double *B, int n, double *A){
	// Input:
	// A is a lower triangular matrix
	//
	// Output:
	// B = inv(A)
	//
	int i,j,k;
	const double EPS=DBL_EPSILON;
	for(i=0;i<n;i++)
		if(fabs(A[i*n+i])<EPS)	return(0);
	for(i=0;i<n;i++)	B[i*n+i] = 1;
	for(j=1;j<n;j++)
		for(i=0;i<j;i++)	B[j*n+i] = 0;

	for(i=n-1;i>=0;i--)//rows
	{
		if(fabs(A[i*n+i]-1)>EPS)
			for(j=i;j<n;j++)
				B[j*n+i] = B[j*n+i]/A[i*n+i];
		if(i>0)
		{
			for(j=i;j<n;j++)// columns
				for(k=0;k<i;k++)// rows
					B[j*n+k] = B[j*n+k] - A[i*n+k]*B[j*n+i];
		}
	}
	return(1);
}

int UpTriangularInv(double *B, int n, double *A)
{ 
	int i,j,k;
	for(i=0;i<n;i++)
		if(fabs(A[i*n+i])<1e-4) return(0);
	for(i=0;i<n;i++) B[i*n+i] = 1;
	for(j=1;j<n;j++)for(i=0;i<j;i++)B[j*n+i] = 0;

	for(i=n-1;i>=0;i--)//rows
	{
		if(A[i*n+i]!=1)
			for(j=i;j<n;j++)
				B[j*n+i] = B[j*n+i]/A[i*n+i];
		if(i>0)
		{
			for(j=i;j<n;j++)// columns
				for(k=0;k<i;k++)// rows
					B[j*n+k] = B[j*n+k] - A[i*n+k]*B[j*n+i];
		}
	}
	return(1);
}

void QRDecompN(double *E, double *R, double *x, int n, int p){
	// Input:
	// X is a p*n matrix
	//
	// Output:
	// R is a p*p lower triangular matrix
	// E is a p*n matrix satisfying E*t(E) = I_p
	//
	double *Z, *znorm;
	double  tmp, tmp1;
	int i,j, k;

	Z = (double*)malloc(sizeof(double)*n*p);
	znorm = (double*)malloc(sizeof(double)*p);

	// calculate the first column
	tmp = 0;
	for(i=0;i<n;i++){
		Z[i] = x[i];
		tmp += Z[i]*Z[i];
	}
	znorm[0] = sqrt(tmp);
	tmp = 0;
	for(i=0;i<n;i++){
		E[i] = x[i]/znorm[0];
		tmp += E[i]*x[i];
	}
	R[0] = tmp;

	//iteration from j=1...p
	for(j=1;j<p;j++){
		for(k=0;k<j;k++){
			tmp=0;	for(i=0;i<n;i++) tmp += E[k*n+i]*x[j*n+i];
			R[j*p+k] = tmp;
		}
		tmp1 = 0;
		for(i=0;i<n;i++){
			tmp = 0; for(k=0;k<j;k++) tmp += R[j*p+k]*E[k*n+i];
			Z[j*n+i] = x[j*n+i] - tmp;
			tmp1 += pow(Z[j*n+i],2);
		}
		znorm[j] = sqrt(tmp1);
		tmp1 = 0;
		for(i=0;i<n;i++) E[j*n+i] = Z[j*n+i]/znorm[j];
		for(i=0;i<n;i++) tmp1 += E[j*n+i]*x[j*n+i];
		R[j*p+j] = tmp1;
	}
	free(Z); free(znorm);
}

double SampleQuantile1(double *z, int n, double q){
	int i, ind;
	double qr;
	double *zs=(double*)malloc(sizeof(double)*n);
	for(i=0;i<n;i++) zs[i] = z[i];
	SortQ(zs, 0, n-1);

	ind = floor(q*n);
	if (ind!=n*q){
		qr = zs[ind];
	}
	else{
		qr = (zs[ind-1] + zs[ind])/2;
	}

	free(zs);
	return qr;
}

void SampleQuantile(double *qr, int m, double *z, int n, double *q)
{
	double *zs=(double*)malloc(sizeof(double)*n);
	int i, ind;
	for(i=0;i<n;i++) zs[i] = z[i];
	SortQ(zs, 0, n-1);
	for(i=0;i<m;i++)
	{
		ind = floor(q[i]*n);
		if (ind!=n*q[i])
			qr[i] = zs[ind];
		else
			qr[i] = (zs[ind-1] + zs[ind])/2;
	}
	free(zs);
}

int MatrixInvSymmetric(double *a,int n){
	int i,j,k,m;
    double w,g,*b;
    b = (double*)malloc(n*sizeof(double));

    for (k=0; k<=n-1; k++){
        w=a[0];
        if (fabs(w)+1.0==1.0){
            free(b); return(-2);
        }
        m=n-k-1;
        for (i=1; i<=n-1; i++){
            g=a[i*n]; b[i]=g/w;
            if (i<=m) b[i]=-b[i];
            for (j=1; j<=i; j++)
                a[(i-1)*n+j-1]=a[i*n+j]+g*b[j];
        }
        a[n*n-1]=1.0/w;
        for (i=1; i<=n-1; i++)
            a[(n-1)*n+i-1]=b[i];
    }
    for (i=0; i<=n-2; i++)
        for (j=i+1; j<=n-1; j++)
            a[i*n+j]=a[j*n+i];
    free(b);
    return(2);
}

int LinearSolveSym(double *B, int n, double *a)
{ 
	int i,j,k,m;
	double w,g,*b;
	b = (double*)malloc(sizeof(double)*n);

	for (k=0; k<=n-1; k++)
		{ w=a[0];
		if (fabs(w)+1.0==1.0)
			{ free(b); return(-2);}
		m=n-k-1;
		for (i=1; i<=n-1; i++)
			{ g=a[i*n]; b[i]=g/w;
			if (i<=m) b[i]=-b[i];
			for (j=1; j<=i; j++)
				a[(i-1)*n+j-1]=a[i*n+j]+g*b[j];
			}
		a[n*n-1]=1.0/w;
		for (i=1; i<=n-1; i++)
			a[(n-1)*n+i-1]=b[i];
		}
	for (i=0; i<n; i++){
		w=0; 
		for(j=i;j<n;j++) w += a[j*n+i]*B[j];
		for(j=0;j<i;j++) w += a[i*n+j]*B[j];
		b[i]=w;
	}
	for(i=0;i<n;i++)B[i]=b[i];
	free(b);
	return(2);
}

void freematrix(double **x, int n)
{
	for (int i = 0; i < n; ++i) {
		free(x[i]);
	}
	free(x);
}

double** requirematrix(int n, int p, int zero)
{
	double **x;
	x = (double**)malloc(sizeof(double*)*n);
	if (zero == 1) {
		for (int i = 0; i < n; ++i) {
			x[i] = (double*)calloc(p, sizeof(double));
			if (x[i] == NULL){
				fprintf(stderr, "Unable to allocate enough memory for array!\n");
			}
		}
	} 
	else {
		for (int i = 0; i < n; ++i) {
			x[i] = (double*)malloc(p*sizeof(double));
			if (x[i] == NULL){
				fprintf(stderr, "Unable to allocate enough memory for array!\n");
			}
		}
	}
	return x;
}

double pythag(double a, double b)
{
    double absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb){
		return absa*sqrt(1.0+(absb/absa)*(absb/absa));
	}
    else{
		return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
	}
}

double SIGN(const double a, const double b){
	return (b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
}

void tred2(double **z, double *d, double *e, int n, int yesvecs)
{
    int l,k,j,i;
    double scale,hh,h,g,f;
    for (i=n-1;i>0;i--) {
        l=i-1;
        h=scale=0.0;
        if (l > 0) {
            for (k=0;k<i;k++)
                scale += fabs(z[i][k]);
            if (scale == 0.0)
                e[i]=z[i][l];
            else {
                for (k=0;k<i;k++) {
                    z[i][k] /= scale;
                    h += z[i][k]*z[i][k];
                }
                f=z[i][l];
                g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
                e[i]=scale*g;
                h -= f*g;
                z[i][l]=f-g;
                f=0.0;
                for (j=0;j<i;j++) {
                    if (yesvecs)
                        z[j][i]=z[i][j]/h;
                    g=0.0;
                    for (k=0;k<j+1;k++)
                        g += z[j][k]*z[i][k];
                    for (k=j+1;k<i;k++)
                        g += z[k][j]*z[i][k];
                    e[j]=g/h;
                    f += e[j]*z[i][j];
                }
                hh=f/(h+h);
                for (j=0;j<i;j++) {
                    f=z[i][j];
                    e[j]=g=e[j]-hh*f;
                    for (k=0;k<j+1;k++)
                        z[j][k] -= (f*e[k]+g*z[i][k]);
                }
            }
        } else
            e[i]=z[i][l];
        d[i]=h;
    }
    if (yesvecs) d[0]=0.0;
    e[0]=0.0;
    for (i=0;i<n;i++) {
        if (yesvecs) {
            if (d[i] != 0.0) {
                for (j=0;j<i;j++) {
                    g=0.0;
                    for (k=0;k<i;k++)
                        g += z[i][k]*z[k][j];
                    for (k=0;k<i;k++)
                        z[k][j] -= g*z[k][i];
                }
            }
            d[i]=z[i][i];
            z[i][i]=1.0;
            for (j=0;j<i;j++) z[j][i]=z[i][j]=0.0;
        } else {
            d[i]=z[i][i];
        }
    }
}

void tqli(double **z, double *d, double *e, int n, int yesvecs)
{
    int m,l,iter,i,k;
    double s,r,p,g,f,dd,c,b;
    const double EPS=DBL_EPSILON;
    for (i=1;i<n;i++) e[i-1]=e[i];
    e[n-1]=0.0;
    for (l=0;l<n;l++) {
        iter=0;
        do {
            for (m=l;m<n-1;m++) {
                dd=fabs(d[m])+fabs(d[m+1]);
                if (fabs(e[m]) <= EPS*dd) break;
            }
            if (m != l) {
                if (iter++ == 50) {
                    // printf("Too many iterations in eigencmp.\n");
					break;
                    // exit(1);
                }
                g=(d[l+1]-d[l])/(2.0*e[l]);
                r=pythag(g,1.0);
                g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
                s=c=1.0;
                p=0.0;
                for (i=m-1;i>=l;i--) {
                    f=s*e[i];
                    b=c*e[i];
                    e[i+1]=(r=pythag(f,g));
                    if (r == 0.0) {
                        d[i+1] -= p;
                        e[m]=0.0;
                        break;
                    }
                    s=f/r;
                    c=g/r;
                    g=d[i+1]-p;
                    r=(d[i]-g)*s+2.0*c*b;
                    d[i+1]=g+(p=s*r);
                    g=c*r-b;
                    if (yesvecs) {
                        for (k=0;k<n;k++) {
                            f=z[k][i+1];
                            z[k][i+1]=s*z[k][i]+c*f;
                            z[k][i]=c*z[k][i]-s*f;
                        }
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l]=g;
                e[m]=0.0;
            }
        } while (m != l);
    }
}

void eigencmp(double **z, double *d, int n, int yesvecs){
    double *e = (double*)calloc(n, sizeof(double));
    tred2(z, d, e, n, yesvecs);
    tqli(z, d, e, n, yesvecs);
    free(e);
}

void _ginv_(double *Sn_inv, double *Sn, int p, double tol){
	int i,j,k,*ind;
	double tmp, **U, *singular;
	ind 	= (int*)malloc(sizeof(int)*p);
	U 		= requirematrix(p, p, 0);
	singular = (double*)malloc(sizeof(double)*p);

	for (j = 0; j < p; ++j) {
		U[j][j] = Sn[j*p+j];
		for (i = j+1; i < p; ++i){
			U[j][i] = U[i][j] = Sn[j*p+i];
		}
	}

	eigencmp(U, singular, p, 1);
	sortN_LS(ind, singular, p, p);

	for (j = 0; j < p; j++) {
		singular[j] = (singular[j]>tol)?(1.0/singular[j]):0.0;
	}
	for (j = 0; j < p; j++) {
		for (i = j; i < p; i++) {
			tmp = 0.0;
			for (k = 0; k < p; k++){
				tmp += U[i][ind[k]]*U[j][ind[k]]*singular[k];
			}
			Sn_inv[i*p+j] = tmp;
			if (i != j){
				Sn_inv[j*p+i] = Sn_inv[i*p+j];
			}
		}
	}

	freematrix(U, p);
	free(singular);
	free(ind);
}

void _eigdecom_(double *Sn_inv, double *eigvects, double *eigvals, double *Sn, int p, double tol){
	int i,j,k,*ind;
	double tmp, **U, *singular;
	ind 	= (int*)malloc(sizeof(int)*p);
	U 		= requirematrix(p, p, 0);
	singular = (double*)malloc(sizeof(double)*p);

	for (j = 0; j < p; ++j) {
		U[j][j] = Sn[j*p+j];
		for (i = j+1; i < p; ++i){
			U[j][i] = U[i][j] = Sn[j*p+i];
		}
	}

	eigencmp(U, singular, p, 1);
	sortN_LS(ind, singular, p, p);

	for (j = 0; j < p; j++) {
		eigvals[j] 	= singular[j];
		singular[j] = (singular[j]>tol)?(1.0/singular[j]):0.0;
	}
	for (j = 0; j < p; j++) {
		for (i = j; i < p; i++) {
			tmp = 0.0;
			for (k = 0; k < p; k++){
				tmp += U[i][ind[k]]*U[j][ind[k]]*singular[k];
			}
			Sn_inv[i*p+j] 	= tmp;
			if (i != j){
				Sn_inv[j*p+i] 	= Sn_inv[i*p+j];
			}
		}
	}

	for (j = 0; j < p; j++) {
		for (i = 0; i < p; i++) {
			eigvects[j*p+i] = U[i][ind[j]];
		}
	}

	freematrix(U, p);
	free(singular);
	free(ind);
}

double calculate_bic(double loss, int *param, double gamma, int model){
	int n       = param[0];
	int p       = param[1];
	int number  = param[2];

	if(model == 1){
		return 2.0*n*(loss) + number*log(n) + 2.0*gamma*lchoose(p, number);
	}
	else if(model == 2){
		return 2.0*n*(loss) + number*log(n) + 2.0*gamma*lchoose(p, number);
	}
	else if(model == 3){
		return 2.0*n*(loss) + number*log(n) + 2.0*gamma*lchoose(p, number);
	}
	else{
		return 2.0*n*(loss) + number*log(n) + 2.0*gamma*lchoose(p, number);
	}
}

void usi(int *s_i, int *s_j, int *tt_b, int *tt_a, int *act, int ns, int iter)
{
	int i;
	s_i += *tt_a;
	s_j += *tt_a;
	for (i = 0; i < ns; ++i)
	{
		s_i[i] = act[i];
		s_j[i] = iter;
	}
	*tt_b   = *tt_a;
	*tt_a  += ns;
}

void preprocess(double *x, int n, int p, double *sdx){
	int i, j;
	double cor, var, *x_j, tmp1;

	for(j=0; j<p; j++){
		cor = 0.0;
		var = 0.0;
		x_j = x + j*n;
		for(i=0; i<n; i++){
			tmp1 = x_j[i];
			cor += tmp1;
			var += tmp1*tmp1;
		}
		tmp1    = cor/n;
		sdx[j]  = sqrt((var-n*tmp1*tmp1)/(n-1));
	}
}

void BsplineBasis(int order, double *knots, int nknots, double *x, int n, double *B, int j)
{
	int i;
	double *b, dd, dn;
	b = (double*)malloc(sizeof(double)*n);


	for(i=0;i<n;i++)b[i]=0;
	if(order>1){
		BsplineBasis(order-1,knots,nknots,x,n,B,j);
		dd = knots[j+order-1] - knots[j];
		if(dd!=0){
			for(i=0;i<n;i++){
				dn = x[i]- knots[j];
				b[i] += B[i]*dn/dd;
			}
		}
		BsplineBasis(order-1,knots,nknots,x,n,B,j+1);
		dd = knots[j+order] - knots[j+1];
		if(dd!=0){
			for(i=0;i<n;i++){
				dn = knots[j+order] - x[i];
				b[i] += B[i]*dn/dd;
			}
		}
	}
	else if(knots[j+1]<knots[nknots-1]){
			for(i=0;i<n;i++){
				if(knots[j]<=x[i]&&knots[j+1]>x[i]) {
					b[i]=1;
				}
			}
	}
	else{
		for(i=0;i<n;i++){
			if(knots[j]<=x[i]) b[i]=1;
		}
	}
	for(i=0;i<n;i++){
		B[i] = b[i];
	}
	free(b);
}

void BsplineMatrix(int order, double *knots, int nknots, double *x, int n, double *B)
{
	int i,j;
	double *b;
	b = (double*)malloc(sizeof(double)*n);

	for(j=0;j<nknots-order;j++){
		BsplineBasis(order, knots, nknots, x, n, b,j);
		for(i=0;i<n;i++) {
			B[j*n+i] = b[i];
		}
	}
	free(b);
}

void BsplineDerivative(int order, double *knots, int nknots, double *x, int n, double *B, int d)
{
	int i,j, nk = nknots - order;
	double *b, *new_knots, tmp1, tmp2, tmp3;

	//if(d>order-1) error("Knot vector values should be nondecreasing.");
	if(d==0){
		BsplineMatrix(order, knots, nknots, x, n, B);
	}
	else if(d==1){
		b = (double*)malloc(sizeof(double)*n*(nk-1));
		new_knots = (double*)malloc(sizeof(double)*(nknots-2));

		for(j=0;j<nknots-2;j++) {
			new_knots[j] = knots[j+1];
		}
		BsplineMatrix(order-1, new_knots, nknots-2, x, n, b);

		tmp1 = -(order-1)/(knots[order]-knots[1]); 	
		for(i=0;i<n;i++) {
			B[i] = b[i]*tmp1;
		}
		for(j=1;j<nk-1;j++){
			tmp1 =  (order-1)/(knots[order+j-1]-knots[j]); 
			tmp2 = (order-1)/(knots[order+j]-knots[j+1]);
			for(i=0;i<n;i++){
				B[j*n+i] = b[(j-1)*n+i]*tmp1 - b[j*n+i]*tmp2;
			}
		}
		tmp2 = (order-1)/(knots[nknots-2]-knots[nk-1]);  
		for(i=0;i<n;i++) {
			B[n*(nk-1)+i] =b[(nk-2)*n+i]*tmp2;
		}
		free(b); 
		free(new_knots);
	}
	else{
		b = (double*)malloc(sizeof(double)*n*(nk-2));
		new_knots = (double*)malloc(sizeof(double)*(nknots-4));

		for(j=0;j<nknots-4;j++) {
			new_knots[j] = knots[j+2];
		}
		BsplineMatrix(order-2, new_knots, nknots-4, x, n, b);

		tmp1 = (order-1)*(order-2)/(knots[order]-knots[2])/(knots[order]-knots[1]); 
		for(i=0;i<n;i++) {
			B[i] = b[i]*tmp1;
		}

		tmp1 = -(order-1)*(order-2)/(knots[order]-knots[2])*( 1.0/(knots[order]-knots[1]) + 1.0/(knots[order+1]-knots[2]) ); 
		tmp2 = (order-1)*(order-2)/(knots[order+1]-knots[3])/(knots[order+1]-knots[2]);
		for(i=0;i<n;i++) {
			B[n+i] = b[i]*tmp1 + b[n+i]*tmp2;
		}

		for(j=2;j<nk-2;j++){
			tmp1 =  (order-1)*(order-2)/(knots[order+j-2]-knots[j])/(knots[order+j-1]-knots[j]);
			tmp2 = -(order-1)*(order-2)/(knots[order+j-1]-knots[j+1])*( 1.0/(knots[order+j-1]-knots[j]) + 1.0/(knots[order+j]-knots[j+1]) );
			tmp3 = (order-1)*(order-2)/(knots[order+j]-knots[j+2])/(knots[order+j]-knots[j+1]);
			for(i=0;i<n;i++)  {
				B[j*n+i] = b[(j-2)*n+i]*tmp1 + b[(j-1)*n+i]*tmp2 + b[j*n+i]*tmp3;
			}
		}

		tmp2 = (order-1)*(order-2)/(knots[nknots-4]-knots[nk-2])/(knots[nknots-3]-knots[nk-2]);  
		tmp3 = -(order-1)*(order-2)/(knots[nknots-3]-knots[nk-1])*( 1.0/(knots[nknots-3]-knots[nk-2]) + 1.0/(knots[nknots-2]-knots[nk-1]) );
		for(i=0;i<n;i++) {
			B[n*(nk-2)+i] = b[(nk-4)*n+i]*tmp2 + b[(nk-3)*n+i]*tmp3;
		}


		tmp3 = (order-1)*(order-2)/(knots[nknots-3]-knots[nk-1])/(knots[nknots-2]-knots[nk-1]);
		for(i=0;i<n;i++) {
			B[n*(nk-1)+i] =b[(nk-3)*n+i]*tmp3;
		}
		free(b); 
		free(new_knots);
	}
}

void ParamEst_Logist(double *hatxi, double *y, double *xstar, int n, int p, int maxstep, double eps)
{
	int i,j,k, step=0, ninner;
	double *d1, *d2, *nxi, *p0, tmp1, expx, ologelr, nlogelr, delta;

	d1 = (double*)malloc(sizeof(double)*p);
	d2 = (double*)malloc(sizeof(double)*p*p);
	nxi = (double*)malloc(sizeof(double)*p);
	p0 = (double*)malloc(sizeof(double)*n);

	while(step<maxstep){
		step++;
		ologelr = 0;
		for(i=0;i<n;i++){
			tmp1 = 0; 
			for(j=0;j<p;j++) {
				tmp1 += xstar[j*n+i]*hatxi[j];
			}
			expx = exp(tmp1);
			ologelr += y[i]*tmp1 - log(1+expx);
			p0[i] = expx/(1+expx);			
		}
		for(j=0;j<p;j++){
			tmp1 = 0;	
			for(i=0;i<n;i++) {
				tmp1 += xstar[j*n+i]*(y[i]-p0[i]);
			}
			d1[j] = tmp1;
			for(k=0;k<p;k++){
				tmp1 = 0;	
				for(i=0;i<n;i++) {
					tmp1 += xstar[j*n+i]*p0[i]*(1-p0[i])*xstar[k*n+i];
				}
				d2[j*p+k] = - tmp1;
			}
		}
		LinearSolveSym(d1, p, d2);

		/* adjusted step for newton algorithm */
		ninner = 16;
		for(k=0;k<16;k++){
			for(j=0;j<p;j++) {
				nxi[j] = hatxi[j] - gradient_wts[k]*d1[j];
			}
			nlogelr = 0;
			for(i=0;i<n;i++){
				expx = 0;	
				for(j=0;j<p;j++) {
					expx += xstar[j*n+i]*nxi[j];
				}
				nlogelr += y[i]*expx - log(1+exp(expx));
			}
			if(nlogelr>ologelr){
				for(j=0;j<p;j++) {
					hatxi[j] = nxi[j];
				}
				ninner = k; break;
			}
		}
		if(ninner==16) break;

		delta = 0.0;
		for(j=0;j<p;j++){
			tmp1 = fabs(gradient_wts[ninner]*d1[j]);
			if(delta<tmp1) {
				delta = tmp1;
			}
		}
		if(delta<eps) break;
	}//end while
	free(p0); 
	free(d1); 
	free(d2); 
	free(nxi);
}

void EstKernel_inner(double *g, double *w, double *y, double *w0, int n, int nw0, double h, double *G)
{
	int i, j;
	double tmp1, tmp2, kern0, k20, k10, k01, k11, k00, G0;

	for(j=0;j<nw0;j++){
		k20 = 0; 
		k10 = 0; 
		k01 = 0; 
		k11 = 0; 
		k00 = 0;
		for(i=0;i<n;i++){
			tmp1 =w[i]-w0[j];
			tmp2 = tmp1*tmp1/h/h;
			kern0 = 0.75*(1-tmp2)*(tmp2<1?1:0)/h;
			G0 = G[i];
			k20 += kern0*tmp1*tmp1*G0*G0;
			k10 += kern0*tmp1*G0*G0;
			k01 += y[i]*kern0*G0;
			k11 += y[i]*kern0*tmp1*G0;
			k00 += kern0*G0*G0;
		}
		tmp1 = k00*k20-k10*k10;
		if(fabs(tmp1)<1.0/n) {
			g[j]=0;
		}
		else {
			g[j] = (k20*k01 - k10*k11)/tmp1;
		}
	}
}

void Est_Kernel(double *hatg, double *y, double *G, double *haty, double *w, int n, int p, int q, int nf, double *h, double *g, double *w0, int nw0)
{
	int i,j, k;
	double *yy, *hatg1, *G1;

	hatg1 	= (double*)malloc(sizeof(double)*nw0);
	G1 		= (double*)malloc(sizeof(double)*n);
	yy 		= (double*)malloc(sizeof(double)*n);

	for(k=0;k<nf;k++){
		for(i=0;i<n;i++){
			G1[i] = G[n*k+i]; 
		}
		for(j=0;j<p;j++){
			for(i=0;i<n;i++) {
				yy[i] = y[i] - haty[i] + g[(k*p+j)*n+i];
			}
			EstKernel_inner(hatg1, w, yy, w0, n, nw0, h[k], G1);
			for(i=0;i<nw0;i++) {
				hatg[(k*p+j)*nw0+i] = hatg1[i];
			}
		}
	}
	free(yy); 
	free(hatg1); 
	free(G1);
}

void OLS(double *beta, double *x, double *y, int n, int p)
{
	int i,j;
	double *Q, *R, *invR, *b, tmp;
	Q 	= (double*)malloc(sizeof(double)*n*p);
	R 	= (double*)malloc(sizeof(double)*p*p);
	b 	= (double*)malloc(sizeof(double)*p);
	invR = (double*)malloc(sizeof(double)*p*p);
	
	QRDecompN(Q, R, x, n, p);
	UpTriangularInv(invR, p, R);

	for(j=0;j<p;j++){
		tmp = 0; 
		for(i=0;i<n;i++) {
			tmp += Q[j*n+i]*y[i];
		}
		b[j] = tmp;
	}
	for(j=0;j<p;j++){
		tmp = 0; 
		for(i=j;i<p;i++) {
			tmp += invR[i*p+j]*b[i];
		}
		beta[j] = tmp;
	}	
	free(Q);
	free(R); 
	free(invR); 
	free(b);
}

void Kernelh(double *x, double *kern, int n, double h0, int type){
	// type = 1: Gaussian kernel
	// type = 2; Epanecknikov kernel
	// type = 3: uniform kernel
	int i;
	double tmp, h;
	Std( &h, x, n, 1, 0);
	h /= h0;
	for(i=0; i<n; i++){
		tmp = x[i]*h;
		tmp *= tmp;
		if(type==1){
			kern[i] = exp(-0.5*tmp)*MPI2*h;
		}
		else if(type==2){
			kern[i] = 0.75*(1-tmp)*(tmp<1?1:0)*h;
		}
		else{
			kern[i] = (tmp>0.5?0:1)*(tmp<-0.5?0:1)*h;
		}
	}
}

void EstLinear(double *beta, double *x, double *z, double *y, double *weight, int n, int p, double h, int type){
	int i,j,k;
	double tmp, *hess, *xy;
	xy 		= (double*)malloc(sizeof(double)*p);
	hess	= (double*)malloc(sizeof(double)*p*p);

	Kernelh(z, weight, n, h, type);
	for(j=0;j<p;j++) xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p;j++){
			xy[j] 	+= x[j*n+i]*y[i]*weight[i];
		}
	}
	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i]*weight[i];
			}
			hess[j*p+k] = tmp;
		}
	}

	MatrixInvSymmetric(hess,p);
	AbyB(beta, hess, xy, p, p, 1);

	free(hess);
	free(xy);
}

void EstLogistic(double *beta, double *x, double *z, double *y, double *weight, int n, int p, double h, int maxstep, double eps, int type){
	int i,j,k, step=0;
	double *beta0, *hess, *qy, *dpsi;
	double tmp, bnorm, yk, expx, wk;

	beta0 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	hess 	= (double*)malloc(sizeof(double)*p*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	Kernelh(z, weight, n, h, type);
	for(j=0;j<p;j++)	beta0[j] 	= 0.0;

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += x[j*n+i]*beta0[j];
			}
			expx = exp(tmp);
			wk 	= expx/(1.0+expx);
			yk 	= wk - y[i];
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk*weight[i];
				dpsi[j*n+i] = x[j*n+i]*wk*(1.0-wk)*weight[i];
			}
		}
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		MatrixInvSymmetric(hess,p);
    	AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp 	=  beta[j];
			bnorm 	+= tmp*tmp;
			beta[j] = beta0[j] - tmp;
		}
		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)
				beta0[j] = beta[j];
		}
	}

	free(beta0);
	free(hess);
	free(dpsi);
	free(qy);
}

void EstPoiss(double *beta0, double *x, double *z, double *y, double *weight, int n, int p, double h, int maxstep, double eps, int type){
	int i,j,k, step=0;
	double *beta, *qy, *dpsi, *hess;
	double tmp, bnorm, yk, wk;

	beta 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);
	hess 	= (double*)malloc(sizeof(double)*p*p);

	Kernelh(z, weight, n, h, type);

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += x[j*n+i]*beta0[j];
			}
			wk 	= exp(tmp);
			yk 	= wk - y[i];
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk*weight[i];
				dpsi[j*n+i] = x[j*n+i]*wk*weight[i];
			}
		}
		for(j=0; j < p; j++){
			for(k=j; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		for(j=1; j < p; j++){
			for(k=0; k < j; k++){
				hess[j*p+k] = hess[k*p+j];
			}
		}
		MatrixInvSymmetric(hess,p);
    	AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp 	=  beta[j];
			bnorm 	+= tmp*tmp;
			beta[j] = beta0[j] - tmp;
		}
		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)
				beta0[j] = beta[j];
		}
	}

	free(beta);
	free(qy);
	free(dpsi);
	free(hess);
}

void EstLinearW(double *beta, double *residual, double *x, double *z, double *y, double *weight, int n, int p, double h, int type){
	int i,j,k;
	double tmp, *hess, *xy;
	xy 		= (double*)malloc(sizeof(double)*p);
	hess	= (double*)malloc(sizeof(double)*p*p);

	Kernelh(z, weight, n, h, type);
	for(j=0;j<p;j++) xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p;j++){
			xy[j] 	+= x[j*n+i]*y[i]*weight[i];
		}
	}
	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i]*weight[i];
			}
			hess[j*p+k] = tmp;
		}
	}

	MatrixInvSymmetric(hess,p);
	AbyB(beta, hess, xy, p, p, 1);

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0; j < p; j++){
			tmp += x[j*n+i]*beta[j];
		}
		residual[i] = y[i] - tmp;
	}

	free(hess);
	free(xy);
}

void EstLogisticW(double *beta, double *residual, double *x, double *z, double *y, double *weight, int n, int p, double h, int maxstep, double eps, int type){
	int i,j,k, step=0;
	double *beta0, *hess, *qy, *dpsi;
	double tmp, bnorm, yk, expx, wk;

	beta0 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	hess 	= (double*)malloc(sizeof(double)*p*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	Kernelh(z, weight, n, h, type);
	for(j=0;j<p;j++)	beta0[j] 	= 0.0;

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += x[j*n+i]*beta0[j];
			}
			expx = exp(tmp);
			wk 	= expx/(1.0+expx);
			yk 	= wk - y[i];
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk*weight[i];
				dpsi[j*n+i] = x[j*n+i]*wk*(1.0-wk)*weight[i];
			}
		}
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		MatrixInvSymmetric(hess,p);
    	AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp 	=  beta[j];
			bnorm 	+= tmp*tmp;
			beta[j] = beta0[j] - tmp;
		}
		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)
				beta0[j] = beta[j];
		}
	}

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++){
			tmp += x[j*n+i]*beta[j];
		}
		expx = exp(tmp);
		wk 	= expx/(1.0+expx);
		yk 	= y[i] - wk;
		residual[i] = yk;
	}

	free(beta0);
	free(hess);
	free(dpsi);
	free(qy);
}

void EstPoissW(double *beta0, double *residual, double *x, double *z, double *y, double *weight, int n, int p, double h, int maxstep, double eps, int type){
	int i,j,k, step=0;
	double *beta, *qy, *dpsi, *hess;
	double tmp, bnorm, yk, wk;

	beta 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);
	hess 	= (double*)malloc(sizeof(double)*p*p);

	Kernelh(z, weight, n, h, type);
	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += x[j*n+i]*beta0[j];
			}
			wk 	= exp(tmp);
			yk 	= wk - y[i];
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk*weight[i];
				dpsi[j*n+i] = x[j*n+i]*wk*weight[i];
			}
		}
		for(j=0; j < p; j++){
			for(k=j; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		for(j=1; j < p; j++){
			for(k=0; k < j; k++){
				hess[j*p+k] = hess[k*p+j];
			}
		}
		MatrixInvSymmetric(hess,p);
    	AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp 	=  beta[j];
			bnorm 	+= tmp*tmp;
			beta[j] = beta0[j] - tmp;
		}
		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)
				beta0[j] = beta[j];
		}
	}

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++){
			tmp += x[j*n+i]*beta[j];
		}
		residual[i] = y[i] - exp(tmp);
	}

	free(beta);
	free(qy);
	free(dpsi);
	free(hess);
}

void _BandMatrix(double *sighalf, double *rho, int p, int T){
	int i,j,n=p+T-1;
	for(i=0;i<n*p;i++){
		sighalf[i] = 0.0;
	}
	for(j=0;j<p;j++){
		for(i=0;i<T;i++){
			sighalf[j*n+j+i] = rho[i];
		}
	}
}

void _Toeplitz(double *sighalf, double *rho, int p){
	int i,j;

	for(j = 0; j<p; j++){
		for(i = j; i<p; i++){
			sighalf[j*p+i] = rho[i-j];
		}
	}

	for(j = 0; j< p-1; j++){
		for(i = j+1; i<p; i++){
			sighalf[i*p+j] = sighalf[j*p+i];
		}
	}
}
