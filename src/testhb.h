#ifndef RBS_H_INCLUDED
#define RBS_H_INCLUDED

#define EPS_PARAM 1E-6
#define EPS_KERNEL 1E-6
#define MEPS 1e-10
#define FPMIN 1.0e-30
#define MPI 3.1415926
#define MPI1 0.1591549   //1.0/(2*MPI)
#define MPI2 0.3989423   //1.0/sqrt(2*MPI)
#define SQRT2 0.7071068  //1.0/sqrt(2.0)
#define LOG2PI 1.837877  //log(2pi)
#define MIN(a, b) (((a)<(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))
#define IDEX(a, b) (((a)<(b))?1:0)

int sgn(double x);

void dist(double *x, double *x0, int n, int p);

void rkhs(double *x, double *x0, int n, int p, int type, double sigma2, int d);

void rkhs2(double *x, double *x0, double *y0, int n, int p, int type, double sigma2, int d);

void Std(double *std, double *x, int n, int p, int flag);

void Standarize(double *y, double *std, int n, int p, double *x, int flag);

void sortN(int *ind0, double *x, int n, int dd);

void sortN0(int *ind, double *x, int n, int dd);

void sortN_LS(int *ind, double *x, int n, int dd);

void sortN_SL(int *ind, double *x, int n, int dd);

void SortQ(double *s, int l, int r);

double minx(double *x, int n);

double maxx(double *x, int n);

double cumsum(double *x, int n);

double meanx(double *x, int n);

void AbyB(double *outVector, double *A, double *v, int n, int p, int q);

void tAbyB(double *outMatrix, const double *A, const double *B, int n, int p, int q);

int LowTriangularInv(double *B, int n, double *A);

int UpTriangularInv(double *B, int n, double *A);

void QRDecompN(double *E, double *R, double *x, int n, int p);

void freematrix(double **x, int n);

double** requirematrix(int n, int p, int zero);

double pythag(double a, double b);

double SIGN(const double a, const double b);

void tred2(double **z, double *d, double *e, int n, int yesvecs);

void tqli(double **z, double *d, double *e, int n, int yesvecs);

void eigencmp(double **z, double *d, int n, int yesvecs);

void matrixprod(double **x, double **y, double **xy, int n, int p, int q, int transpose);

void _ginv_(double *Sn_inv, double *Sn, int p, double tol);

void _eigdecom_(double *Sn_inv, double *eigvects, double *eigvals, double *Sn, int p, double tol);

double SampleQuantile1(double *z, int n, double q);

void SampleQuantile(double *qr, int m, double *z, int n, double *q);

int MatrixInvSymmetric(double *a,int n);

int LinearSolveSym(double *B, int n, double *a);

double calculate_bic(double loss, int *param, double gamma, int model);

void usi(int *s_i, int *s_j, int *tt_b, int *tt_a, int *act, int ns, int iter);

void preprocess(double *x, int n, int p, double *sdx);

void BsplineBasis(int order, double *knots, int nknots, double *x, int n, double *B, int j);

void BsplineMatrix(int order, double *knots, int nknots, double *x, int n, double *B);

void BsplineDerivative(int order, double *knots, int nknots, double *x, int n, double *B, int d);

void ParamEst_Logist(double *hatxi, double *y, double *xstar, int n, int p, int maxstep, double eps);

void EstKernel_inner(double *g, double *w, double *y, double *w0, int n, int nw0, double h, double *G);

void Est_Kernel(double *hatg, double *y, double *G, double *haty, double *w, int n, int p, int q, int nf, double *h, double *g, double *w0, int nw0);

void OLS(double *beta, double *x, double *y, int n, int p);

void Kernelh(double *x, double *kern, int n, double h0, int type);

void EstLinear(double *beta, double *x, double *z, double *y, double *weight, int n, int p, double h, int type);

void EstLogistic(double *beta, double *x, double *z, double *y, double *weight, int n, int p, double h, int maxstep, double eps, int type);

void EstPoiss(double *beta0, double *x, double *z, double *y, double *weight, int n, int p, double h, int maxstep, double eps, int type);

void EstLinearW(double *beta, double *residual, double *x, double *z, double *y, double *weight, int n, int p, double h, int type);

void EstLogisticW(double *beta, double *residual, double *x, double *z, double *y, double *weight, int n, int p, double h, int maxstep, double eps, int type);

void EstPoissW(double *beta0, double *residual, double *x, double *z, double *y, double *weight, int n, int p, double h, int maxstep, double eps, int type);

void _BandMatrix(double *sighalf, double *rho, int p, int T);

void _Toeplitz(double *sighalf, double *rho, int p);

#endif // RBS_H_INCLUDED
