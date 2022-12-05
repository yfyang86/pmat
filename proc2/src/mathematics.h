#ifndef _MATHEMATICS_H_
#define _MATHEMATICS_H_

#include <stdlib.h>
#include <stdio.h>
#include "basic-types.h"
#include <math.h>
#include "mm2.h"

#ifdef __linux
#define MARKER_st1_randgam
#endif

/*typedef struct
{
    size_t k;
    gsl_matrix *A;
    gsl_matrix *dB;
} gsl_bspline_deriv_workspace;
*/

/* important for compiling with cygwin */
#undef log2

#define RANDUNIT ((rand() / (RAND_MAX + 1.0)))
//#define RANDINT(MAX)                    ((Uint)round(((double)MAX) * (rand()/((double)RAND_MAX + 1.0))))
#define RANDINT(MAX) ((Uint)round((((double)MAX + 1) * (rand() / ((double)RAND_MAX + 1.0))) - 0.5))

#define MATRIX2D(X, NCOL, I, J) X[(I) * (NCOL) + (J)]
#define MATRIX3D(X, DIM_M, DIM_N, I, J, K) X[((I) * (DIM_M) + (J)) * DIM_N + (K)]
#define MATRIX4D(X, DIM_M, DIM_N, DIM_Z, I, J, K, L) X[(((I) * (DIM_M) + (J)) * DIM_N + (K)) * DIM_Z + (L)]
#define VECTOR(X, I) ((X)->elements[I])

#define INITMATRIX(X, SIZE, TYPESIZE) initArray(X, SIZE, TYPESIZE)
#define INITMATRIX2D(X, M, N, TYPESIZE) initArray(X, (M) * (N), TYPESIZE)
#define INITMATRIX3D(X, M, N, L, TYPESIZE) initArray(X, (M) * (N) * (L), TYPESIZE)
#define INITVECTOR(V)     \
    (V)->elements = NULL; \
    (V)->length = 0

#define RESIZEVEC(V, N) initArray(V, N, sizeof(vectorelem))
#define LENGTHVEC(V) ((V)->length)
#define SWEEPVEC(V)                        \
    {                                      \
        int m;                             \
        for (m = 0; m < LENGTHVEC(V); m++) \
            VECTOR(V, mi) = 0;             \
    }
#define APPENDVEC(S, V, E) appendvector(S, V, E)
#define SWAPVEC(I, J, V)                     \
    {                                        \
        vectorelem msv;                      \
        msv = (V)->elements[I];              \
        (V)->elements[I] = (V)->elements[J]; \
        (V)->elements[J] = msv;              \
    }
#define DESTRUCTVEC(S, V) destructVector((S), (V));
#define SWAPUINT(X, A, B) \
    {                     \
        Uint msu;         \
        msu = X[A];       \
        X[A] = X[B];      \
        X[B] = msu;       \
    }

#define REVERSEVEC(A, B, V)              \
    {                                    \
        Uint mi;                         \
        for (mi = 0; mi < (B - A); mi++) \
        {                                \
            SWAPVEC(A + mi, B - mi, V);  \
        }                                \
    }
#define EMPTYVEC(V) (V->elements == NULL)
#define MAX(A, B) (((A) >= (B)) ? (A) : (B))
#define DMAX(A, B) (double)MAX(A, B)

#define MAX3(A, B, C) MAX(MAX((A), (B)), (C))
#define MAX4(A, B, C, D) MAX(MAX(MAX((A), (B)), (C)), (D))
#define MAX5(A, B, C, D, E) MAX(MAX(MAX(MAX((A), (B)), (C)), (D)), (E))

#define MIN(A, B) (((A) <= (B)) ? (A) : (B))
#define DMIN(A, B) (double)MIN(A, B)
#define MIN3(A, B, C) MIN((MIN((A), (B))), (C))

#define MAX_DOUBLE 1e200
#define ABS(X) ((X) >= 0 ? (X) : -(X))
#define DABS(X) (double)ABS(X)
#define OVERLAP(A, B, C, D) (((A) >= (C) && (B) <= (D)) || \
                             ((A) <= (D) && (B) >= (D)) || \
                             ((B) >= (C) && (A) <= (C)))
#define CLOSEDINTERVAL(X, A, B) ((A) >= (X) && (X) <= (B))
#define OPENINTERVAL(X, A, B) ((A) > (X) && (X) < (B))
#define ISELEM(X, A, B) ((A) <= (X) && (X) <= (B))
#define expmax log(MAX_DOUBLE)

#ifndef FLT_EPSILON
#define FLT_EPSILON 1.192092896e-06
#define DBL_EPSILON 2.2204460492503131e-016
#endif

#ifndef FLT_MIN
#define FLT_MIN 1.175494e-38
#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

#ifndef M_LOGSQRT2PI
#define M_LOGSQRT2PI 0.918938533204673
#endif

#ifndef M_SQRT2PI
#define M_SQRT2PI 2.50662827463100005
#endif

#define EPSILON4 1.0e-4
#define EPSILON3 1.0e-3
#define EPSILON8 1.0e-8
#define EPSILON2 1.0e-14

#ifndef M_SQRT2
#define M_SQRT2 1.4142135623730950488016887
#endif

#ifndef ALLOCMEMORY
#include "memory.h"
#endif

#ifndef VECTORELEM
typedef Uint vectorelem;
#endif

typedef struct
{
    vectorelem *elements;
    Lint length;
} vector_t;

typedef struct
{
    Uint noofbreaks;
    Uint *breaks;
    double RSS;
    double LL;
    double BIC;
} breakpoints_t;

typedef struct
{
    double *l;
    Uint n;
} ecdf_t;

// folded normal: [author](0x9527eWZ5YW5n)
#ifdef MARKER_st1_randgam
double st1_randgam(double alpha, double beta);
#endif
inline double logit_f_ex(double x, double a)
{
    return log((x + a) * (1 - a) / (1 - (x + a) * (1 - a)));
}
double foldednomal(double *a, int m, double *b, int n);
double foldednomalpvalue(double *a, int m, double *b, int n, bool logit, double * samParam, double logittuning = 0.00001, bool nonZero=true, bool useBartlette=false);

// Orginal meteline:
void *initArray(void *, int, size_t);
double *coldel(void *, double *, Uint, Uint, Uint);
double *rowdel(void *, double *, Uint, Uint, Uint);
double scalar(double *x, double *y, Uint m);
double var(double *x, Uint n);
double rho(void *space, double *x, double *y, Uint m);
double gammaln(double x);
double kscdf(double x);
char IsFiniteNumber(double x);
double mannwhitneyPvalue(double u, Uint m, Uint n, double ***CDF, Uint maxm, Uint maxn, int nonZeroPairs, bool useFoldedNormal);
double mannwhitney(double *a, Uint m, double *b, Uint n, int *nonZeroPairs);
double mannwhitney(double *a, Uint m, double *b, Uint n, int *nonZeroPairs, bool useBartlett, double * samParam);
double *generateMannWhitneyCDF(Uint m, Uint n);
double ***generateMannWhitneyCDFMatrix(Uint maxm, Uint maxn);
void destructMannWhitneyCDFMatrix(double ***CDF, Uint m, Uint n);
double uniform_rand();
void randgauss(double *r1, double *r2);
double randgam(double alpha, double beta);
double randbeta(double alpha, double beta);
double rbeta_mv(double mu, double var);

#endif