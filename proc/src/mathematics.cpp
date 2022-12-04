
/*
 * mathematics.c
 * implemtation of various mathematical functions
 *
 * @author Steve Hoffmann
 * @date Wed 22 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 54 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-09-10 22:13:30 +0200 (Wed, 10 Sep 2008) $
 *
 *  Id: $Id: mathematics.c 54 2008-09-10 20:13:30Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/trunk/libs/mathematics.c $
 *
 */

#define __MFOLDED_NORMAL__
#define __FOLDED_NORMAL__

#include "mathematics.h"
#include "sort.h"
#include "mm2.h"

#include <float.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <vector>
#include <iostream>

 void *initArray(void *space, int size, size_t datatype) {
	void *ptr=NULL;

	/*dirty trick: sizeof(char) == 1*/
	ptr = ALLOCMEMORY(space, ptr, char, size*datatype);
	return ptr;
 }

/*---------------------------------- coldel ----------------------------------
 *
 * @brief delete column for matrix
 * @author Steve Hoffmann
 *
 */

double*
coldel (void *space, double *a, Uint m, Uint n, Uint d) {

	double *t;
	Uint	i,
			j=-1,
			k=0,
			l=0;

  t = (double*) INITMATRIX2D(space, m, (n-1), sizeof(double));

  for(i=0; i < m*n; i++) {
	if(i % n == 0) {
	  j++; k=0; l=0;
	}
	if(k++ != d) {
	  MATRIX2D(t, n-1, j, l++) = a[i];
	}
  }

  FREEMEMORY(space, a);
  return t;
}


/*---------------------------------- rowdel ----------------------------------
 *
 * @brief delete row from matrix
 * @author Steve Hoffmann
 *
 */

double*
rowdel (void *space, double *a, Uint m, Uint n, Uint d) {

	double *t;
	Uint	i,
			j=-1,
			k=0,
			l=-1;

  t = (double*) INITMATRIX2D(space, (n-1), m, sizeof(double));

  for(i=0; i < m*n; i++) {
	if(i % n == 0) {
	  j++; k=0;
	  l = (j != d) ? l+1 : l;
	}
	if(j != d) {
	  MATRIX2D(t, n, l, k++) = a[i];
	}
  }

  FREEMEMORY(space, a);
  return t;
}

/*----------------------------------- add ------------------------------------
 *
 * @brief componentwise addition of a to a vector of length m
 * @author Steve Hoffmann
 *
 */

double*
add(double *x, Uint m, double a) {
  Uint i;

  for(i=0; i < m; i++) {
    x[i] += a;
    //fprintf(stdout, "add: %f -> %f\n", x[i]-a, x[i]);
  }
  return x;
}


/*----------------------------------- mean -----------------------------------
 *
 * @brief calculate the arithmetic mean for a vector of length m
 * @author Steve Hoffmann
 *
 */

double
mean (double *x, Uint m) {
  Uint i;
  double sum=0;

  for (i=0; i < m; i++) {
    sum += x[i];
  }

  return sum/m;
}

/*---------------------------------- scalar ----------------------------------
 *
 * @brief calculate the scalar product of two vectors of length m
 * @author Steve Hoffmann
 *
 */

double
scalar (double* x, double *y, Uint m) {
  double  p=0;
  Uint 	i;

  for (i=0; i < m; i++) {
    //fprintf(stdout, "scal: %f*%f + %f\n", x[i], y[i], p);
    p += x[i]*y[i];
  }

  return p;
}


/*----------------------------------- cov ------------------------------------
 *
 * @brief get the covariance matrix (2x2) for two vectors of length m
 * @author Steve Hoffmann
 *
 */

double*
cov (void *space, double *x, double *y, Uint m) {
  double *c,
         xm,
         ym;

  c = (double*) INITMATRIX2D(space, 2, 2, sizeof(double));
  xm = mean(x, m);
  ym = mean(y, m);

  //fprintf(stdout, "meanx: %f\n", xm);
  //fprintf(stdout, "meany: %f\n", ym);

  /*center*/
  add(x, m, (-1)*xm);
  add(y, m, (-1)*ym);

  MATRIX2D(c, 2, 0, 0) = (double) scalar(x,x,m)/(m-1);
  MATRIX2D(c, 2, 0, 1) = MATRIX2D(c, 2, 1, 0) = (double) scalar(x,y,m)/(m-1);
  MATRIX2D(c, 2, 1, 1) = (double) scalar(y,y,m)/(m-1);

  return c;
}

/*----------------------------------- var ------------------------------------
 *
 * @brief get the variance
 * @author Steve Hoffmann
 *
 */

double
var (double *x, Uint n)
{
    int i;
    double m, r, sum=0;

    m=mean(x, n);
    for (i=0; i < n; i++) {
      r = x[i]-m;
      sum += (r*r);
    }

	return sum/n;
}


/*----------------------------------- rho ------------------------------------
 *
 * @brief calculate correlation $\rho$ for two vectors of length m
 * @author Steve Hoffmann
 *
 */

double
rho (void *space, double *x, double *y, Uint m) {
  double *cv;
  double ret;

  cv = cov(space, x, y, m);
  ret = (MATRIX2D(cv, 2, 0, 1)/sqrt(MATRIX2D(cv, 2, 0, 0)*MATRIX2D(cv, 2, 1, 1)));
  //fprintf(stdout, "%f %f %f -> %f\n",MATRIX2D(cv, 2, 0, 1), MATRIX2D(cv, 2, 0, 0), MATRIX2D(cv, 2, 1, 1), ret);
  FREEMEMORY(NULL, cv);

  return ret;
}


/*-------------------------------- lancozs gamma -------------------------------
 *
 * @brief lancozs approximation for the gamma function G=5 N=6+1
 * @author Steve Hoffmann
 *
 */

double
gammaln (double x)
{
    Uint i;
    double y,tmp,sum;
    double base;
    double G = 5.0;
    int N = 7;

    double lancozscoef[7] = { 1.000000000190015,
                             76.18009172947146,
                            -86.50532032941677,
                             24.01409824083091,
                             -1.231739572450155,
                              0.1208650973866179e-2,
                             -0.5395239384953e-5};
    y = x;

    sum = lancozscoef[0];
    for(i=1; i < N; i++) {
      sum += lancozscoef[i]/++y;
    }

    base = x + G + 0.5;
    tmp = (base) - log(base)*(x+0.5);
    return -tmp+log(M_SQRT2PI*sum/x);
}


/*-------------------------------- kscdf ---------------------------------
 *
 * @brief kscdf returns 1-x of cumulative distribution function of
 * K distribution
 * @author Steve Hoffmann
 *
 */

double kscdf (double x)
{

  unsigned int k;
  double sum = 0.0, tmp=0.0, old=0.0, val;
  double base, coeff = 1.0;

  base = -2.0*(x*x);

  for(k=1; k <= 100; k++) {
    tmp = exp(base*k*k);
    sum += coeff*tmp;
    if (tmp <= EPSILON3*old || tmp <= EPSILON8*sum) {
      val = (2.0*sum);
      return val;
    }
    coeff *= -1;
    old = tmp;
  }

  return 1.0;
}


/*------------------------------ IsFiniteNumber ------------------------------
 *
 * @brief check if x is a finite number
 * @author Steve Hoffmann
 *
 */


char
IsFiniteNumber(double x) {
  return (x <= DBL_MAX && x >= -DBL_MAX);
}


double
mannwhitneyGamma(double a, double x, double eps, int iter) {
  double an = 1.0 / a;
  double psum = an;
  double n = 0.0;
  while ( n < iter && fabs(an/psum) > eps) {
    n += 1;
    an *= x/(a+n);
    psum += an;
  }
  n = psum * exp((a * log(x)) - gammaln(a) - x );
  return n;
}

double
mannwhitneyPvalue(
  double u, Uint m, Uint n, double ***CDF, Uint maxm, Uint maxn,
  int nonZeroPairs
  ) {

double p;

#ifndef __FOLDED_NORMAL__
  double e, z, x;
  double u1;

  if(CDF != NULL && m <= maxm && n <= maxn) {
    if(m<n){
      p = CDF[m][n][u];
    }
    else {
      p = CDF[n][m][u];
    }
    p=1-p;
  } else {
    u1 = u;
    // z = (u1 - (m * n / 2.0)) / sqrt(m * n * (m + n + 1) / 12.0);
    //printf("mannwhitneypvalue start");
    // z = u1 / sqrt(nonZeroPairs * (nonZeroPairs+1) * (2 * nonZeroPairs + 1) / 6.0);
    z = (u1 -  (nonZeroPairs * (nonZeroPairs + 1)/4)) / sqrt(nonZeroPairs * (nonZeroPairs + 1) * (2 * nonZeroPairs + 1) / 24.0);
    //printf("%d\n",z);
    x = z / sqrt(2.0);
    if (fabs(x) > 40) {
      if(x > 0) {e=1;} else {e=-1;}
    }
    else {
      e = mannwhitneyGamma(0.5, x * x, 1.0e-15, 10000);
      if(x < 0) {e=-e;}
    }
    if (fabs(z) > 40) {
      if(z<0) p=0.0; else p = 0.5;
    }
    else {
      p = 0.5 * (1 + e);
    }
//     p = 1-p;
  }

//TODO: please check one-sided and two-sided test with p*=2 there are p-values > 1!
    p*=2;
#endif // __FOLDED_NORMAL__

/*in metseg.cpp:
                double ua = mannwhitney (cpg->groupA, noA, cpg->groupB, noB, nonZeroPairs);
		            double p= mannwhitneyPvalue(ua, noA, noB, nfo.MWU, MAXM, MAXN, *nonZeroPairs);
  Hence we only need to bypass the `ua` value as th p-value.
*/
p = (double) u;
//to be safe ;D
p = MIN(p, 1.0 - 1e-5); //TODO: for mannwhitneyGamma, please check one-sided and two-sided test with p*=2 there are p-values > 1!
    //printf("mannwhitney end");
return p;
}


/*------------------------------- mannwhitney --------------------------------
 *
 * @brief mann whitney u test
 * @author Steve Hoffmann
 *
 */
double mannwhitney(double *a, Uint m, double *b, Uint n, int *nonZeroPairs)
{
  double ans;

  if (m != n){
        printf("mannwhitney end\n");
        return 0;
    }
  #ifndef __FOLDED_NORMAL__
    // printf("mannwhitney start\na\n");
    // for(int i=0;i<m;i++){
    //   printf("%lf\n",a[i]);
    // }
    // printf("%d\n",m);
    // printf("mannwhitney start\nb\n");
    // for(int i=0;i<n;i++){
    //   printf("%lf\n",b[i]);
    // }
    // printf("%d\n",n);
    Uint  *sorted;
    double *ranks;
    double sumNegative = 0.0;
    double sumPositive = 0.0;
    double *c;

    double *temp, *positiveTemp;
    int posTemp;
    temp = (double *) ALLOCMEMORY(NULL, NULL, double, n);
    positiveTemp = (double *) ALLOCMEMORY(NULL, NULL, double, n);
    // 将a-b的结果存入temp数组中，丢弃0值，同时记录下temp的绝对值放入positiveTemp数组中
    posTemp = 0;
    for (int i = 0; i < n; i++)
    {
        temp[posTemp] = a[i] - b[i];
        positiveTemp[posTemp] = temp[posTemp] > 0 ? (temp[posTemp]) : (temp[posTemp] * -1);
        if (temp[posTemp] != 0)
        {
            posTemp++;
        }
    }
    // 因为丢弃了包含0的元素，更新数组长度
    n = posTemp;
    // printf("%d\n",n);
    //返回非0对个数
    *nonZeroPairs = posTemp;
    c = (double *) ALLOCMEMORY(NULL, NULL, double, n);
    ranks = (double *) ALLOCMEMORY(NULL, NULL, double, n);
    // 将c中的元素替换为temp的绝对值，开始排序，并对相同元素计算平均序号
    for (int i = 0; i < n; i++)
    {
        c[i] = positiveTemp[i];
    }

    sorted = quickSort(NULL, c, n, cmp_dbl, NULL);

    Uint j = -1;
    double entry = NAN;
    for (int i = 0; i < n; i++)
    {
        double rank = (double)i + 1;
        if (c[sorted[i]] != entry)
        {
            // 第j + 1个元素和第i个元素之间都是相同的值，求j~i之间的序号平均值作为序号值
            if (j < i - 1)
            {
                double m = 0;
                if (j > 0)
                {
                    m = j - 1;
                    m = (j * j + j) / 2;
                }
                double oldrank = i;
                oldrank = (oldrank * oldrank + oldrank) / 2;
                oldrank -= m;
                oldrank /= i - j;
                for (int z = j; z < i; z++)
                    ranks[sorted[z]] = oldrank;
            }
            j = i;
            ranks[sorted[i]] = rank;
        }
        entry = c[sorted[i]];
    }
    if (j < n - 1)
    {
        int i = n;
        double m = 0;
        if (j > 0)
        {
            m = j - 1;
            m = (j * j + j) / 2;
        }
        double oldrank = i;
        oldrank = (oldrank * oldrank + oldrank) / 2;
        oldrank -= m;
        oldrank /= i - j;
        for (int z = j; z < i; z++)
            ranks[sorted[z]] = oldrank;
    }


    // 根据temp中的元素的正负分别计算为正元素的rank累计和和父元素的rank累计和
    for (int i = 0; i < n; i++)
    {
        if (temp[i] > 0)
        {
            sumPositive += ranks[i];
        }
        if (temp[i] < 0)
        {
            sumNegative += ranks[i];
        }
    }
    //结果取绝对值
    //ans = sumNegative - sumPositive;
    ans = sumPositive;
    // printf("%lf\n",ans);
    ans = (ans < 0) ? (ans * -1) : ans;
    // for(int i=0;i<posTemp;i++){
    //   printf("%lf,",ranks[i]);
    // }
    //释放空间
    FREEMEMORY(NULL, c);
    FREEMEMORY(NULL, ranks);
    FREEMEMORY(NULL, sorted);
    FREEMEMORY(NULL, temp);
    FREEMEMORY(NULL, positiveTemp);
    // printf("mannwhitney end\n");
    //printf("%lf\n",ans);
#endif
#ifdef __FOLDED_NORMAL__
ans = foldednomalpvalue(a, m, b, n, true, 0.00001, true);
#endif // __FOLDED_NORMAL__
    return ans;
}



double*
generateMannWhitneyCDF(Uint m, Uint n) {
  double *C;
  long long unsigned int *S;
//  long long unsigned int minU, maxU, i, j, k, iter=0, u=0;
  long long unsigned int maxU, i, j, k, iter=0, u=0;

  assert(m <=n );
  maxU = m*n + (m*(m+1))/2;

//replace Umax with Umin, from here maxU <- Umin
//  minU = m*n - maxU;
//  maxU = MIN(maxU,minU);

  //setup matrix
  C = (double *) ALLOCMEMORY(NULL, NULL, double, maxU);
  S = (long long unsigned int *) ALLOCMEMORY(NULL, NULL, long long unsigned int, m);
  memset(C, 0, sizeof(double)*maxU);
  C[0] = 1;

  for(i=0; i < m; i++) {
    S[i]=i+1;
  }

  i =0;
  j = m-1;

  while(1) {

    if(iter > 0) {
      Uint maxrank = m+n;
      j = m-1;

      while(j > 0 && S[j] == maxrank) {
        maxrank--;
        j--;
      }
      //the highest pos has its max rank (=m+1)
      //all U are calculated
      if(j==0 && S[j] == maxrank)
        break;

      k = m-j-1;
      j = m-1;
      maxrank = m+n;

      while(j > 0 && S[j]==maxrank) {
        //this is the lowest element possible
        S[j] = S[j-k]+k+1;
        j--;
        k--;
        maxrank--;
      }
      S[j]++;
    }

    for(u=0, i=0; i < m; i++) {
      u += S[i];
    }

    u = maxU - u;
    C[u]++;
    iter++;
  }

  C[maxU-1] = C[maxU-1]/iter;

  for(i=maxU-2; i > 0; i--) {
    C[i] = C[i+1] + C[i]/iter;
  }

  C[0] = C[1];

  FREEMEMORY(NULL, S);

  return C;
}


/*----------------------- generateMannWhitneyCDFMatrix -----------------------
 *
 * @brief generate the CDF matrix
 * @author Steve Hoffmann
 *
 */

double***
generateMannWhitneyCDFMatrix(Uint maxm, Uint maxn)
{
  double ***M;
  Uint i, j;

  M = (double ***) ALLOCMEMORY(NULL, NULL, double**, maxm+1);

  for(i=0; i <= maxm; i++) {
    M[i] = (double **) ALLOCMEMORY(NULL, NULL, double*, maxn+1);
    for(j=i; j <= maxn; j++) {
      if(i >= 1) {
	//      fprintf(stderr, "generating matrix for %d,%d\n", i, j);
      M[i][j] = generateMannWhitneyCDF(i,j);
      }
    }
  }
  return M;
}


/*----------------------- destructMannWhitneyCDFMatrix -----------------------
 *
 * @brief destruct the MWU matrix
 * @author Steve Hoffmann
 *
 */

void
destructMannWhitneyCDFMatrix ( double ***CDF, Uint m, Uint n)
{
  Uint i,j;

  for(i=0; i <= m; i++){
    for(j=i; j <= n; j++){
      if(i >= 1)
      FREEMEMORY(NULL, CDF[i][j]);
    }
    FREEMEMORY(NULL, CDF[i]);
  }

  FREEMEMORY(NULL, CDF);
  return ;
}

double uniform_rand_range(double rangeLow, double rangeHigh) {

    double myRand = rand()/(1.0 + RAND_MAX);
    double range = rangeHigh - rangeLow ;
    double myRand_scaled = (myRand * range) + rangeLow;
    return myRand_scaled;
}

double uniform_rand() {
    return uniform_rand_range(0.0,1.0);
}

void randgauss(double *r1, double *r2)
{

  double u1, u2, r, f;

  do {
//    u1 = 2.0 * ((double)rand() / (double)RAND_MAX) -1.0;
//    u2 = 2.0 * ((double)rand() / (double)RAND_MAX) -1.0;
    u1 = 2.0 * uniform_rand()  -1.0;
    u2 = 2.0 * uniform_rand()  -1.0;

    r = u1*u1 + u2*u2;
  } while(r >= 1.0 || r ==0.0);

  f = sqrt(-2.0*log(r)/r);
  *r1 = f*u1;
  *r2 = f*u2;

  return ;
}

/*--------------------------------- randgam ----------------------------------
 *
 * @brief generate random gamma
 * @author Steve Hoffmann
 * http://www.hongliangjie.com/2012/12/19/how-to-generate-gamma-random-variables/
 *
 */

double randgam (double alpha, double beta)
{

  double x, z1=0, z2=0, u, v, d, c;
  char phase = 0;

  assert(alpha > 0);

  if(alpha >= 1.0) { // ? > or >=

    d = alpha-1.0/3.0;
    c = 1.0/sqrt(9.0*d);
    while(1) {

      if(phase == 0) {
        randgauss(&z1, &z2);
        phase = 1;
      } else {
        z1 = z2;
        phase = 0;
      }

      if (z1 > -1.0/c) {
        v = (1+c*z1);
        v = v*v*v;
//        u = (double)rand() / (double)RAND_MAX; // U(0,1)
        u = uniform_rand();
        if(log(u) < 0.5*(z1*z1)+d-d*v+d*log(v)) break;
      }
    }
    x = d*v/beta;

  } else {

    x = randgam(alpha+1, beta);
//    u = (double)rand() / (double)RAND_MAX; // U(0,1)
    u = uniform_rand();

    x = x * pow(u, 1.0/beta);
  }

  return x;
}


/*--------------------------------- randbeta ---------------------------------
 *
 * @brief random beta
 * @author Steve Hoffmann
 *
 */

double randbeta (double alpha, double beta)
{

  double x,y,u,v;

  if (alpha <= 1.0 && beta <= 1.0){
    while(1) {
//      u= (double)rand() / (double)RAND_MAX;
//      v= (double)rand() / (double)RAND_MAX;
      u= uniform_rand();
      v= uniform_rand();
      x = pow(u, 1.0/alpha);
      y = pow(v, 1.0/beta);
      if(x+y <= 1.0) {
        return x/(x+y);
      }
    }
  } else {
    x = randgam(alpha, 1);
    y = randgam(beta,  1);
    return x/(x+y);
  }
}

double rbeta_mv(double mu, double var) {
    if(mu <= 0.000001)
        mu=0.000001;
    if(var <= 0.000001)
        //var=0.0000000001;
    //if( var <= 0.0000000001)
        return mu;
//  fprintf(stdout,"no return\n");
  double x = ((mu * (1-mu))/var)-1;
  double alpha = x * mu;
  double beta = (1-mu) * x;
//  fprintf(stderr,"#alpha: %f\tbeta: %f\n",alpha,beta);

  if(alpha <= 0.000001)
        alpha=0.000001;
    if(beta <= 0.000001)
        beta=0.000001;
 // fprintf(stderr,"#alpha: %f\tbeta: %f\n",alpha,beta);

  //double alpha = ((1 - mu) / var - 1 / mu) * (mu *mu);
  //double beta = alpha * (1 / mu - 1);
  //fprintf(stdout,"MuVar\t%f\t%f\n",mu,var);
 // fprintf(stdout,"AlphaBeta\t%f\t%f\n",alpha,beta);
//  double ret = rbeta(alpha, beta);
  double ret = randbeta(alpha, beta);

  return ret;
}



/*
********************************************************     *******
*       Folded normal: [author](0x9527eWZ5YW5n)    ******* *******
********************************************************     *******
*/

/**
 * folder normal -2LLR test stats
 * @param a a double array: samples A
 * @param m an integer: sample size A
 * @param b a double array: samples B
 * @param n an integer: sample size B
 * @return p-value (0 is error)
 */
double foldednomal(double *a, int m, double *b, int n ){
    if (m != n){return 0.;}
    std::vector<double> x(m);
    for (size_t i = 0; i < x.size(); i++) x[i] = fabs(a[i] - b[i]);
    double re = LRFstatLLR(x);
#ifdef DEBUG
    std::cout<<"-2LLR = "<<re<<std::endl;
#endif // DEBUG
    return re;
}

/**
 * folder normal -2LLR test, P-value(FN, abs(a-b))
 * @param a a double array: samples A
 * @param m an integer: sample size A
 * @param b a double array: samples B
 * @param n an integer: sample size B
 * @param logit t/f: use Logit transformation or not:  log((x+a)*(1-a)/(1-(x+a)*(1-a)))
 * @param logittuning double:  Logit transformation tuning size
 * @return p-value (0 is error)
 */
double foldednomalpvalue(double *a, int m, double *b, int n, bool logit=true, double logittuning=0.00001, bool nonZero = true){
    bool useBartlette = true;
    // config.txt
    double samParam[6]={0.60105772, 0., -4.0224, 0.89, 1., 1.};
    //double samParam[6]={0.5, 0., 0, 0., 1., 1.};

    if (m != n)
    {
        return 0.;
    }
    std::vector<double> x(m);
    int nonZeroSize = 0;
    double s;
    double ABSadifb;
    double corr_v[6] = {0., 0., 0., 0., 0., 0.};
    for (size_t i = 0; i < x.size(); i++)
    {
        ABSadifb = fabs(a[i] - b[i]);
        s = (logit ? fabs(logit_f_ex(a[i], logittuning) - logit_f_ex(b[i], logittuning)) : ABSadifb);
        if (nonZero)
        {
            if (ABSadifb > 1e-15)
            {
                x[nonZeroSize] = s;
                nonZeroSize++;
            }
        }
        else
        {
            x[i] = s;
        }
    }
    if (nonZero)
        x.resize(nonZeroSize);
    std::vector<double> re;
    if (useBartlette){
            for (size_t i = 0; i < m ; i++){
                corr_v[0] += a[i] * b[i];
                corr_v[1] += a[i] * a[i];
                corr_v[2] += b[i] * b[i];
                corr_v[3] += a[i] ;
                corr_v[4] += b[i] ;
            }
            corr_v[5] = corr_v[0] / sqrt( (corr_v[1] * corr_v[1] - corr_v[3]*corr_v[3]/((double) n )) * (corr_v[2] * corr_v[2] - corr_v[4]*corr_v[4]/((double) n ))  );
    }
    int indx_tunning_param[4] = {-1, 0 , 1, 2};
    /** TODO:
     * 0:       Use EM as the initial value. If there are extreme values in the input vector, this may fail.
     * -1:      negative initial value tuning
     * 1/2/...: Positive  initial value tuning
     */
    for (int i = 0; i <= 3; i++)
    {
        if (useBartlette){
        re = LRFstatNative(x, indx_tunning_param[i], corr_v[5], useBartlette, samParam);
        }else{
            re = LRFstatNative(x, indx_tunning_param[i]);
        }

        if (re[3] > 0 && !(isnan(re[2])))
        {
            // std::cerr << "p-value:\t" << i << "\t" << re[2] << std::endl;
            break;
        }
    }
    if (re[3] < 0)
    {
        re[2] = 1.;
        // std::cerr << "EM-Algorithm Not Converge!!!" << std::endl;
    }

    // re = LRFstat(x);
#ifdef DEBUG
    std::cout << "Result: \nMLE\t[" << re[0] << ", " << re[1] << "p-value:\t" << re[2] << std::endl;
#endif // DEBUG
    return re[2];

}