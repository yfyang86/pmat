#include <stdio.h>
#include <math.h>
#include "mm2.h"  //add by xp
#include <iostream>  //add by xp
#include <gsl/gsl_cdf.h>  //add by xp
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#ifdef USE_ARMADILLO
#include <armadillo>
#endif

#define RETRY 100
#define MIN_ITER 50
#define MAX_ITER 10000



double sam_function(double * y, double * para, int n, double sy2){
    double s = 0;
    double mu_pi = 3.1415926535897932;
    size_t i = 0;
    // n = len(y)
    for (;i<n; i++){
        s = s + log(cosh(para[0] * y[i] / para[1]));
    }
    s = -((double) n)* log(2.0/mu_pi/para[1])/2. +
        (double) n * para[0] * para[0] / 2.0/para[1] + sy2/2.0/para[1] -
        s;
    return s;
}

// In GSL, param is the constant variable one uses in the likelihood
// const gsl_vector * y, meanwhile, is your parameters to be updated ...
double gsl_sam_function(const gsl_vector * y, void * param){
    double para[2] = {0.,0.};
    para[0] = gsl_vector_get(y,0); //mean
    para[1] = gsl_vector_get(y,1); //sd
    return(
        sam_function(
                ( (struct params *) param)->y,
                para,
                ( (struct params *) param)->n,
                ( (struct params *) param)->sy2)
        );
}

void solver_sam_function(double * y, int n, double * param_init, double *result, size_t control){
    gsl_set_error_handler_off();
    double sy2=0;
    for(size_t ii =0;ii<n;ii++) sy2+=y[ii]*y[ii];
    size_t iter=0;
    double size=0;
    int status;
    size_t status_loop = 0;
    struct params p = {y, n, sy2};
    gsl_multimin_function f = {&gsl_sam_function, 2, &p};
    gsl_vector *ss, *x;
    double x_init[2] = {param_init[0], param_init[1]};
    x = gsl_vector_alloc(2);
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;//gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    // Initialization
    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);
    /* Set initial step sizes to 1 */
    ss = gsl_vector_alloc (2);
    gsl_vector_set_all (ss, 1.0);

    s = gsl_multimin_fminimizer_alloc (T, 2);
    gsl_multimin_fminimizer_set (s, &f, x, ss);
    #ifdef DEBUG
        int steps__=1;
    #endif
        try{
            do{ status_loop=0;
                iter++;
                status = gsl_multimin_fminimizer_iterate(s);

                size = gsl_multimin_fminimizer_size (s);
                status = gsl_multimin_test_size (size, 1e-6);
                if (control==1){
                    if (status == GSL_SUCCESS) printf ("converged to minimum at\n");
                    printf ("%10d %10.6e %10.6e f() = %10.6f size = %.6f\n",
                            (int)iter,
                            gsl_vector_get(s->x, 0),
                            gsl_vector_get(s->x, 1),
                            s->fval, size);
                }
                if (iter <= MIN_ITER) status_loop = 1;
                if (iter < MAX_ITER && status == GSL_CONTINUE) status_loop = 1;
                if (iter > MAX_ITER) status_loop = 0;
                
                #ifdef DEBUG
                steps__ ++;
                #endif
                } while ( status_loop == 1);
        }
        catch (...){
            printf ("Infinite");
        }

  result[0] = gsl_vector_get(s->x, 0);
  result[1] = gsl_vector_get(s->x, 1);
  // result[2] = (status == GSL_SUCCESS ? 1 : -1);
  result[2] = 1; // bypass
  if (status != GSL_SUCCESS) std::cerr<<"Cov-Issue"<<std::endl;

#ifdef DEBUG
std::cout <<"[DEBUG] Loops: "<< steps__ << std::endl;
#endif // DEBUG
  

    gsl_vector_free (x);
    gsl_vector_free (ss);
    gsl_multimin_fminimizer_free (s);
}

/**
 * folder normal likelihood function: $FoldedN(\mu, \sigma^2)$
 * @param x a double array: samples
 * @param mu double constant: normal mean $\mu$
 * @param sigma2 double constant: normal $sigma^2$
 * @param n an integer: sample size
 * @return The test results
 */
double LRF(double *x, double mu, double sigma2, size_t n){
    double pi=3.1415926535897932;
    double re0=0.,re1=0.,re=0;
    size_t i = 0;
    if (sigma2<1e-14) {
        sigma2=1;
    }
    for(;i<n;i++){
        re0 = re0 + x[i] * x[i];
        re1 = re1 + log(cosh(mu*x[i]/sigma2));
    }
    re = (double)n/2.*log(2/pi) -  (double)n/2.*log(sigma2) - (re0 + n*mu*mu)/sigma2/2 + re1;
return re;
}

/**
 * folder normal -2LLR test, P-value
 * @param y a double array: samples
 * @param n an integer: sample size
 * @return p-value
 */
double LRFstat(double* y, size_t n) {
    double mean = 0.;
    double sd = 0.;
    double mle[3] = { 0., 0., 0. };
    double lrt;
    for (size_t i = 0; i < n; i++) mean += y[i];
    mean = mean / n;
    for (size_t i = 0; i < n; i++) sd += (y[i] - mean) * (y[i] - mean);
    sd = sqrt(sd / (n + 1.));
    double param_init[2] = { mean, sd };
    try {
        solver_sam_function(y, n, param_init, mle, 0);
    }
    catch (...) {
    }
    if (mle[2] > -0.5) {
        double meany2 = 0;
        for (size_t i = 0; i < n; i++) {
            meany2 += y[i] * y[i];
        }
        meany2 = meany2 / n;
        lrt = 2 * (LRF(y, mle[0], mle[1], n) - LRF(y, 0., meany2, n));
        lrt = lrt < 0 ? 0 : lrt;
    }
    else {
        lrt = -1.;
    }
   // return lrt;
    double chisq_Q = 1 - chisq_mix_P(lrt, 0.5);
    return chisq_Q;
}

/**
 * folder normal -2LLR test, P-value
 * @param y a std::vector<double> array: samples
 * @return std::vector<double> = { MLE[0], MLE[1], p-value}
 */
std::vector<double> LRFstat(std::vector<double> input_array){
    std::vector<double> result(3,0.);
    double * y = &input_array[0];
    int n = input_array.size();
        double mean = 0.;
    double sd = 0.;
    double mle[3] = {0., 0., 0.};
    double lrt;
    for(size_t i = 0; i<n; i++) mean += y[i];
    mean = mean/n;
    for(size_t i = 0; i<n; i++) sd += (y[i]-mean)*(y[i]-mean);
    sd = sqrt(sd/(n+1.));
    double param_init[2] = {mean, sd};
        try{
            solver_sam_function(y, n, param_init, mle, 0);
        }catch(...){
        }
    if(mle[2] > -0.5){
        double meany2 = 0;
        for (size_t i = 0; i<n;i++){
            meany2 += y[i]*y[i];
        }
        meany2 = meany2/n;
        lrt = 2 *(LRF(y, mle[0], mle[1], n)-LRF(y, 0., meany2, n));
        if (lrt<0.) std::cerr<< "LRT issue" << std::endl; 
        lrt= lrt<0.?0.:lrt;
    }else{
        lrt = -1.;
    }
    result[0] = mle[0];
    result[1] = mle[1];
    result[2] = 1 - chisq_mix_P(lrt, 0.5);
    return result;
}

/**
 * folder normal -2LLR test
 * @param y a std::vector<double> array: samples
 * @return double: -2LLR Test stats (dim=1)
 */
double LRFstatLLR(std::vector<double> input_array){
    double * y = &input_array[0];
    int n = input_array.size();
        double mean = 0.;
    double sd = 0.;
    double mle[3] = {0., 0., 0.};
    double lrt;
    for(size_t i = 0; i<n; i++) mean += y[i];
    mean = mean/n;
    for(size_t i = 0; i<n; i++) sd += (y[i]-mean)*(y[i]-mean);
    sd = (n >1 ? sqrt(sd/(n-1.)) : 1.0e-10);
    double param_init[2] = {mean, sd};
        try{
            solver_sam_function(y, n, param_init, mle, 0);
        }catch(...){
        }
    if(mle[2] > -0.5){
        double meany2 = 0;
        for (size_t i = 0; i<n;i++){
            meany2 += y[i]*y[i];
        }
        meany2 = meany2/n;
        lrt = 2 *(LRF(y, mle[0], mle[1], n)-LRF(y, 0., meany2, n));
        lrt= lrt<0?0:lrt;
    }else{
        lrt = -1.;
    }
    return lrt;
}


/*
f(x) MUST be monotonic
*/
double f_uniroot(double f(double,std::vector<double>), 
                 double start, 
                 double end, 
                 std::vector<double> paramas, 
                 bool debug){
    // size_t MAX_ITER = 2000;
    double split_lambda = 0.618;
    double split_lambda2 = 0.382;
    double tol = 1e-12;
    double f_l, f_r, f_l0;
    size_t iter =0;
    double tmp = end; 
    if (end < start) {end = start; start = tmp;}
    f_l = f(start, paramas);
    f_r = f(end, paramas);
    do{        
        if ( (MSIGN(f_l) ) * (MSIGN(f_r)) > 0.5) break;
        tmp = split_lambda*start + split_lambda2*end;
        f_l0 = f(tmp, paramas);
        if (fabs(f_l0) < tol) {iter = 10000; break;}
        if ( ( MSIGN(f_l0))  * (MSIGN(f_r)) > 0.5) {
            end = tmp;
            f_r = f_l0;
            if (debug) printf("+\t");
        }else{
            start = tmp;
            f_l = f_l0;
            if (debug) printf("-\t");
        }
        if (debug) printf("Iter %3d f(%4f,%4f)=(%4f, %4f).\n", (int) iter, start, end, f_l, f_r);
        iter++;
    }while (iter < MAX_ITER);
    if (iter < 2) printf("Same Sign, choose another initial start point.\n");
    return(end);
}

inline double chisq_mix_P(double x, double ratio){
    x = (x < 0 ? 0:x);
    ratio = (ratio > 1 ? 1:ratio);
    ratio = (ratio < 0 ? 0:ratio);
    // The cumulative distribution function for the lower tail P(x) 
    // is defined by the integral
    return( x < 0 ? 0: (ratio + (1 - ratio) * gsl_cdf_chisq_P(x, 1.)) );
}

void FoldedNormNativeSolverL1(double *param, double *Df, double *targetV, double *Vf)
{
    double s = param[0];
    double sigma = param[1];
    // sigma = (sigma > 0 ? sigma : -sigma);
    double exps2 = exp(-s * s / 2.);
    double PT1 = 0.797884560802865 * exps2;
    double PT2 = 1. - erfc(s * 0.707106781186547);
    double ds = 0.398942280401433 * exps2;            // dnorm(s)
    double vp_slot_1 = sigma * PT1 + s * sigma * PT2; // muY
    double a = -s * sigma * PT1 + sigma * PT2 + 2 * s * sigma * ds;
    double b = PT1 + s * PT2;
    double c = 2. * s * sigma * sigma - 2. * vp_slot_1 * a;
    double d = 2. * (s * s + 1.) * sigma - 2. * vp_slot_1 * b;
    double detDf = a * d - b * c;
    // detDf = (fabs(detDf) < 1e-15 ? 1.:detDf);
    //  D (target - f)
    Df[0] = -d / detDf;
    Df[1] = b / detDf;
    Df[2] = c / detDf;
    Df[3] = -a / detDf;

    Vf[0] = targetV[0] - vp_slot_1;
    Vf[1] = targetV[1] - ((1. + s * s) * sigma * sigma - vp_slot_1 * vp_slot_1);
}

/**
 *@title: Newton routine to solve (muy, sy2) = (u/s, s) for Folded Normal
 *@describe: Moment Estimator!
 *@param: iter int iteration
 *@param: v1 muy
 *@param: v2 sy
 *@return array [u, s]m i.e. [u/s * s, s]
 *@test: [bypass]
 **/
std::vector<double> FoldedNormNativeSolver(int iter, double v1, double v2)
{
    double targetV[2] = {0., 0.};
    double param[2] = {0., 0.};
    double Df[4] = {0., 0., 0., 0.};
    double Vf[2] = {0., 0.};

    targetV[0] = v1;
    targetV[1] = v2;

    param[1] = sqrt(v2);
    param[1] = (param[1] > 1e-12 ? param[1] : 1);
    param[0] = v1 / param[1];

    iter = (iter > 0 ? iter : 0);

    while (iter > 0)
    { // FoldedNormNativeSolverL1(double *param, double *Df, double *targetV, double *Vf)
        FoldedNormNativeSolverL1(param, Df, targetV, Vf);
#ifdef __DEBUG__L1__
        std::cerr << "cMEAN:\t" << param[0] * param[1] << ",\tcSTD.ERR\t" << param[1] << std::endl;
#endif
        param[0] = param[0] - (Df[0] * Vf[0] + Df[1] * Vf[1]);
        param[1] = param[1] - (Df[2] * Vf[0] + Df[3] * Vf[1]);
        iter--;
    }

    std::vector<double> re(2);
    re[0] = param[0] * param[1];
    re[1] = param[1];

    return re;
}


std::vector<double> LRFstatNative(std::vector<double> y, int init_scale, double rho, bool useBartlett, double * samParam){
    std::vector<double> re = LRFstatNative(y,init_scale);
    if (!useBartlett){
        return re;
    }

    double nsize = (double) y.size();
    double mix = samParam[0]+samParam[1]*rho+samParam[2]*pow(nsize, -samParam[3]);
    re[2] = 1 - pchibarsq(re[4], 1., 1. - mix, true, false);

    return re;

}


/**
 * folder normal -2LLR test, P-value
 * @param y a std::vector<double> array: samples
 * @return std::vector<double> = { MLE[0], MLE[1], p-value}
 */
std::vector<double> LRFstatNative(std::vector<double> x, int init_scale){
double tol = 1e-09;
size_t n = (size_t) x.size();
double m = 0.;
for(size_t i =0; i < x.size(); i++) m += x[i];
m = m/(double)n;
std::vector<double> x2(n);
for(size_t i =0; i < n; i++) x2[i] = x[i] * x[i];
double sx2 = 0.;
for(size_t i =0; i < x2.size(); i++) sx2+=x2[i];

double a = sx2/(double) n - m*m;
double fc01 = (0.95 + 0.9765*(double) init_scale);
std::vector<double> param_inner;
std::vector<double> y(n);
std::vector<double> tanhy(n);
std::vector<double> coshy2(n);
std::vector<double> aold(2);
std::vector<double> anew(2);
double derm ;
double ders ;
double derm2;
double ders2;
double derms;
double sum_x_tanhy = 0.;
double sum_y_tanhy = 0.;
double sum_y_x_coshy2 = 0.;
double sum_x2_coshy2 = 0.;
double sum_y2_coshy2 = 0.;
/*Newton Raphson Iterarion Settings: lr = 1.*/
int loop_counter = 0;
const int MINI_ITER = 5;
const int MAX_ITER__ = 500;
double lr = 1.;
double lrBASE = 1.;
int cnt_retry = 1;

fc01 = (fc01 > 0.5) ? fc01:0.5;

/*Init parameters:
sx2:  sum(x^2)
m  :  mean(x)
a  :  mean(x^2) or sd(x)^2;  
for H0: mu = 0; \har{\sigma^2} = mean(x^2)
without constraint, a proper estimation is (mean, sd^2)
*/

if (init_scale == 0)
{
        param_inner = FoldedNormNativeSolver(10, m, a);
        m = param_inner[0];
        a = param_inner[1];
#ifdef __DEBUG__L1__
        std::cerr << "MEAN:\t" << m << ", STD.ERR\t" << a << std::endl;
#endif
}

/*
if ( init_scale != 0){
    m = m * fc01;
    a = fc01 * sx2/(double) n;
}
*/




for(size_t i =0; i < n; i++) {
    y[i] = m *x[i]/a;
    tanhy[i] = tanh(y[i]);
    coshy2[i] = cosh(y[i]);  coshy2[i] = 1./coshy2[i]/coshy2[i];
}

/**
 * Update Routine: Newton Raphson
 * x[n] <- x[n-1] - H(f(x[n-1]))  \Lambda f(x[n-1])
 */
for(size_t i =0; i < n; i++) {
    sum_x_tanhy += x[i] * tanhy[i]; 
    sum_y_tanhy += y[i] * tanhy[i]; 
    sum_x2_coshy2 += x2[i] * coshy2[i];
    sum_y2_coshy2 += y[i] * y[i] * coshy2[i];
    sum_y_x_coshy2 += y[i] * x[i] * coshy2[i];
}
derm  = -(double)n * m/a + sum_x_tanhy/a;
ders  = -(double)n/2./a + ( sx2 + (double) n * m*m) /2./a/a - sum_y_tanhy / a;
derm2 = -(double) n /a + sum_x2_coshy2 /a/a; 
ders2 = (double) n/2./a/a - (sx2 + (double) n * m * m )/a/a/a + 2.* sum_y_tanhy /a/a + sum_y2_coshy2/a/a;
derms = (double) n * m /a/a - sum_x_tanhy /a/a - sum_y_x_coshy2 /a/a;
aold[0] = m; 
aold[1] = a;
anew[0] = aold[0] - lrBASE*lr*(ders2 * derm- derms * ders)/(derm2 * ders2 - derms*derms) ;
anew[1] = aold[1] - lrBASE*lr*((-derms * derm + derm2 * ders)/(derm2 * ders2 - derms*derms));


while(loop_counter < MAX_ITER__ && (fabs(anew[0] - aold[0]) + fabs(anew[1] - aold[1])  > tol || loop_counter > MINI_ITER)) {

if (isnan(anew[0]) || isnan(anew[1]) ) break;

#ifdef __DEBUG__L3
        std::cerr<<"DERM: "<<derm<<'\t'<<ders<<'\t'<<derm2<<'\t'<<ders2<<std::endl;
#endif

#ifdef __DEBUG__L2
        std::cerr<<"iter: "<<loop_counter<<'\t'<<anew[0]<<'\t'<<anew[1]<<std::endl;
#endif
        loop_counter++;
        m = anew[0];
        a = anew[1];
        aold = anew;
        for(size_t i =0; i < n; i++) {
            y[i] = m * x[i]/a;
            tanhy[i] = tanh(y[i]);
            coshy2[i] = cosh(y[i]);
            coshy2[i] = 1./coshy2[i]/coshy2[i];
        }
        
        sum_x_tanhy = 0.;
        sum_y_tanhy = 0.;
        sum_y_x_coshy2 = 0.;
        sum_x2_coshy2 = 0.;
        sum_y2_coshy2 = 0.;
        for(size_t i =0; i < n; i++) {
            sum_x_tanhy += x[i] * tanhy[i]; 
            sum_y_tanhy += y[i] * tanhy[i]; 
            sum_x2_coshy2 += x2[i] * coshy2[i];
            sum_y2_coshy2 += y[i] * y[i] * coshy2[i];
            sum_y_x_coshy2 += y[i] * x[i] * coshy2[i];
        }
            
        derm  = -(double)n * m/a + sum_x_tanhy/a;
        ders  = -(double)n/2./a + ( sx2 + (double) n * m*m) /2./a/a - sum_y_tanhy / a;
        derm2 = -(double) n /a + sum_x2_coshy2 /a/a; 
        ders2 = (double) n/2./a/a - (sx2 + (double) n * m * m )/a/a/a + 2.* sum_y_tanhy /a/a + sum_y2_coshy2/a/a;
        derms = (double) n * m /a/a - sum_x_tanhy /a/a - sum_y_x_coshy2 /a/a;
    
        if (loop_counter % 10 == 0)lr = (1. - (double) loop_counter /(double) MAX_ITER__)/2. + 0.5;

        anew[0] = aold[0] - lrBASE*lr* (ders2 * derm- derms * ders)/(derm2 * ders2 - derms*derms) ;
        anew[1] = aold[1] - lrBASE*lr* ((-derms * derm + derm2 * ders)/(derm2 * ders2 - derms*derms));

        if (isnan(anew[0]) || isnan(anew[1])) {
            //std::cerr<<"Try new Learning rate"<<std::endl;
            lrBASE = lrBASE /1.618;
            anew[0] = aold[0];
            anew[1] = aold[1];
            if (++cnt_retry>5) break;
        }
    }

    m = anew[0];
    a = anew[1];
    double loglik = 0.;
    double sum_x2 = 0.;
    for (size_t i = 0; i < n; i++){
        sum_x2 += x2[i];
        loglik += log(cosh(y[i]));
    }

    loglik += (double) n/2. *log(2./a/3.1415926535897932) - 0.5 * (double) n * m * m /a - 0.5 * sum_x2/a;
    double loglik_null = (double) n*0.5 *log(2. * (double)n / 3.1415926535897932 /sx2) - (double) n / 2.;
    double lrt;

    lrt = 2. * (loglik - loglik_null);
    //iters = i, loglik = loglik, param = anew
    std::vector<double> result(5);
    result[0] = anew[0];
    result[1] = anew[1];

    result[2] = 1 - chisq_mix_P(lrt, 0.5);
    result[3] = 1.0;
    if (lrt >= -1.e-7) lrt = fabs(lrt);
    result[4] = lrt; 

    if (lrt < -1.e-7  || isnan(anew[0]) || isnan(anew[1])) {
        result[2] = 0.99999;
        result[3] = -1.;
        //std::cerr<< "LRT issue" << std::endl; 
        }
#ifdef __DEBUG__
        std::cout<<"LogLik\t"<<loglik<<'\n'<<"LogLik_NULL\t"<<'\t'<<loglik_null<<"\n-2llr:\t"<<lrt<<std::endl;
#endif
    return result;
}

std::vector<double>  pchibarsq(std::vector<double> pvec, double df, double mix, bool lower_tail, bool log_p){
    std::vector<double> result_vec = pvec;
    std::vector<double> c1 = pvec;
    std::vector<double> c2 = pvec;
    bool df_eq_1 = fabs(df - 1.0) < 0.1;
    if (df_eq_1){
        if (lower_tail) {
            for (size_t i =0; i<c1.size();i++) c1[i] = 1.;
            } else {
                for (size_t i =0; i<c1.size();i++) c1[i] = 0.;
                }
    }else{
        for (size_t i =0; i<c1.size();i++) c1[i] = ( lower_tail ? gsl_cdf_chisq_P(pvec[i], df - 1.) : 1. - gsl_cdf_chisq_P(pvec[i], df - 1.));
    }

    for (size_t i =0; i<c2.size();i++) c2[i] = ( lower_tail ? gsl_cdf_chisq_P(pvec[i], df) : 1. - gsl_cdf_chisq_P(pvec[i], df));

    for (size_t i =0; i<c2.size();i++) result_vec[i] = mix * c1[i] + (1.-mix) * c2[i];

    if (log_p)  for (size_t i =0; i<c2.size();i++) result_vec[i] = log(result_vec[i]);

return result_vec;
}

double pchibarsq(double pvec, double df, double mix, bool lower_tail, bool log_p){

    double result_vec = pvec;
    double c1 = pvec;
    double c2 = pvec;
    bool df_eq_1 = fabs(df - 1.0) < 0.1;
    if (df_eq_1){
        c1 = (lower_tail?1.:0. );
    }else{
        c1 = ( lower_tail ? gsl_cdf_chisq_P(pvec, df - 1.) : 1 - gsl_cdf_chisq_P(pvec, df - 1.));
    }

    c2  = ( lower_tail ? gsl_cdf_chisq_P(pvec, df) : 1 - gsl_cdf_chisq_P(pvec, df));

    result_vec = mix * c1 + (1.-mix) * c2;

    if (log_p) result_vec = log(result_vec);

return result_vec;
}

#ifdef USE_ARMADILLO

std::vector<double> LRFstatNative_arma(std::vector<double> inputVec, int init_scale){
const arma::mat& x(inputVec);
int n = x.n_elem;
double m = arma::accu(x) / n;
arma::mat x2 = arma::pow(x, 2.0);
double sx2 = arma::accu(x2);
double a = sx2 / n - std::pow(m, 2.0);

if ( init_scale > 0){
    m = 0.1*(double) init_scale;
    a = sx2/(double) n*(0.9 + m/2.);
}

arma::mat y = m * x / a;
arma::mat tanhy = arma::tanh(y);
arma::mat coshy2 = 1.0 / arma::pow(arma::cosh(y), 2.0);
double derm = -n * m / a + arma::accu(x % tanhy) / a;
double ders = -n / 2.0 / a + (sx2 + n * std::pow(m, 2.0)) / 2.0 / std::pow(a, 2.0) - arma::accu(y % tanhy) / a;
double derm2 = -n / a + arma::accu(x2 % coshy2) / std::pow(a, 2.0);
double ders2 = n / 2.0 / std::pow(a, 2.0) - (sx2 + n * std::pow(m, 2.0)) / std::pow(a, 3.0) + 2.0 * arma::accu(y % tanhy) / std::pow(a, 2.0) + arma::accu(arma::pow(y, 2.0) % coshy2) / std::pow(a, 2.0);
double derms = n * m / std::pow(a, 2.0) - arma::accu(x % tanhy) / std::pow(a, 2.0) - arma::accu(y % x % coshy2) / std::pow(a, 2.0);
arma::vec aold(2, arma::fill::zeros);
aold(0) = m;
aold(1) = a;
arma::vec tmp = aold;
tmp(0) = ders2 * derm - derms * ders;
tmp(1) = -derms * derm + derm2 * ders;
tmp = tmp / (derm2 * ders2 - std::pow(derms, 2.0));
arma::mat anew = aold - tmp;
int i = 2;
const int MINI_ITER = 10;
const int MAX_ITER__ = 500;
double lr = 1.;
double lrBASE = 0.5;
double tol = 1e-8;
int cnt_retry = 0;

while ( arma::accu(arma::abs(anew - aold)) > tol || i < MINI_ITER){
    if (isnan(anew[0]) || isnan(anew[1]) ) break;
    i = i + 1;
    m = anew(0);
    a = anew(1);
    aold = anew;
    y = m * x / a;
    tanhy = arma::tanh(y);
    coshy2 = 1.0 / arma::pow(arma::cosh(y), 2.0);
    derm = -n * m / a + arma::accu(x % tanhy) / a;
    ders = -n / 2.0 / a + (sx2 + n * std::pow(m, 2.0)) / 2.0 / std::pow(a, 2.0) - arma::accu(y % tanhy) / a;
    derm2 = -n / a + arma::accu(x2 % coshy2) / std::pow(a, 2.0);
    ders2 = n / 2.0 / std::pow(a, 2.0) - (sx2 + n * std::pow(m, 2.0)) / std::pow(a, 3.0) + 2.0 * arma::accu(y % tanhy) / std::pow(a, 2.0) + arma::accu(arma::pow(y, 2.0) % coshy2) / std::pow(a, 2.0);
    derms = n * m / std::pow(a, 2.0) - arma::accu(x % tanhy) / std::pow(a, 2.0) - arma::accu(y % x % coshy2) / std::pow(a, 2.0);
    tmp(0) = ders2 * derm - derms * ders;
    tmp(1) = -derms * derm + derm2 * ders;
    tmp = tmp / (derm2 * ders2 - std::pow(derms, 2.0));
    anew = aold - lrBASE*lr* tmp;
    if (i>MAX_ITER__) break;
    if (i % 10 == 0)lr = (1. - (double) i /(double) MAX_ITER__)/2. + 0.5;
    if (isnan(anew[0]) || isnan(anew[1])) {
        //std::cerr<<"Try new Learning rate"<<std::endl;
        lrBASE = lrBASE /5.;
        anew[0] = aold[0];
        anew[1] = aold[1];
        if (++cnt_retry>5) break;
    }
}

m = anew(0);
a = anew(1);
double loglik = n / 2.0 * std::log(2.0 / arma::datum::pi / a) - 0.5 * n * std::pow(m, 2.0) / a - 0.5 * arma::accu(x2) / a + arma::accu(arma::log(arma::cosh(y)));
double loglik_null = (double) n*0.5 *log(2. * (double)n / 3.1415926535897932 /sx2) - (double) n / 2.;
double lrt;

lrt = 2. * (loglik - loglik_null);
//if (lrt<0.) std::cerr << "LRT-error"<<std::endl;

std::vector<double> result(3);
result[2] = 1 - chisq_mix_P(lrt, 0.5);
result[0] = anew(0);
result[1] = anew(1);
return result;
}
#endif