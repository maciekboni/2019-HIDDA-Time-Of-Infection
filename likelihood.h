#include "essentials.h"
#include "prms.h"


// MUST HAVE these two arguments: gsl_vector (pararms over which we're optimizing) and a void* that points to the other parameters 
double minus_loglikelihood_Senterica_re_topt(const gsl_vector *x, void *params);



// integrates the product of these two normals from -inf to +inf
// double convolve_two_normals(double mu1, double sigma1, double mu2, double sigma2);
double convolve_two_normals_by_trapezoid(double mu1, double sigma1, double mu2, double sigma2, double numsd=7.0, double dx=0.01);
double convolve_two_normals_by_trapezoid_finitelowerbound(double mu1, double sigma1, double mu2, double sigma2, double lb, double numsd=7.0, double dx=0.01);
double convolve_normalpdf_with_normalcdf_by_trapezoid(double mu1, double sigma1, double critval, double mu2, double sigma2, double numsd=7.0, double dx=0.01);
double like_4points_randomeffect( double xi[4], double pred[4], double sigma, double limit_of_detection, double re_mean, double re_sigma, double numsd=7.0, double db=0.01 );
