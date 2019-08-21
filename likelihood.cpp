#include "essentials.h"
#include "prms.h"
#include "likelihood.h"


// GLOBAL STRUCTURE THAT HOLDS THE DATA
extern vector< vector<double> > G_DATA;

extern double G_CLO_DX;
extern double G_EXTRA_POINT;
extern bool G_B_INCLUDE_EXTRA_POINT;

// double convolve_two_normals(double mu1, double sigma1, double mu2, double sigma2)
// {
//     //gsl_integration_fixed_workspace *w;
//     //const gsl_integration_fixed_type *T = gsl_integration_fixed_hermite;
//     
//     return 0.0;   
// }

double convolve_two_normals_by_trapezoid(double mu1, double sigma1, double mu2, double sigma2, double numsd, double dx)
{
    double smaller_mean = (mu1<mu2) ? mu1 : mu2;
    double larger_mean  = (mu1<mu2) ? mu2 : mu1;
    double larger_sigma = (sigma1<sigma2) ? sigma2 : sigma1;
    
    double start_x= smaller_mean - numsd*larger_sigma;
    double end_x  = larger_mean  + numsd*larger_sigma;
    double x = start_x;

    double pdfval1 = ( 1.0/(sigma1*SQRT2PI) ) * exp( (mu1-x)*(x-mu1)/(2.0*sigma1*sigma1) );
    double pdfval2 = ( 1.0/(sigma2*SQRT2PI) ) * exp( (mu2-x)*(x-mu2)/(2.0*sigma2*sigma2) );
    
    double ya = pdfval1*pdfval2;
    double yb;
    
    double sum=0.0;
    
    for( x=start_x+dx; x<=end_x; x+=dx )
    {
        pdfval1 = ( 1.0/(sigma1*SQRT2PI) ) * exp( (mu1-x)*(x-mu1)/(2.0*sigma1*sigma1) );
        pdfval2 = ( 1.0/(sigma2*SQRT2PI) ) * exp( (mu2-x)*(x-mu2)/(2.0*sigma2*sigma2) );
        yb = pdfval1*pdfval2;
        sum += (ya+yb)*dx; // multiply by 0.5 at the end to make it more efficient
        ya=yb;
    }
    
    return 0.5*sum;   
}

double convolve_two_normals_by_trapezoid_finitelowerbound(double mu1, double sigma1, double mu2, double sigma2, double lb, double numsd, double dx)
{
    double smaller_mean = (mu1<mu2) ? mu1 : mu2;
    double larger_mean  = (mu1<mu2) ? mu2 : mu1;
    double larger_sigma = (sigma1<sigma2) ? sigma2 : sigma1;
    
    double start_x= lb;
    double end_x  = larger_mean  + numsd*larger_sigma;
    double x = start_x;

    double pdfval1 = ( 1.0/(sigma1*SQRT2PI) ) * exp( (mu1-x)*(x-mu1)/(2.0*sigma1*sigma1) );
    double pdfval2 = ( 1.0/(sigma2*SQRT2PI) ) * exp( (mu2-x)*(x-mu2)/(2.0*sigma2*sigma2) );
    
    double ya = pdfval1*pdfval2;
    double yb;
    
    double sum=0.0;
    
    for( x=start_x+dx; x<=end_x; x+=dx )
    {
        pdfval1 = ( 1.0/(sigma1*SQRT2PI) ) * exp( (mu1-x)*(x-mu1)/(2.0*sigma1*sigma1) );
        pdfval2 = ( 1.0/(sigma2*SQRT2PI) ) * exp( (mu2-x)*(x-mu2)/(2.0*sigma2*sigma2) );
        yb = pdfval1*pdfval2;
        sum += (ya+yb)*dx; // multiply by 0.5 at the end to make it more efficient
        ya=yb;
    }
    
    return 0.5*sum;   
}


double convolve_normalpdf_with_normalcdf_by_trapezoid(double mu1, double sigma1, double critval, double mu2, double sigma2, double numsd, double dx)
{
//     double smaller_mean = (mu1<mu2) ? mu1 : mu2;
//     double larger_mean  = (mu1<mu2) ? mu2 : mu1;
//     double larger_sigma = (sigma1<sigma2) ? sigma2 : sigma1;
    
    double start_x= mu2 - numsd*sigma2;
    double end_x  = mu2 + numsd*sigma2;
    double x = start_x;

    double cdfval1 = 1.0 - gsl_cdf_gaussian_P( critval - (x+mu1), sigma1 );
    double pdfval2 = ( 1.0/(sigma2*SQRT2PI) ) * exp( (mu2-x)*(x-mu2)/(2.0*sigma2*sigma2) );
    
    double ya = cdfval1*pdfval2;
    double yb;
    
    double sum=0.0;
    
    for( x=start_x+dx; x<=end_x; x+=dx )
    {
        cdfval1 = 1.0 - gsl_cdf_gaussian_P( critval - (x+mu1), sigma1 );
        pdfval2 = ( 1.0/(sigma2*SQRT2PI) ) * exp( (mu2-x)*(x-mu2)/(2.0*sigma2*sigma2) );
        yb = cdfval1*pdfval2;
        sum += (ya+yb)*dx; // multiply by 0.5 at the end to make it more efficient
        ya=yb;
    }
    
    return 0.5*sum;   
}



// convolve four normals or cdfs (each function will be either a normal pdf *or* or a normal cdf) and integrate against a normal pdfval
// this is likelihood, in a random effects framework, of four data points
//
//
// xi are the four titer data points
//
// pred are the four predicted values of the data points WITHOUT THE INTERCEPT ... i.e. assuming that the intercept is zero
//
// the intercept "b" will be integrated over in the function below
//
double like_4points_randomeffect( double xi[4], double pred[4], double sigma, double limit_of_detection, double re_mean, double re_sigma, double numsd, double db )
{
    // assume that xi[0] is the largest mean, and that xi[3] is the smallest mean
    // WARNING you can assume that they are monotonic in this analysis, but the search will not guarantee that xi[0] will be the largest

    double start_b = -3.0; // xi[3] - numsd*2.0; // WARNING assuming here that a sigma of 2.0 is big enough
    double end_b   =  3.0;  // xi[0] + numsd*2.0;
    double b = start_b;
    
    double ya=1.0;
    double yb=1.0;
    
    //fprintf(stderr, "\n\t integrating from %1.3f to %1.3f", start_b, end_b );
    
    // compute the first ya value at the left boundary
    for(int i=0; i<4; i++)
    {
        // multiply out the four pdf or cdf values for the four data points, depending on whether we have censoring or not
        if( xi[i] >= limit_of_detection )
        {
            // the probability that a normal with mean b + t*m is greater than limit_of_detection
            // which is the same as the prob that a normal with mean zero, is greater than (limit_of_detection - (b + t*m) )
            ya *= (1.0 - gsl_cdf_gaussian_P( limit_of_detection - (b+pred[i]), sigma ));
            printf("\n\t censoring");
        }
        else if( xi[i] > 0.0 ) // meaning that it's non-missing
        {
            // pdf value of xi[i] in a normal with mean b + pred[i]
            // which is pdf value in a normal with mean zero of (xi[i] - (b+pred[i]))
            ya *= ( 1.0/(sigma*SQRT2PI) ) * exp( (-1.0)*(xi[i] - (b+pred[i]))*(xi[i] - (b+pred[i]))/(2.0*sigma*sigma) );
            //printf("\n\t evaluating pdf");
        }
    }

    // multiply the result by the pdf of the random effects density
    ya *= (( 1.0/(re_sigma*SQRT2PI) ) * exp( (re_mean-b)*(b-re_mean)/(2.0*re_sigma*re_sigma) ));
        
        
    double sum=0.0;
        
    // now begin looping and computing all the yb's, i.e. the right-hand sides of the trapezoids, and then the trapezoid areas (ta)
    for( b=start_b+db; b<=end_b; b+=db )
    {
        
        yb=1.0;
        
        for(int i=0; i<4; i++)
        {
            // multiply out the four pdf or cdf values for the four data points, depending on whether we have censoring or not
            if( xi[i] >= limit_of_detection )
            {
                yb *= (1.0 - gsl_cdf_gaussian_P( limit_of_detection - (b+pred[i]), sigma ));
                printf("\n\t censoring");
            }
            else if( xi[i] > 0.0 ) // meaning that it's non-missing
            {
                yb *= ( 1.0/(sigma*SQRT2PI) ) * exp( (-1.0)*(xi[i] - (b+pred[i]))*(xi[i] - (b+pred[i]))/(2.0*sigma*sigma) );
            }
        }

        // multiply the result by the pdf of the random effects density
        yb *= (( 1.0/(re_sigma*SQRT2PI) ) * exp( (re_mean-b)*(b-re_mean)/(2.0*re_sigma*re_sigma) ));

        //fprintf(stderr, "\n\t\t Inside function like_4points_randomeffect, b=%1.3f , ya=%1.6e , yb=%1.6e", b, ya, yb ); fflush(stderr);
        
        sum += (ya+yb)*db; // multiply by 0.5 at the end to make it more efficient
        ya=yb;
        
    }
    
    return 0.5*sum;       
    
}















//
//
double minus_loglikelihood_Senterica_re_topt(const gsl_vector *x, void *params)
{
    prms* ppc = (prms*)params;
    //int n = ppc->num_10FL_patients;
    
    if( x != NULL )
    {
        ppc->absorb(x); // x is absorbed into ppc, so that the parameters contained in x 
			// overwrite any parameter values set in ppc; suppose the ppc object keeps
			// track of the parameters a0, a1, a2, .. , a9, and suppose that there
			// vector x is simply (a3, a7); then the absorb function overwrites the a3 and
			// a7 values in the ppc object with the values in the gsl_vector x
    }
    
    // loop through all the patients in G_10FL_SUBSET and compute the neg-log-like for each patient
    double sum_minus_log_like=0.0;
    
    vector< vector<double> >::iterator it;
    
    //BEGIN - LOOPING THROUGH THE SALMONELLA PATIENTS
    for( it=G_DATA.begin(); it != G_DATA.end(); it++ )
    {
        vector<double> vdata = *it; // this is the data for the current patient we're looping through
        // in the vdata above, index 0 is the patient ID, indices 1, 3, 5, 7 are the days the measurements 
        // were taken, and indices 2, 4, 6, 8 are the OD values from the ELISA
        assert( vdata[1] > 0.0 );
        
        //if( it != G_DATA.begin() ) break;
        
        // define the time intervals for the four observations
        double tt[4]; // these are the four time intervals from the first follow-up to the other follow-ups
        tt[0] = 0.0;
        tt[1] = ( vdata[3]-vdata[1] ); //  
        tt[2] = ( vdata[5]-vdata[1] ); //  
        tt[3] = ( vdata[7]-vdata[1] ); //  
        
        // shortcuts for three of the parameters
        double sd     = ppc->v[i_sigma];
        double A_H3   = ppc->v[i_mean_starttiter_H3];
        double sdA_H3 = ppc->v[i_sigma_starttiter_H3];
        //double slope
        
        // the data points .. keep -99 as the 'missing' value
        double xi[4];
//      xi[0] = vdata[1] < -50.0 ? -99.0 : log2( vdata[1]/10.0 );
//      xi[1] = vdata[2] < -50.0 ? -99.0 : log2( vdata[2]/10.0 );
//      xi[2] = vdata[3] < -50.0 ? -99.0 : log2( vdata[3]/10.0 );
//      xi[3] = vdata[4] < -50.0 ? -99.0 : log2( vdata[4]/10.0 );
        xi[0] = vdata[2] < MISSING_VALUE_BOUNDARY ? MISSING_VALUE_BOUNDARY-1 : log2( vdata[2] );
        xi[1] = vdata[4] < MISSING_VALUE_BOUNDARY ? MISSING_VALUE_BOUNDARY-1 : log2( vdata[4] );
        xi[2] = vdata[6] < MISSING_VALUE_BOUNDARY ? MISSING_VALUE_BOUNDARY-1 : log2( vdata[6] );
        xi[3] = vdata[8] < MISSING_VALUE_BOUNDARY ? MISSING_VALUE_BOUNDARY-1 : log2( vdata[8] );
        
        // the predicted values (excluding the intercept b) according to the "model" which is y = mx+b
        double pred[4];
        pred[0] = ppc->v[i_slope]*tt[0];
        pred[1] = ppc->v[i_slope]*tt[1];
        pred[2] = ppc->v[i_slope]*tt[2];
        pred[3] = ppc->v[i_slope]*tt[3];
        
        //fprintf(stderr, "\nabout to call new likelihood function %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f", xi[0], xi[1], xi[2], xi[3], pred[0], pred[1], pred[2], pred[3] ); fflush(stderr);
        
        double prob_of_data_integrated_over_b = like_4points_randomeffect( xi, pred, sd, 7.0, A_H3, sdA_H3, 7.0, G_CLO_DX );
        
        //fprintf(stderr, "\n\t\treturned from new likelihood function %1.16f ", prob_of_data_integrated_over_b ); fflush(stderr);
        
        sum_minus_log_like = sum_minus_log_like - log ( prob_of_data_integrated_over_b );
    }
    //END - LOOPING THROUGH THE SALMONELLA  PATIENTS

    
    //BEGIN - NOW COMPUTE THE LIKELIHOOD OF THE ONE EXTRA POINT, GIVEN THE FOUR PARAMS AND TOPT (TIME OF PEAK TITER)
    if( G_B_INCLUDE_EXTRA_POINT )
    {
        
        double slope = ppc->v[i_slope];
        double sd    = ppc->v[i_sigma];
        double A     = ppc->v[i_mean_starttiter_H3];;
        double sdA   = ppc->v[i_sigma_starttiter_H3];
        double topt  = ppc->v[i_topt]; // a positive number indicating that this many days ago was your peak titer
        
        // for the convolution below, you have to start the integral at b=G_EXTRA_POINT, and you have to renormalize
        // the probability afterwards
        
        // take the log of the OD of the extra point
        double log_new_titer = log2( G_EXTRA_POINT );
        
        // convolve_two_normals_by_trapezoid_finitelowerbound            (double mu1, double sigma1, double mu2, double sigma2, double lb, double numsd, double dx)
        double prob2;
        prob2 = convolve_two_normals_by_trapezoid_finitelowerbound( log_new_titer - slope*topt, sd, A, sdA, log_new_titer, 7.0, G_CLO_DX);
        if(prob2==0.0) prob2 = 1e-31;
        double normalizing_factor = 1.0 - gsl_cdf_gaussian_P( log_new_titer, sdA );
        sum_minus_log_like = sum_minus_log_like - log ( prob2/normalizing_factor ); 
    }
    //END     - COMPUTING THE LIKELIHOOD OF THE ONE EXTRA POINT, GIVEN THE FOUR PARAMS AND TOPT (TIME OF PEAK TITER)
    

    return sum_minus_log_like;
}




















