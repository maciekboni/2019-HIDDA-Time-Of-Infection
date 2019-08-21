#ifndef PRMS
#define PRMS

#include <vector>
#include <math.h>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h> // random number generators from Gnu Scientific Library
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>  // this includes the gsl_ran_gaussian function


using namespace std;

enum parameter_index { i_slope, i_sigma, i_mean_starttiter_H3, i_sigma_starttiter_H3, i_a, i_topt, num_params }; 
typedef enum parameter_index prm_index;

extern gsl_rng *G_RNG;	

class prms
{   
public:    
    explicit prms( void );
    ~prms();         	// destructor

    //static string parameter_name[3] = {"slope","sigma","mean_starttiter"};
    
    int dim;		     // the number of parameters
    int dim_search;	     // the number of parameters we are searching over

    vector<double> v;			// this holds all the parameters; they are indexed by the enums above
    vector<double> v_initial_values;	// this holds initial values for all the parameters (when doing a maxlike search)

    vector<double> v_stepsizes;                 // for initial steps in Nelder-Mead search
    vector<double> v_sigmas;
    vector<double> v_leftbound;
    vector<double> v_rightbound;
    
    vector<double> v_stepsize_for_profile;      // this is the stepsize to move forward by when profiling this parameter

    

    // ### --- MOST IMPORTANT PART OF THIS CLASS TO UNDERSTAND
    //
    //         The vector below stores an ascending sequence of unique integers, e.g. something like 2, 4, 7, 9 that tells
    //         you which parameters -- according to the 'enum parameter_index' or 'prm_index' -- are being searched over
    //
    //         Hence, the vector below allows you to perform manipulations on the search parameters only, without affecting 
    //         all of the other parameters in the class that you would like to stay fixed (permanently or temporarily)
    //
    vector<int> v_indices_search_params;
    //
    //
    // ### --- MOST IMPORTANT PART OF THIS CLASS TO UNDERSTAND


    
    vector< vector<double> > vv_stored_search_params;	// each row here is a vector of search params that has been
							// stored for some reason, probably because they are optimal
    
    bool include_penalty; // include the penalty calculation inside the likelihood function; if you set this 
			  // to false, then the likelihood function returns the pure likelihood which has nice
			  // statistical properties that can be used for profiling, ci, and like-ratio tests

    // these are pointers to GSL ODE structures that integrate the ODE system
    // gsl_odeiv_step* 	    os;
    // gsl_odeiv_control*   oc;
    // gsl_odeiv_evolve*    oe;

    
    
    // member functions
    double get_val_search_param( int index );	//  
    double get_iv_search_param( int index );	// iv = inital value
    
    void multiply_sigmas_by( double f );
    
    void absorb( const gsl_vector* x );
    void absorb_stored_params( int r );
    
    void copy_to_iv_vector( void );
    void copy_from_iv_vector( void );
    void copy_all_params_from( prms* p );
    
};

#endif // PRMS
