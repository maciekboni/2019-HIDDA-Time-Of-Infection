#include "essentials.h"
#include "readdata.h"
#include "prms.h"
#include "likelihood.h"


// GLOBAL RANDOM NUMBER GENERATOR
gsl_rng *G_RNG;		

// GLOBAL VARIABLES
int G_CLO_NUMITER = 5000; 	// number of iterations in the Nelder-Mead minimizer routine
bool G_CLO_PROFILE = false;     // determines whether a likelihood profile should be done (in addition to the MLE)
int G_CLO_NUM_IC = 5;           // number of initial conditions to try in NM-searches; i.e. number of Nelder-Mead
                                // searches that will be done

double G_CLO_DX = 0.02;                     // integration step size in the likelihood function                                

double G_EXTRA_POINT = 2.0;                 // this is the value of the OD whose 'peak titer time' will be inferred
bool   G_B_INCLUDE_EXTRA_POINT = false;     // this tells us whether we should include the likelihood of this extra point in the likelihood function


// GLOBAL STRUCTURES THAT HOLDS DATA
//vector< vector<double> > G_02FL_DATA;
vector< vector<double> > G_DATA;        // this is just the longitudinal data for Salmonella enterica IgG
//ssedata G_SSEDATA;

// GLOBAL STRUCTURE THAT HOLDS THE PARTICULAR 10FL DATA WE ARE ANALYZING NOW
//map< int,vector<double> > G_10FL_SUBSET;
// int n; // this is the number of patients in the above subset


// FUNCTION DECLARATIONS
void ParseArgs(int argc, char **argv);
vector<double> DoNelderMead( prms* ppc, bool profile );





int main(int argc, char* argv[])
{
    ParseArgs(argc,argv);
    
    // get random number seed from current time
    struct timeval pTV[1];
    gettimeofday( pTV, NULL );
    int seed = ((int)pTV->tv_sec) + ((int)pTV->tv_usec);  // this adds seconds to microseconds

    // make random number generator (RNG) the Mersenne Twister which has period 2^19937 - 1
    const gsl_rng_type *TT_RAND = gsl_rng_mt19937;
    G_RNG = gsl_rng_alloc(TT_RAND);
    gsl_rng_set( G_RNG, seed ); // seed the RNG    

    //string str1("1E"); if( isFloat(str1) ) cout << str1 << " is a float.\n"; else cout << str1 << " is not a float.\n";
    //string str2("T0455"); if( isFloat(str2) ) cout << str2 << " is a float.\n"; else cout << str2 << " is not a float.\n";
    //string str3("VN004"); if( isFloat(str3) ) cout << str3 << " is a float.\n"; else cout << str3 << " is not a float.\n";
    
    //    
    // ###  1.  READ IN THE SALMONELLA IGG DATA IN TAB-DELIMITED FORMAT
    //

    string fn("salmonella_enterica_longitudinal_IgG_subset138.20181216.tdl");
    readdata( fn, G_DATA, 9, true ); // the true means that there is a header
    
    
    

    
    
    //
    // ###  2.  ALLOCATE SPACE FOR A PARAMETERS CLASS -- ppc=pointer to a parameters class
    //
    
    prms* ppc = new prms;  assert( ppc );

    // initialize the ODE system structures
    // const gsl_odeiv_step_type* T = gsl_odeiv_step_rkf45;
    // ppc->os = gsl_odeiv_step_alloc(T, DIM);
    // ppc->oc = gsl_odeiv_control_y_new(1e-6, 0.0);
    // ppc->oe = gsl_odeiv_evolve_alloc(DIM);
    




    
    
    //
    // ###  3.  INITIALIZE THE PARAMETERS, THEIR INITIAL STEP-SIZES IN A SEARCH, AND THEIR LEFT AND RIGHT BOUNDS
    //
    
    ppc->v[ i_slope ] = -0.0027159;                 // this is the per-day decline in the log2(OD) antibody measurement
    ppc->v[ i_sigma ] = 0.3944;                     // this is the sd of the normal distribution that is the basic error function in the linear regression
    ppc->v[ i_mean_starttiter_H3 ]  = 1.2462;       // this is the mean of the start log2(OD) (i.e. the log2(OD) at day 21)  
    ppc->v[ i_sigma_starttiter_H3 ] = 0.3742;       // 
    ppc->v[ i_topt ] = 200.0;                       // this is the time (ago, in days) of peak titer for a new individual with a particular OD measurement

    ppc->v_stepsizes[ i_slope ] = 0.02;
    ppc->v_stepsizes[ i_sigma ] = 0.02;
    ppc->v_stepsizes[ i_mean_starttiter_H3 ] = 0.25;
    ppc->v_stepsizes[ i_sigma_starttiter_H3 ] = 0.2;
    ppc->v_stepsizes[ i_topt ] = 50.0;
    
    // these bounds left and right are for choosing random starting values in the search
    // and, (TODO check) the left and right endpoints of the prior distribution (a uniform distribution) in a Bayesian search
    ppc->v_leftbound[ i_slope ] = -0.004;
    ppc->v_leftbound[ i_sigma ] = 0.35;
    ppc->v_leftbound[ i_mean_starttiter_H3 ] = 0.0;
    ppc->v_leftbound[ i_sigma_starttiter_H3 ] = 0.30;
    ppc->v_leftbound[ i_topt ] = 200;
    
    ppc->v_rightbound[ i_slope ] = -0.001;
    ppc->v_rightbound[ i_sigma ] = 0.50;
    ppc->v_rightbound[ i_mean_starttiter_H3 ] = 2.0;
    ppc->v_rightbound[ i_sigma_starttiter_H3 ] = 0.45;
    ppc->v_rightbound[ i_topt ] = 600;

// TODO IMPLEMENT THESE BELOW    
//     ppc->v_stepsize_for_profile[ i_slope ] = 0.02;
//     ppc->v_stepsize_for_profile[ i_sigma ] = 0.02;
//     ppc->v_stepsize_for_profile[ i_mean_starttiter_H3 ] = 0.25;
//     ppc->v_stepsize_for_profile[ i_sigma_starttiter_H3 ] = 0.2;
//     ppc->v_stepsize_for_profile[ i_topt ] = 50.0;

    
    
    
    
    
    
    
    //
    // ###  4.  DEFINE THE PARAMETERS THAT MAKE UP THE SEARCH SPACE
    //
    
    ppc->v_indices_search_params.push_back( i_slope );
    ppc->v_indices_search_params.push_back( i_sigma );      // this is the sd of the error in the observation function
    ppc->v_indices_search_params.push_back( i_mean_starttiter_H3 );
    ppc->v_indices_search_params.push_back( i_sigma_starttiter_H3 );
    // ppc->v_indices_search_params.push_back( i_topt );    // NOTE:    when just inferring the slope and peak titer, do not search over this parameter; i.e. comment it out
    // G_B_INCLUDE_EXTRA_POINT = true;                      // NOTE:    if you uncomment the line above, uncomment this one as well so that the likelihood function
                                                            //          includes the extra titer value and the time of peak titer (topt)
    ppc->dim_search = ppc->v_indices_search_params.size();
    
    // populate x_params_init with the values in ppc->v (only the ones that correspond to parameters being searched over)
    gsl_vector* x_params_init = gsl_vector_alloc( ppc->dim_search );
    for( int k=0; k < ppc->dim_search; k++ )
    {
	gsl_vector_set( x_params_init, k, ppc->get_val_search_param(k) );
    } 
    double mll_inf_lb, mll_fin_lb;

    //ppc->alternate_likelihood = false;
    
    
    
    
    
    
    
    //
    // ###  5.  BASIC NELDER-MEAD LIKELIHOOD SEARCH WITH 'G_CLO_NUM_IC' INITIAL CONDITIONS IN THE SEARCH
    //
    
    if( !G_CLO_PROFILE )
    {
        cout << "\n  To perform likelihood profiling, use command-line option '-profile'\n\n  See 'ParseArgs' function in the main cpp file for other command-line options.\n";
    }
    cout << "\n\n  ::::::::::::::::::::::::::::::::::::::::::::::::::::\n  ::\n  ::  ANALYSIS 1 -- OBTAIN MLE FOR SLOPE OF TITER DECLINE, USING A RANDOM EFFECT FOR THE PEAK TITER  \n  ::\n  ::::::::::::::::::::::::::::::::::::::::::::::::::::\n";
    cout << "\n  " << G_DATA.size() << " patients included in this analysis.\n";
    
    fprintf(stdout,"\n  min-log-like \t converged? \t num_iter ");
    for( int k=0; k < ppc->dim_search; k++ )
    {
        fprintf(stdout,"\t param[%d] ", ppc->v_indices_search_params[k]);
    }
    fprintf(stdout,"\t \t best neg-log-like so far ");
    
    cout << "\n  ------------ \t ---------- \t -------- ";
    for( int k=0; k < ppc->dim_search; k++ ) fprintf(stdout,"\t ---------- " );
    fprintf(stdout,"\t \t ---------- ");
    
    
    double min_loglike_sofar = 100000000.0;
    double loglike_mle;

    for(int k=0; k<G_CLO_NUM_IC; k++)
    {
        
        //###    ###    ###    ###
        //
        //
        vector<double> v = DoNelderMead( ppc, false );  // profiling = false
        //
        //
        //###    ###    ###    ###
        
        bool bClose=false;
        if( fabs(v[0] - min_loglike_sofar) < 0.01 ) bClose=true;

        if( v[0] < min_loglike_sofar ) min_loglike_sofar = v[0];
        
        fprintf(stdout,"\n  %5.5f \t %8d \t %8d ", v[0], ((int)v[1]), ((int)v[2]) );
        for( int k=0; k < ppc->dim_search; k++ ) fprintf(stdout,"\t %1.7f ", v[k+3] );
        fprintf(stdout,"\t \t %5.5f ", min_loglike_sofar );
        
        if( bClose ) printf("*");
    }
    loglike_mle = min_loglike_sofar;
    
    printf("\n\n");

    
    
    
    
    
    

    //
    // ###  6.  LIKELIHOOOD PROFILING OF THE SLOPE PARAMETER
    //

    if(G_CLO_PROFILE){
    
    bool bLeftBoundDone = false;
    bool bRightBoundDone = false;
    bool bInsideInterval = false;
    double lb = -99.0;
    double rb = -99.0;
    
    // NOTE this is the line where you declare the index of the parameter to be profiled
    int i_profiled_parameter = i_slope;  
    
    //BEGIN removal of one parameter from the "v_indices_search_params" vector inside ppc
    
    /*printf("\nvalues before removal\n\n");
    for(int j=0; j<ppc->dim_search; j++ )
    {        
        printf("\t%1.7f", ppc->v[ ppc->v_indices_search_params[j] ] );   
    }*/
    
    int index_to_remove = -1;
    for( int j=0; j<ppc->dim_search; j++ )
    {
        int fi = ppc->v_indices_search_params[j];
        if( fi == i_profiled_parameter )    
        {
            index_to_remove = j;
            break;
        }
    }
    assert( index_to_remove >= 0 );
        
    ppc->v_indices_search_params.erase( ppc->v_indices_search_params.begin() + index_to_remove );
    ppc->dim_search = ppc->v_indices_search_params.size();
    
    /*printf("\nvalues after removal\n\n");
    for(int j=0; j<ppc->dim_search; j++ )
    {        
        printf("\t%1.7f", ppc->v[ ppc->v_indices_search_params[j] ] );   
    }*/
    
    //END removal of one parameter from the "v_indices_search_params" vector inside ppc
    
    
    
    
    
    
    
    
    
    
    cout << "\n\n  ::::::::::::::::::::::::::::::::::::::::::::::::::::\n  ::\n  ::  ANALYSIS 2 -- LIKELIHOOD PROFILE \n  ::\n  ::::::::::::::::::::::::::::::::::::::::::::::::::::\n\n";
    cout << "\n  param-val  \t min-log-like \t diff-loglike \t  ";
    cout <<   "\n  ---------    \t ------------  \t ------------  ";

    
    
    double prm        =  ppc->v_leftbound[  i_profiled_parameter ];
    double prm_end    =  ppc->v_rightbound[ i_profiled_parameter ];
    double prm_step   =  0.0002;

    double previous_loglike_difference = -99.0;
    double current_loglike_difference = -99.0;
    
    //BEGIN PROFILE
    while( prm <= prm_end && G_CLO_PROFILE )
    {
        ppc->v[ i_profiled_parameter ] = prm;
        
        min_loglike_sofar = 100000000.0;
        for(int k=0; k<G_CLO_NUM_IC; k++)
        {
        
            //###    ###    ###    ###
            //
            //
            vector<double> v = DoNelderMead( ppc, false);  // always put false for the second argument
            //
            //
            //###    ###    ###    ###
            
            if( v[0] < min_loglike_sofar ) 
            {
                min_loglike_sofar = v[0];
                //best_slope_sofar = v[3];
            }
            //printf("\n%1.5f \t %1.5f \t %1.5f", slope, min_loglike_sofar, v[3] ); // v[3] is the sigma estimate here
        }

        current_loglike_difference = loglike_mle-min_loglike_sofar;
        if( current_loglike_difference > -1.92 ) bInsideInterval = true;
        else bInsideInterval = false;
        
        if( bInsideInterval && !bLeftBoundDone )
        {
            lb = prm - ((current_loglike_difference+1.92)/(current_loglike_difference-previous_loglike_difference))*prm_step;
            bLeftBoundDone = true;
        }
        if( !bInsideInterval && bLeftBoundDone && !bRightBoundDone )
        {
            rb = (prm-prm_step) + ((previous_loglike_difference+1.92)/(previous_loglike_difference-current_loglike_difference))*prm_step;
            bRightBoundDone = true;
        }
        
        
        printf("\n  %1.5f \t %1.5f \t %1.5f  ", prm, min_loglike_sofar, loglike_mle-min_loglike_sofar ); 
        if( loglike_mle-min_loglike_sofar > -1.92 ) printf(" +");
        
        previous_loglike_difference = current_loglike_difference;
        prm += prm_step;
    }
    printf("\n\n  Confidence Interval is ( %1.6f , %1.6f )\n\n", lb, rb);
    }
    printf("\n\n");
    //END PROFILE
    

    // free memory
    // gsl_odeiv_evolve_free( ppc->oe );
    // gsl_odeiv_control_free( ppc->oc );
    // gsl_odeiv_step_free( ppc->os );
    

    delete ppc;
    return 0;
}






vector<double> DoNelderMead( prms* ppc, bool profile )
{

    vector<double> v;
    int k, fi;
    
    int sdim = ppc->dim_search;
    
    // 
    // ###  (A) SET UP THE GSL VECTOR "x_params_init" -- it holds the initial conditions for the param space search
    //
    gsl_vector* x_params_init = gsl_vector_alloc( sdim );
    for( k=0; k < sdim; k++ )
    {
        // fi is the full index; the index in the full parameter vector of all parameters
        fi = ppc->v_indices_search_params[k];
        
        // choose a random starting value for each parameter
        double rand_start = gsl_ran_flat(G_RNG, ppc->v_leftbound[fi], ppc->v_rightbound[fi]);
	
        gsl_vector_set( x_params_init, k, rand_start );
    } 
    
    //printf("\n\nAbout to set up minimization structures Nelder-Mead function\n\n"); fflush(stdout);
    
    
    
    //
    // ###  (B) SET UP The Various Minimization Structures
    //
    const gsl_multimin_fminimizer_type* TT = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss; // step size in the minimizer
    gsl_multimin_function minex_func;
    size_t iter = 0;
    int status;
    double size;    

     //printf("\n\nAbout to set step sizes in Nelder-Mead function\n\n"); fflush(stdout);
    
    
    
    //
    // ###  (C) SET INITIAL STEP SIZES IN THE SEARCH
    //
    // 2014-12-22 : these step sizes appear to be in *absolute* units
    ss = gsl_vector_alloc( sdim );
    gsl_vector_set_all(ss, 0.1);         

    
    for(k=0;k<sdim;k++)
    {
        gsl_vector_set( ss, k, ppc->v_stepsizes[ ppc->v_indices_search_params[k] ] );
    }

    
    // Initialize the objective function to be minimized
    minex_func.n = sdim;
    minex_func.f = minus_loglikelihood_Senterica_re_topt;
    minex_func.params = ppc;

    // allocate space for the minimizer and give it the objective function and initial search conditions
    s = gsl_multimin_fminimizer_alloc(TT, sdim);
    gsl_multimin_fminimizer_set(s, &minex_func, x_params_init, ss); 

    //
    // ###  (D) PERFORM THE MAIN NELDER-MEAD ITERATION TO GET MAXIMUM LIKELIHOOD ESTIMATES
    //

    
//     printf("\n\nAbout the begin Nelder-Mead Loop\n\n"); fflush(stdout);
    
    //BEGIN ITERATION OF THE NELDER-MEAD MINIMIZER
    int nConverged=0;
    do
    {
        iter++;
	//printf("\n\tIteration %d",(int)iter); fflush(stdout);
	if( iter%10000==0  )
        {
            printf("%5d f()=%7.5f size = %.7f", (int)iter, s->fval, size);
	    for(int i=0; i < sdim; i++) printf(" %1.2f ", gsl_vector_get(s->x, i) );
	    printf("\n"); fflush(stdout);
        }

        status = gsl_multimin_fminimizer_iterate(s);
           
        if(status) break;
     
        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, 1e-5);
     
        if (status == GSL_SUCCESS) { /*printf ("converged to minimum at\n");*/ nConverged=1;}
     
		
    } while (status == GSL_CONTINUE && iter < G_CLO_NUMITER);
    //END ITERATION OF THE NELDER-MEAD MINIMIZER       
    
    
    
    int i=0;
    map< int,vector<double> >::iterator it; // not declared inside this function
    /*for( it=G_10FL_SUBSET.begin(); it != G_10FL_SUBSET.end(); it++ )
    {
        int id = it->first;
        vector<double> v = it->second;
        double logttr = gsl_vector_get(s->x, i);
        printf("\n\t%d\t%1.1f\t%1.1f\t%1.1f\t%1.1f\t%1.1f \t     %1.2f  %1.1f  \t%1.5f \t%1.3f", id, v[0], v[1], v[2], v[3], v[4], logttr, 10*pow(2.0,logttr), gsl_vector_get(s->x, n+i_slope), gsl_vector_get(s->x, n+i_sigma) );
        i++;
    }
    printf("\n");*/
    
    // this is the neg-log-likelihood
    v.push_back( s->fval );
    
    // this is the convergence status
    v.push_back( nConverged );
    
    // this is the number of iterations
    v.push_back( (double)iter );
    
    for(int i=0; i<sdim; i++)
    {
        v.push_back( gsl_vector_get(s->x, i) );  // push back all of the gsl-vector's optimal search params into "v", which the function will return
    }

    
    
    // free the minimizer structure memory
    gsl_vector_free( ss );
    gsl_multimin_fminimizer_free(s);
    gsl_vector_free( x_params_init );	

    return v;
        
}








// parses command line arguments
void ParseArgs(int argc, char **argv)
{
    int i, start;
    start=1;

    /*if( argc<start )
    { 
        PrintUsageModes(); 
        exit(-1);
    }
        
    if( argv[1][0]=='n' && argv[1][1]=='o' && argv[1][2]=='n' && argv[1][3]=='e' && argv[1][4]==0 )
    {
        //fprintf(stderr, "\n\n\tnot printing to Outfile\n\n");
    }
    else 
    {
        Outfile = fopen( argv[1], "w" );
    }
    
    prm_intro_day = atof( argv[2] );*/
    
    string str;
    for(i=start; i<argc; i++)
    {
        str = argv[i];
             if( str == "-iter" )		G_CLO_NUMITER  		= atoi( argv[++i] );
	else if( str == "-numic" ) 		G_CLO_NUM_IC  		= atoi( argv[++i] );
        else if( str == "-dx" ) 		G_CLO_DX  		= atof( argv[++i] );
        else if( str == "-extrapoint" )	{       G_EXTRA_POINT	        = atof( argv[++i] ); G_B_INCLUDE_EXTRA_POINT=true;}
        else if( str == "-profile" ) 		G_CLO_PROFILE  		= true;
	else
        {
            fprintf(stderr, "\n\tUnknown option [%s] on command line.\n\n", argv[i]);
            exit(-1);
        }
    }

    return;
}
