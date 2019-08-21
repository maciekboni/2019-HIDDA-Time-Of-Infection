#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <stdlib.h>
#include <ctime>
#include <sys/time.h>

#include <gsl/gsl_rng.h> // random number generators from Gnu Scientific Library
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include <set>
#include <map>
#include <vector>
#include <algorithm>

#define MISSING_VALUE_BOUNDARY -98

#define NUMAC 1

// this is the square root of 2*pi
#define SQRT2PI 2.5066282746309994


//#define HAVE_INLINE 	// this will make GSL vector operations faster



#ifndef ESSENTIALS
#define ESSENTIALS

using namespace std;

//
//
// this function contains the ode system
//
// int func(double t, const double y[], double f[], void *params);


// void* jac;	// do this for C-compilation
//
// for C++ compilation we are replacing the void* declaration above with
// the inline dummy declaration below
inline int jac(double a1, const double* a2, double* a3, double* a4, void* a5)
{
    return 0;	
};


inline void write_to_file( const char* szFilename, vector< vector<double> >& vvDATA )
{
    FILE* fp = fopen( szFilename, "w" );
    int nr = vvDATA.size();	// number of rows
    int nc = vvDATA[0].size();
    
   for(int rr=0;rr<nr;rr++)
    {
	for(int cc=0;cc<nc;cc++)
	{
	    fprintf(fp, "%1.3f \t", vvDATA[rr][cc] );	
	}
	fprintf(fp, "\n");
    }

    fclose(fp);
    return;
}

inline void write_to_csvfile( const char* szHeader, const char* szFilename, vector< vector<double> >& vvDATA )
{
    FILE* fp = fopen( szFilename, "w" );
    
    fprintf(fp, "%s\n", szHeader );
    
    int nr = vvDATA.size();	// number of rows
    int nc = vvDATA[0].size();
    
   for(int rr=0;rr<nr;rr++)
    {
	for(int cc=0;cc<nc-1;cc++)
	{
	    fprintf(fp, "%1.6f, ", vvDATA[rr][cc] );	
	}
	fprintf(fp, "%1.6f \n", vvDATA[rr][nc-1] );	
	//fprintf(fp, "\n");
    }

    fclose(fp);
    return;
}



#endif // ESSENTIALS
