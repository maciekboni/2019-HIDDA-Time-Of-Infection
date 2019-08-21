//#include <iostream>
//#include <string>
//#include <cstdlib>

#include "assert.h"
#include "prms.h"
#include "likelihood.h"


// constructor
prms::prms( void )
{
    v.insert(                v.begin(),                num_params, 0.0 );
    v_initial_values.insert( v_initial_values.begin(), num_params, 0.0 );
    v_stepsizes.insert(      v_stepsizes.begin(),      num_params, 0.0 );
    v_sigmas.insert(         v_sigmas.begin(),         num_params, 0.0 );
    v_leftbound.insert(      v_leftbound.begin(),      num_params, 0.0 );
    v_rightbound.insert(     v_rightbound.begin(),     num_params, 0.0 );
    
    assert( v.size()==num_params );
    assert( v_initial_values.size()==num_params );
    assert( v_indices_search_params.size()==0 );
    
    dim = num_params;
    dim_search=0;
    
    include_penalty = false;
    
    //os=NULL;
    //oc=NULL;
    //oe=NULL;
}



// destructor
prms::~prms()
{

    
}

void prms::multiply_sigmas_by( double f )
{
    for( int i=0; i<num_params; i++ )
    {
        v_sigmas[i] = v_sigmas[i] * f;    
    }
}

double prms::get_val_search_param( int index )
{
    assert( index >= 0 && index < dim_search );
    return v[ v_indices_search_params[index] ];
}

double prms::get_iv_search_param( int index )
{
    assert( index >= 0 && index < dim_search );
    return v_initial_values[ v_indices_search_params[index] ];
}


void prms::absorb( const gsl_vector* x )
{
    assert( x->size == dim_search );	
    assert( x->size == v_indices_search_params.size() );
    
    for( int k=0; k<dim_search; k++ )
    {
        v[ v_indices_search_params[k] ] = gsl_vector_get( x, k );
    }

    return;
}

void prms::absorb_stored_params( int r )
{
    int num_stored_sets = vv_stored_search_params.size(); // gives you the number of rows
    
    if( r >= 0 && r < num_stored_sets )
    {
        int np = vv_stored_search_params[r].size();
	assert( np==dim_search );
	
        for( int k=0; k<dim_search; k++ )
        {
            v[ v_indices_search_params[k] ] = vv_stored_search_params[r][k];
        }	
    }
    else
    {
        fprintf( stderr, "\n\tERROR INSIDE prms.cpp :: out of range in vv_stored_search_params\n" );      
    }
   
    return;
}



void prms::copy_to_iv_vector( void )
{
    for( int k=0; k<num_params; k++ )
    {
	v_initial_values[k] = v[k];    
    }
}

void prms::copy_from_iv_vector( void )
{
    for( int k=0; k<num_params; k++ )
    {
	v[k] = v_initial_values[k];    
    }
}


void prms::copy_all_params_from( prms* p )
{
    for( int k=0; k<num_params; k++ )
    {
	v[k] = p->v[k];    
    }
}


