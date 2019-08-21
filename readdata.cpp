#include "readdata.h"

// struct tm string_to_tmstruct( string str, int format_type );
// int daysdiff( struct tm a, struct tm b );


bool isFloat( string myString ) 
{
    istringstream iss(myString);
    float f;
    iss >> noskipws >> f; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    return iss.eof() && !iss.fail(); 
}

vector< vector<double> > column_filter( vector< vector<double> > &DATA, int col_index, double value )
{
    vector< vector<double> > vv;
    
    for( int r=0; r<DATA.size(); r++ )
    {
        if( DATA[r][col_index] == value )
        {
            vv.push_back( DATA[r] );   
        }
    }
    return vv;
}



void readdata( string filename, vector< vector<double> > &DATA, int numcol, bool header  )
{
    DATA.clear();
    assert( DATA.size() == 0 );

    ifstream infile( filename.c_str(), ifstream::in );
    
    //TODO need some error handling here if the file above is not found
	
    string strHeader;
    if( header ) getline(infile, strHeader);

    int col=0;
    int row=0;
    
    int counter=0;
    
    while( true )
    {
        counter++;
	  
	string str("");
        infile >> str;
	
	//fprintf( stderr, "\n\t%d \t%d \t%d \t %1.5f", (int)str.length(), row, col, atof( str.c_str() ) ); fflush(stderr);
	
        if( str.length()==0 && infile.eof() ) break;
	
	// quick & dirty check to make sure the string is a float
	//assert( isdigit( str[0] ) || str[0]=='-' || str[0]=='.' );
	
        // if this is the first column, then it's a new row and we need
        // to allocate a new vector
	if( col==0 )
	{
	    vector<double> vd(numcol);
	    DATA.push_back( vd );
	}

	if(isFloat(str))
        {
	    DATA[row][col] = atof( str.c_str() );
        }
        else
        {
            DATA[row][col] = -101.0;
        }
	
	
	if( col==numcol-1 ) // if the previous statement just read in the last column
	{
	    col=0;
	    row++;
	}
	else
	{
	    col++;
	}
	
    }
    
    infile.close();
    return;
}



