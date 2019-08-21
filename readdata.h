#include "essentials.h"

bool isFloat( string myString );
vector< vector<double> > column_filter( vector< vector<double> > &DATA, int col_index, double value );

//
// this is a generic function that reads in a textfile into a vector< vector<double> >
//
void readdata( string filename, vector< vector<double> > &DATA, int numcol, bool header  );

