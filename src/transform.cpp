// a list of commonly used transformation function
// RCPP files
#include <Rcpp.h>
#include <math.h>
#include <stdio.h>

using namespace Rcpp;


// [[Rcpp::export]]
NumericVector f01( NumericVector x)
{
   for(int i = 0; i < x.size(); i++)
      x[i] = min(NumericVector::create(x[i],100));
   for(int i = 0; i < x.size(); i++)
      x[i] = max(NumericVector::create(x[i],-100));
   return( wrap( exp(x)/(1 + exp(x)) ) );
}


// [[Rcpp::export]]
NumericVector fpos( NumericVector x)
{
   return( wrap( exp(x) ) );
}


// [[Rcpp::export]]
NumericVector fint( NumericVector x )
{
   x = exp(x);
   for(int i = 0; i < x.size(); i++)
   x[i] = round(x[i]);
   return( x );
}

// [[Rcpp::export]]
NumericVector fid( NumericVector x)
{
   return( x);
}

// [[Rcpp::export]]
NumericVector ftransform( NumericVector x,  std::string name  )
{

   if( name.compare("01") == 0 ) return(f01(x));
   if( name.compare("pos") == 0 ) return(fpos(x));
   if( name.compare("int") == 0 ) return(fint(x));
   if( name.compare("id") == 0 ) return(x);

}   // ftransform
