#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector funoptim( NumericVector thetavar,
                     NumericVector thetafixed,
                     CharacterVector transforms,
                     List Xmu1,
                     List Xmu2,
                     List Xsigma1,
                     List Xsigma2,
                     NumericVector idxvar,
                     NumericVector idxfixed
                     )
{
NumericMatrix theta( thetavar.size()+thetafixed.size(), 1);

for( int i = 0; i < idxvar.size(); i++)  // copy idxvar iteself
{
   theta[idxvar[i]] = idxvar[i];
}

if(   Rf_isNull(idxfixed)  == 0 )
{
   for( int i = 0; i < thetafixed.size(); i++)
   theta[as<int>(idxfixed),1]  = thetafixed[as<int>(idxfixed)];
}

   return( 0 );
}
