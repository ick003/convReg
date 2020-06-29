#include <Rcpp.h>
// [[Rcpp::(depends(Rcpp, compoisson))]]
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
//
// to compile and test in R.
// Rcpp::sourceCpp('c++/test.cpp', verbose = TRUE, rebuild = TRUE, dryRun = TRUE)
//

// [[Rcpp::export]]
int timesTwo(int x) {
   return x * 2;
}

// [[Rcpp::export]]
double indicator(double x){
   double res = 0.0;
   if(x == 0) res = 1;
   return res;
}


double factorial(int x) {
   double res = 1.0;
   for (double d = 1.0; d <= x; ++d) {
		res *= d;
	}
	return res;
}

double W(double lam, double nu, int sumTo) {
	double sum = 0.0;
	double factorial = 1.0;
	double lamPower = lam;
	for (int i = 1; i <= sumTo; ++i) {
		factorial *= i;
		sum += lamPower * log(factorial) / pow(factorial, nu);
		lamPower *= lam;
	}
	return sum;
}

double Y(double lam, double nu, int sumTo) {
	double sum = 0.0;
	double factorial = 1.0;
	double lamPower = lam;
	for (int i = 1; i <= sumTo; ++i) {
		factorial *= i;
		sum += lamPower * i / pow(factorial, nu);
		lamPower *= lam;
	}
	return sum;
}

double YY(double lam, double nu, int sumTo) {
   double sum = 0.0;
   double factorial = 1.0;
   double lamPower = lam;
   for (int i = 1; i <= sumTo; ++i) {
      factorial *= i;
      sum += lamPower * i * i / pow(factorial, nu);
      lamPower *= lam;
   }
   return sum;
}

double Z(double lam, double nu, int sumTo) {
	double sum = 1.0;
	double factorial = 1.0;
	double lamPower = lam;
	for (int i = 1; i <= sumTo; ++i) {
		factorial *= i;
		sum += lamPower / pow(factorial, nu);
		lamPower *= lam;
	}
	return sum;
}

// [[Rcpp::export]]
NumericVector Indicator(NumericVector xx) {
   int size = xx.size();
   NumericVector out(size);
   for (int i = 0; i < size; ++i) {
      out[i] = indicator(xx[i]);
   }
   return out;
}

// [[Rcpp::export]]
NumericVector W(NumericVector lam, NumericVector nu, int sumTo) {
	int size = lam.size();
	NumericVector out(size);
	for (int i = 0; i < size; ++i) {
		out[i] = W(lam[i], nu[i], sumTo);
	}
	return out;
}

// [[Rcpp::export]]
NumericVector Y(NumericVector lam, NumericVector nu, int sumTo) {
	int size = lam.size();
	NumericVector out(size);
	for (int i = 0; i < size; ++i) {
		out[i] = Y(lam[i], nu[i], sumTo);
	}
	return out;
}

// [[Rcpp::export]]
NumericVector YY(NumericVector lam, NumericVector nu, int sumTo) {
   int size = lam.size();
   NumericVector out(size);
   for (int i = 0; i < size; ++i) {
      out[i] = YY(lam[i], nu[i], sumTo);
   }
   return out;
}

// [[Rcpp::export]]
NumericVector Z(NumericVector lam, NumericVector nu, int sumTo) {
	int size = lam.size();
	NumericVector out(size);
	for (int i = 0; i < size; ++i) {
		out[i] = Z(lam[i], nu[i], sumTo);
	}
	return out;
}

double dcomp0(int y, double lam, double nu, int sumTo) {
	if (y < 0) {
		return 0.0;
	} else {
		return pow(lam, y) / (pow(factorial(y), nu) * Z(lam, nu, sumTo));
	}
}

double dzip0(int y, double lam, double pi) {
   if (y < 0) {
      return 0.0;
   } else {
      return R::dpois(y,lam, FALSE) * (1 - pi) + indicator(y) * pi;
   }
}


double dtpois0(int y, double lam) {
   if (y <= 0) {
      return 0.0;
   } else {
      return R::dpois(y,lam, FALSE) / (1-R::dpois(0,lam, FALSE));
   }
}

double dhp0(int y, double lam, double pi) {
   if (y < 0) {
      return 0.0;
   } else {
      return dtpois0(y,lam) * (1 - pi) + indicator(y) * pi;
   }
}

NumericVector dcomp1(NumericVector y, double lam, double nu, int sumTo) {
   int n = y.size();
   NumericVector yvec(n);
   for(int i=0; i<y.size();i++){
      yvec[i] = pow(lam, y[i]) / (pow(factorial(y[i]), nu) * Z(lam, nu, sumTo));
   }
		return yvec;
}

NumericVector dhp1(NumericVector y, double lam, double pi,bool logP = false) {

   int size = y.size();

   NumericVector ans = NumericVector(size);

   for (int i = 0; i < y.size(); ++i) {
      ans[i] = dhp0(y[i], lam, pi);
   }

   if (logP) {
      return log(ans);
   } else {
      return ans;
   }
}

double dcompoissongauss(double y, double mu, double size, double mug, double sigma) {
   double sum = 0.0;
	for (int i = 0; i <= 100; ++i) {
		sum += R::dnorm(y-i,mug,sigma,0) * dcomp0(i,mu,size,100);
	}
	return sum;
}

double dkcompoissongauss(double y, double mu, double size, double mug, double sigma, int k) {
   double sum = 0.0;
   sum = R::dnorm(y-k,mug,sigma,0) * dcomp0(k,mu,size,100);

	return sum;
}

void checkInputs(NumericVector lam, NumericVector nu) {

	int lamSize = lam.size();
	int nuSize = nu.size();

	for (int i = 0; i < lamSize; ++i) {
		if (lam[i] < 0) {
			throw exception("input 'lam' should be >= 0");
		}
	}

	for (int i = 0; i < nuSize; ++i) {
		if (nu[i] < 0) {
			throw exception("input 'nu' should be >= 0");
		}
	}
}




//' Zero-inflated Poisson distribution
//' @param y observations
//' @param lam Mean Poisson paramter
//' @param pi Inflation probability
//' @param sumTo integer value for calculating the distribution. Default to 100.
//' @param logP Logical. Should the value be on log scale.
//' @keywords zero inflated Poisson
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dZIP(NumericVector y, NumericVector lam, NumericVector pi,
                          int sumTo = 100, bool logP = false) {

   checkInputs(lam, pi);

   int size = max(NumericVector::create(y.size(), lam.size(), pi.size()));

   y = rep_len(y, size);
   lam = rep_len(lam, size);
   pi = rep_len(pi, size);

   NumericVector ans = NumericVector(size);

   for (int i = 0; i < size; ++i) {
      ans[i] = dzip0(y[i], lam[i], pi[i]);
   }

   if (logP) {
      return log(ans);
   } else {
      return ans;
   }
}

//' Hurdle Poisson distribution
//' @param y observations
//' @param lam Mean Poisson paramter
//' @param pi Hurdle probability
//' @param logP Logical. Should the value be on log scale.
//' @keywords hurdle Poisson
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dHP(NumericVector y, NumericVector lam, NumericVector pi,
                        bool logP = false) {


   int size = max(NumericVector::create(y.size(), lam.size(), pi.size()));

   y = rep_len(y, size);
   lam = rep_len(lam, size);
   pi = rep_len(pi, size);

   NumericVector ans = NumericVector(size);

   for (int i = 0; i < size; ++i) {
      ans[i] = dhp0(y[i], lam[i], pi[i]);
   }

   if (logP) {
      return log(ans);
   } else {
      return ans;
   }
}

//' Conway Maxwell Poisson distribution
//' @param y observations
//' @param lam Mean Poisson paramter
//' @param nu Variance parameter
//' @param sumTo integer value for calculating the distribution. Default to 50.
//' @param logP Logical. Should the value be on log scale.
//' @keywords Conway Maxwell Poisson
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dcomp(NumericVector y, NumericVector lam, NumericVector nu,
   	int sumTo = 50, bool logP = false) {

	checkInputs(lam, nu);

	int size = max(NumericVector::create(y.size(), lam.size(), nu.size()));

	y = rep_len(y, size);
	lam = rep_len(lam, size);
	nu = rep_len(nu, size);

	NumericVector ans = NumericVector(size);

	for (int i = 0; i < size; ++i) {
		ans[i] = dcomp0(y[i], lam[i], nu[i], sumTo);
	}

	if (logP) {
		return log(ans);
	} else {
		return ans;
	}
}


//' Zero-truncated Poisson distribution
//' @param y observations
//' @param lam Mean Poisson paramter
//' @param logP Logical. Should the value be on log scale.
//' @keywords zero truncated Poisson
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dtpois(NumericVector y, NumericVector lam,
                         bool logP = false) {


   int size = max(NumericVector::create(y.size(), lam.size()));

   y = rep_len(y, size);
   lam = rep_len(lam, size);

   NumericVector ans = NumericVector(size);

   for (int i = 0; i < size; ++i) {
      ans[i] = dtpois0(y[i], lam[i]);
   }

   if (logP) {
      return log(ans);
   } else {
      return ans;
   }
}



//' Binomial Gaussian convolution
//' @param x observations
//' @param prob Binomial probability
//' @param size Binomial size
//' @param mu Gaussian mean
//' @param sigma Gaussian variance
//' @param log Logical. Should the value be on log scale.
//' @keywords Binomial Gaussian convolution
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dBinomGauss(Rcpp::NumericVector x,
Rcpp::NumericVector prob,
Rcpp::NumericVector size,
Rcpp::NumericVector mu,
Rcpp::NumericVector sigma,
bool  log = false)
{

   int n = x.size();
   NumericVector res(n);
   int s = max(size);
   int yy = mu.size();

   NumericVector rdp(s);
   for( int i=0; i < s; i++ )
   { rdp[i] = i; }

   // we use Rcpp sugar expression doing wrapping
   for(int i=0; i < x.size(); i++)
   {
      if(yy == 1)
      res[i] = sum ( dnorm(x[i]-rdp,mu[0],sigma[0]) * Rcpp::dbinom(rdp, size[0], prob[0] ) ); else
      res[i] = sum ( dnorm(x[i]-rdp,mu[i],sigma[i]) * Rcpp::dbinom(rdp, size[i], prob[i] ) );
      if( res[i] < pow(10, -300) )
         res[i] = pow(10, -300);
   }


   // doing log transformation for the last
   if( log ) return(wrap(Rcpp::log(res)));
   return(wrap(res));

}

//' Binomial Lognormal convolution
//' @param x observations
//' @param prob Binomial probability
//' @param size Binomial size
//' @param mu Log Gaussian mean
//' @param sigma Log Gaussian variance
//' @param log Logical. Should the value be on log scale.
//' @keywords Binomial Log Gaussian convolution
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dBinomLnorm(Rcpp::NumericVector x,
Rcpp::NumericVector prob,
Rcpp::NumericVector size,
Rcpp::NumericVector mu,
Rcpp::NumericVector sigma,
bool  log = false)
{

   int n = x.size();
   NumericVector res(n);
   int s = max(size);
   int yy = mu.size();

   NumericVector rdp(s);
   for( int i=0; i < s; i++ )
   { rdp[i] = i; }

   // we use Rcpp sugar expression doing wrapping
   for(int i=0; i < x.size(); i++)
   {
      if(yy == 1)
      res[i] = sum ( dlnorm(x[i]-rdp,mu[0],sigma[0]) * Rcpp::dbinom(rdp, size[0], prob[0] ) ); else
      res[i] = sum ( dlnorm(x[i]-rdp,mu[i],sigma[i]) * Rcpp::dbinom(rdp, size[i], prob[i] ) );
      if( res[i] < pow(10, -300) )
         res[i] = pow(10, -300);
   }


   // doing log transformation for the last
   if( log ) return(wrap(Rcpp::log(res)));
   return(wrap(res));

}

//' Negative Binomial Gaussian convolution
//' @param x observations
//' @param mu Negative Binomial mean
//' @param size Negative Binomial size
//' @param mug Gaussian mean
//' @param sigmag Gaussian variance
//' @param log Logical. Should the value be on log scale.
//' @keywords Negative Binomial Gaussian convolution
//' @export
// [[Rcpp::export]]
NumericVector dNbinomGauss(NumericVector x,
                           NumericVector mu,
                           NumericVector size,
                           NumericVector mug,
                           NumericVector sigmag,
                           bool log = false
                          )
{
                           #include <Rcpp.h>
                           using namespace Rcpp;

                           int n  = x.size();
                           int yy = mu.size();

                           NumericVector res(n);
                           int s = 200;

                           NumericVector rdp(s);

                           // rep


                           for( int i=0; i < s; i++ )
                           { rdp[i] = i; }


                           // we use Rcpp sugar expression doing wrapping
for(int i=0; i < x.size(); i++)
{

 if(yy==1)
 res[i] = sum ( dnorm(x[i]-rdp,mug[0],sigmag[0]) * Rcpp::dnbinom_mu(rdp, size[0], mu[0]) ) ;   else
 res[i] = sum ( dnorm(x[i]-rdp,mug[i],sigmag[i]) * Rcpp::dnbinom_mu(rdp, size[i], mu[i]) ) ;

                           if( res[i] < pow(10, -300) )
                           res[i] = pow(10, -300);
                           if( res[i] > pow(10, 100 ) )
                           res[i] = pow(10,100);
} // for loop

                           if (log){
                              res = Rcpp::log(res);
                           }
                           return(wrap(res));

}

//' Zero inflated Poisson Gaussian convolution
//' @param x observations
//' @param lam ZIP mean
//' @param pi ZIP probability
//' @param mug Gaussian mean
//' @param sigmag Gaussian variance
//' @param log Logical. Should the value be on log scale.
//' @keywords ZIP Gaussian convolution
//' @export
// [[Rcpp::export]]
NumericVector dZIPGauss(NumericVector x,
                         NumericVector lam,
                         NumericVector pi,
                         NumericVector mug,
                         NumericVector sigmag,
                         bool log= false)
{
#include <Rcpp.h>
   using namespace Rcpp;

   int n  = x.size();
   int yy = lam.size();

   NumericVector res(n);
   int s = 100;

   NumericVector rdp(s);

   // rep


   for( int i=0; i < s; i++ )
   { rdp[i] = i; }


   // we use Rcpp sugar expression doing wrapping
   for(int i=0; i < x.size(); i++)
   {

      if(yy==1){
         res[i] = sum ( dnorm(x[i]-rdp,mug[0],sigmag[0]) * (dpois(rdp, lam[0])*(1-pi[0]) + pi[0] * Indicator(rdp)));
      } else {
            res[i] = sum ( dnorm(x[i]-rdp,mug[i],sigmag[i]) * (dpois(rdp, lam[i])*(1-pi[i]) + pi[i] * Indicator(rdp))) ;
               }

         if( res[i] < pow(10, -300) )
            res[i] = pow(10, -300);
         if( res[i] > pow(10, 100 ) )
            res[i] = pow(10,100);
   } // for loop

   if (log){
      res = Rcpp::log(res);
   }
   return(wrap(res));

}

//' Hurdle Poisson Gaussian convolution
//' @param x observations
//' @param lam Hurdle mean
//' @param pi Hurdle probability
//' @param mug Gaussian mean
//' @param sigmag Gaussian variance
//' @param log Logical. Should the value be on log scale.
//' @keywords ZIP Gaussian convolution
//' @export
// [[Rcpp::export]]
NumericVector dHPGauss(NumericVector x,
                        NumericVector lam,
                        NumericVector pi,
                        NumericVector mug,
                        NumericVector sigmag,
                        bool log= false)
{
#include <Rcpp.h>
   using namespace Rcpp;

   int n  = x.size();
   int yy = lam.size();

   NumericVector res(n);
   int s = 100;

   NumericVector rdp(s);

   // rep


   for( int i=0; i < s; i++ )
   { rdp[i] = i; }


   // we use Rcpp sugar expression doing wrapping
   for(int i=0; i < x.size(); i++)
   {

      if(yy==1){
         res[i] = sum ( dnorm(x[i]-rdp,mug[0],sigmag[0]) * (dhp1(rdp, lam[0],pi[0])));
      } else {
         res[i] = sum ( dnorm(x[i]-rdp,mug[i],sigmag[i]) * (dhp1(rdp, lam[i],pi[i])));
      }

      if( res[i] < pow(10, -300) )
         res[i] = pow(10, -300);
      if( res[i] > pow(10, 100 ) )
         res[i] = pow(10,100);
   } // for loop

   if (log){
      res = Rcpp::log(res);
   }
   return(wrap(res));

}


//' Negative Binomial Log Gaussian convolution
//' @param x observations
//' @param mu Negative Binomial mean
//' @param size Negative Binomial size
//' @param mug Log Gaussian mean
//' @param sigmag Log Gaussian variance
//' @param log Logical. Should the value be on log scale.
//' @keywords Negative Binomial Log Gaussian convolution
//' @export
// [[Rcpp::export]]
NumericVector dNbinomLnorm(NumericVector x,
                           NumericVector mu,
                           NumericVector size,
                           NumericVector mug,
                           NumericVector sigmag,
                           bool log = false)
{
                           #include <Rcpp.h>
                           using namespace Rcpp;

                           int n  = x.size();
                           int yy = mu.size();

                           NumericVector res(n);
                           int s = 100;

                           NumericVector rdp(s);

                           // rep


                           for( int i=0; i < s; i++ )
                           { rdp[i] = i; }


                           // we use Rcpp sugar expression doing wrapping
for(int i=0; i < x.size(); i++)
{

 if(yy==1)
 res[i] = sum ( dlnorm(x[i]-rdp,mug[0],sigmag[0]) * Rcpp::dnbinom_mu(rdp, size[0], mu[0]) ) ;   else
 res[i] = sum ( dlnorm(x[i]-rdp,mug[i],sigmag[i]) * Rcpp::dnbinom_mu(rdp, size[i], mu[i]) ) ;

                           if( res[i] < pow(10, -300) )
                           res[i] = pow(10, -300);
                           if( res[i] > pow(10, 100 ) )
                           res[i] = pow(10,100);
} // for loop

                           if (log){
                              res = Rcpp::log(res);
                           }
                           return(wrap(res));

}

//' Conway Maxwell Poisson Gaussian convolution
//' @param x observations
//' @param mu CoMP mean
//' @param size CoMP size
//' @param mug Gaussian mean
//' @param sigmag Gaussian variance
//' @param log Logical. Should the value be on log scale.
//' @keywords CoMP Gaussian convolution
//' @export
// [[Rcpp::export]]
NumericVector dCoMPoissonGauss2(NumericVector x,
                           NumericVector mu,
                           NumericVector size,
                           NumericVector mug,
                           NumericVector sigmag,
                           bool log = false
                          )
{
                           #include <Rcpp.h>
                           using namespace Rcpp;

                           int n  = x.size();
                           int yy = mu.size();

                           NumericVector res(n);
                           int s = 100;

                           NumericVector rdp(s);

                           // rep


                           for( int i=0; i < s; i++ )
                           { rdp[i] = i; }


                           // we use Rcpp sugar expression doing wrapping
for(int i=0; i < x.size(); i++)
{

 if(yy==1) {
    res[i] = sum ( dnorm(x[i]-rdp,mug[0],sigmag[0]) * dcomp1(rdp, mu[0], size[0], 100) ) ;
    }
 else{
    res[i] = sum ( dnorm(x[i]-rdp,mug[i],sigmag[i]) * dcomp1(rdp, mu[i], size[i],100) ) ;
 }
                           if( res[i] < pow(10, -300) )
                           res[i] = pow(10, -300);
                           if( res[i] > pow(10, 100 ) )
                           res[i] = pow(10,100);
} // for loop
                           if (log){
                              res = Rcpp::log(res);
                           }
                           return(wrap(res));

}

//' Poisson Gaussian convolution
//' @param x observations
//' @param lam Poisson mean
//' @param mug Gaussian mean
//' @param sigmag Gaussian variance
//' @param log Logical. Should the value be on log scale.
//' @keywords Poisson Gaussian convolution
//' @export
// [[Rcpp::export]]
NumericVector dPoisGauss(NumericVector x,
                           NumericVector lam,
                           NumericVector mug,
                           NumericVector sigmag,
                           bool log= false)
{
                           #include <Rcpp.h>
                           using namespace Rcpp;

                           int n  = x.size();
                           int yy = lam.size();

                           NumericVector res(n);
                           int s = 100;

                           NumericVector rdp(s);

                           // rep


                           for( int i=0; i < s; i++ )
                           { rdp[i] = i; }


                           // we use Rcpp sugar expression doing wrapping
for(int i=0; i < x.size(); i++)
{

 if(yy==1)
 res[i] = sum ( dnorm(x[i]-rdp,mug[0],sigmag[0]) * dpois(rdp, lam[0]) ) ;   else
 res[i] = sum ( dnorm(x[i]-rdp,mug[i],sigmag[i]) * dpois(rdp, lam[i]) ) ;

                           if( res[i] < pow(10, -300) )
                           res[i] = pow(10, -300);
                           if( res[i] > pow(10, 100 ) )
                           res[i] = pow(10,100);
} // for loop

                           if (log){
                              res = Rcpp::log(res);
                           }
                           return(wrap(res));

}

//' Poisson Log Gaussian convolution
//' @param x observations
//' @param lam Poisson mean
//' @param mu Log Gaussian mean
//' @param sigma Log Gaussian variance
//' @param log Logical. Should the value be on log scale.
//' @keywords Poisson Log Gaussian convolution
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dPoisLnorm(Rcpp::NumericVector x,
                              Rcpp::NumericVector lam,
                              Rcpp::NumericVector mu,
                              Rcpp::NumericVector sigma,
                              bool  log= false)
{

   int n = x.size();
   NumericVector res(n);
   int s = 100;
   int yy = lam.size();

   NumericVector rdp(s);
   for( int i=0; i < s; i++ )
   { rdp[i] = i; }

   // we use Rcpp sugar expression doing wrapping
   for(int i=0; i < x.size(); i++)
   {
      if(yy == 1)
      res[i] = sum ( dlnorm(x[i]-rdp,mu[0],sigma[0]) * Rcpp::dpois(rdp, lam[0]) ); else
      res[i] = sum ( dlnorm(x[i]-rdp,mu[i],sigma[i]) * Rcpp::dpois(rdp, lam[i]) );
      if( res[i] < pow(10, -300) )
         res[i] = pow(10, -300);
   }


   // doing log transformation for the last
   if( log ) return(wrap(Rcpp::log(res)));
   return(wrap(res));

}

// [[Rcpp::export]]
NumericVector my_dnorm( NumericVector x, NumericVector means, NumericVector sds){
    int n = x.size() ;
    NumericVector res(n) ;
    for( int i=0; i<n; i++) res[i] = R::dnorm( x[i], means[i], sds[i],0) ;
    return res ;
}

// [[Rcpp::export]]
NumericVector dkNbinomGauss(NumericVector x,
                           NumericVector mu,
                           NumericVector size,
                           NumericVector mug,
                           NumericVector sigmag,
                           int k)
{
                           #include <Rcpp.h>
                           using namespace Rcpp;

                           int n  = x.size();
                           int yy = mu.size();

                           NumericVector res(n);
                           int s = 100;

                           NumericVector rdp(s);

                           // rep


                           for( int i=0; i < s; i++ )
                           { rdp[i] = i; }


                           // we use Rcpp sugar expression doing wrapping

for(int i=0; i < x.size(); i++)
{

 if(yy==1)
 res[i] =  (dnorm(x[i]-rdp,mug[0],sigmag[0]) * Rcpp::dnbinom_mu(rdp, size[0], mu[0] ))[k];   else
 res[i] =  (dnorm(x[i]-rdp,mug[i],sigmag[i]) * Rcpp::dnbinom_mu(rdp, size[i], mu[i] ))[k];

                           if( res[i] < pow(10, -300) )
                           res[i] = pow(10, -300);
                           if( res[i] > pow(10, 100) )
                           res[i] = pow(10, 100);
} // for loop

                           return(wrap(res));

}

// [[Rcpp::export]]
NumericVector dkNbinomLnorm(NumericVector x,
                           NumericVector mu,
                           NumericVector size,
                           NumericVector mug,
                           NumericVector sigmag,
                           int k)
{
                           #include <Rcpp.h>
                           using namespace Rcpp;

                           int n  = x.size();
                           int yy = mu.size();

                           NumericVector res(n);
                           int s = 100;

                           NumericVector rdp(s);

                           // rep


                           for( int i=0; i < s; i++ )
                           { rdp[i] = i; }


                           // we use Rcpp sugar expression doing wrapping

for(int i=0; i < x.size(); i++)
{

 if(yy==1)
 res[i] =  (dlnorm(x[i]-rdp,mug[0],sigmag[0]) * Rcpp::dnbinom_mu(rdp, size[0], mu[0] ))[k];   else
 res[i] =  (dlnorm(x[i]-rdp,mug[i],sigmag[i]) * Rcpp::dnbinom_mu(rdp, size[i], mu[i] ))[k];

                           if( res[i] < pow(10, -300) )
                           res[i] = pow(10, -300);
                           if( res[i] > pow(10, 100) )
                           res[i] = pow(10, 100);
} // for loop

                           return(wrap(res));

}


// [[Rcpp::export]]
NumericMatrix dkBinomGauss(NumericVector x,
                           NumericVector size,
                           NumericVector prob,
                           NumericVector mug,
                           NumericVector sigmag
                           )
{
                           #include <Rcpp.h>
                           using namespace Rcpp;

                           int n  = x.size();
                           int yy = size.size();


                           int s = max(size)+1;
                           NumericMatrix res(n,s);

                           NumericVector rdp(s);

                           // rep


                           for( int i=0; i < s; i++ )
                           { rdp[i] = i; }


                           // we use Rcpp sugar expression doing wrapping
for(int i=0; i < x.size(); i++)
{
   for(int j=0; j < s; j++){
      if(yy==1)
      res(i,j) = R::dnorm(x[i]-rdp[j],mug[0],sigmag[0],0) * R::dbinom(rdp[j], size[0], prob[0],0);   else
      res(i,j) = R::dnorm(x[i]-rdp[j],mug[i],sigmag[i],0) * R::dbinom(rdp[j], size[i], prob[i],0);

                           if( res(i,j) < pow(10, -300) )
                           res(i,j) = pow(10, -300);
                           if( res(i,j) > pow(10, 100) )
                           res(i,j) = pow(10, 100);
   }


} // for loop

                           return(wrap(res));

}

// [[Rcpp::export]]
NumericMatrix dkPoisGauss(NumericVector x,
                           NumericVector lam,
                           NumericVector mug,
                           NumericVector sigmag
                           )
{
                           #include <Rcpp.h>
                           using namespace Rcpp;

                           int n  = x.size();
                           int yy = lam.size();


                           int s = 100;
                           NumericMatrix res(n,s);

                           NumericVector rdp(s);

                           // rep


                           for( int i=0; i < s; i++ )
                           { rdp[i] = i; }


                           // we use Rcpp sugar expression doing wrapping
for(int i=0; i < x.size(); i++)
{
   for(int j=0; j < s; j++){
      if(yy==1)
      res(i,j) = R::dnorm(x[i]-rdp[j],mug[0],sigmag[0],0) * R::dpois(rdp[j], lam[0],0);   else
      res(i,j) = R::dnorm(x[i]-rdp[j],mug[i],sigmag[i],0) * R::dpois(rdp[j], lam[i],0);

                           if( res(i,j) < pow(10, -300) )
                           res(i,j) = pow(10, -300);
                           if( res(i,j) > pow(10, 100) )
                           res(i,j) = pow(10, 100);
   }


} // for loop

                           return(wrap(res));

}

// [[Rcpp::export]]
NumericMatrix dkBinomLnorm(NumericVector x,
                           NumericVector size,
                           NumericVector prob,
                           NumericVector mug,
                           NumericVector sigmag
                           )
{
                           #include <Rcpp.h>
                           using namespace Rcpp;

                           int n  = x.size();
                           int yy = size.size();


                           int s = max(size)+1;
                           NumericMatrix res(n,s);

                           NumericVector rdp(s);

                           // rep


                           for( int i=0; i < s; i++ )
                           { rdp[i] = i; }


                           // we use Rcpp sugar expression doing wrapping
for(int i=0; i < x.size(); i++)
{
   for(int j=0; j < s; j++){
      if(yy==1)
      res(i,j) = R::dlnorm(x[i]-rdp[j],mug[0],sigmag[0],0) * R::dbinom(rdp[j], size[0], prob[0],0);   else
      res(i,j) = R::dlnorm(x[i]-rdp[j],mug[i],sigmag[i],0) * R::dbinom(rdp[j], size[i], prob[i],0);

                           if( res(i,j) < pow(10, -300) )
                           res(i,j) = pow(10, -300);
                           if( res(i,j) > pow(10, 100) )
                           res(i,j) = pow(10, 100);
   }


} // for loop

                           return(wrap(res));

}

// [[Rcpp::export]]
NumericMatrix dkPoisLnorm(NumericVector x,
                           NumericVector lam,
                           NumericVector mug,
                           NumericVector sigmag
                           )
{
                           #include <Rcpp.h>
                           using namespace Rcpp;

                           int n  = x.size();
                           int yy = lam.size();


                           int s = 100;
                           NumericMatrix res(n,s);

                           NumericVector rdp(s);

                           // rep


                           for( int i=0; i < s; i++ )
                           { rdp[i] = i; }


                           // we use Rcpp sugar expression doing wrapping
for(int i=0; i < x.size(); i++)
{
   for(int j=0; j < s; j++){
      if(yy==1)
      res(i,j) = R::dlnorm(x[i]-rdp[j],mug[0],sigmag[0],0) * R::dpois(rdp[j], lam[0], 0);   else
      res(i,j) = R::dlnorm(x[i]-rdp[j],mug[i],sigmag[i],0) * R::dpois(rdp[j], lam[i], 0);

                           if( res(i,j) < pow(10, -300) )
                           res(i,j) = pow(10, -300);
                           if( res(i,j) > pow(10, 100) )
                           res(i,j) = pow(10, 100);
   }


} // for loop

                           return(wrap(res));

}


