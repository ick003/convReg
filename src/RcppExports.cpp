// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// timesTwo
int timesTwo(int x);
RcppExport SEXP _convReg_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}
// indicator
double indicator(double x);
RcppExport SEXP _convReg_indicator(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(indicator(x));
    return rcpp_result_gen;
END_RCPP
}
// Indicator
NumericVector Indicator(NumericVector xx);
RcppExport SEXP _convReg_Indicator(SEXP xxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    rcpp_result_gen = Rcpp::wrap(Indicator(xx));
    return rcpp_result_gen;
END_RCPP
}
// W
NumericVector W(NumericVector lam, NumericVector nu, int sumTo);
RcppExport SEXP _convReg_W(SEXP lamSEXP, SEXP nuSEXP, SEXP sumToSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type sumTo(sumToSEXP);
    rcpp_result_gen = Rcpp::wrap(W(lam, nu, sumTo));
    return rcpp_result_gen;
END_RCPP
}
// Y
NumericVector Y(NumericVector lam, NumericVector nu, int sumTo);
RcppExport SEXP _convReg_Y(SEXP lamSEXP, SEXP nuSEXP, SEXP sumToSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type sumTo(sumToSEXP);
    rcpp_result_gen = Rcpp::wrap(Y(lam, nu, sumTo));
    return rcpp_result_gen;
END_RCPP
}
// YY
NumericVector YY(NumericVector lam, NumericVector nu, int sumTo);
RcppExport SEXP _convReg_YY(SEXP lamSEXP, SEXP nuSEXP, SEXP sumToSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type sumTo(sumToSEXP);
    rcpp_result_gen = Rcpp::wrap(YY(lam, nu, sumTo));
    return rcpp_result_gen;
END_RCPP
}
// Z
NumericVector Z(NumericVector lam, NumericVector nu, int sumTo);
RcppExport SEXP _convReg_Z(SEXP lamSEXP, SEXP nuSEXP, SEXP sumToSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type sumTo(sumToSEXP);
    rcpp_result_gen = Rcpp::wrap(Z(lam, nu, sumTo));
    return rcpp_result_gen;
END_RCPP
}
// dZIP
Rcpp::NumericVector dZIP(NumericVector y, NumericVector lam, NumericVector pi, int sumTo, bool logP);
RcppExport SEXP _convReg_dZIP(SEXP ySEXP, SEXP lamSEXP, SEXP piSEXP, SEXP sumToSEXP, SEXP logPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< int >::type sumTo(sumToSEXP);
    Rcpp::traits::input_parameter< bool >::type logP(logPSEXP);
    rcpp_result_gen = Rcpp::wrap(dZIP(y, lam, pi, sumTo, logP));
    return rcpp_result_gen;
END_RCPP
}
// dHP
Rcpp::NumericVector dHP(NumericVector y, NumericVector lam, NumericVector pi, bool logP);
RcppExport SEXP _convReg_dHP(SEXP ySEXP, SEXP lamSEXP, SEXP piSEXP, SEXP logPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< bool >::type logP(logPSEXP);
    rcpp_result_gen = Rcpp::wrap(dHP(y, lam, pi, logP));
    return rcpp_result_gen;
END_RCPP
}
// dcomp
Rcpp::NumericVector dcomp(NumericVector y, NumericVector lam, NumericVector nu, int sumTo, bool logP);
RcppExport SEXP _convReg_dcomp(SEXP ySEXP, SEXP lamSEXP, SEXP nuSEXP, SEXP sumToSEXP, SEXP logPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type sumTo(sumToSEXP);
    Rcpp::traits::input_parameter< bool >::type logP(logPSEXP);
    rcpp_result_gen = Rcpp::wrap(dcomp(y, lam, nu, sumTo, logP));
    return rcpp_result_gen;
END_RCPP
}
// dtpois
Rcpp::NumericVector dtpois(NumericVector y, NumericVector lam, bool logP);
RcppExport SEXP _convReg_dtpois(SEXP ySEXP, SEXP lamSEXP, SEXP logPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< bool >::type logP(logPSEXP);
    rcpp_result_gen = Rcpp::wrap(dtpois(y, lam, logP));
    return rcpp_result_gen;
END_RCPP
}
// dBinomGauss
Rcpp::NumericVector dBinomGauss(Rcpp::NumericVector x, Rcpp::NumericVector prob, Rcpp::NumericVector size, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, bool log);
RcppExport SEXP _convReg_dBinomGauss(SEXP xSEXP, SEXP probSEXP, SEXP sizeSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dBinomGauss(x, prob, size, mu, sigma, log));
    return rcpp_result_gen;
END_RCPP
}
// dBinomLnorm
Rcpp::NumericVector dBinomLnorm(Rcpp::NumericVector x, Rcpp::NumericVector prob, Rcpp::NumericVector size, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, bool log);
RcppExport SEXP _convReg_dBinomLnorm(SEXP xSEXP, SEXP probSEXP, SEXP sizeSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dBinomLnorm(x, prob, size, mu, sigma, log));
    return rcpp_result_gen;
END_RCPP
}
// dNbinomGauss
NumericVector dNbinomGauss(NumericVector x, NumericVector mu, NumericVector size, NumericVector mug, NumericVector sigmag, bool log);
RcppExport SEXP _convReg_dNbinomGauss(SEXP xSEXP, SEXP muSEXP, SEXP sizeSEXP, SEXP mugSEXP, SEXP sigmagSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mug(mugSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmag(sigmagSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dNbinomGauss(x, mu, size, mug, sigmag, log));
    return rcpp_result_gen;
END_RCPP
}
// dZIPGauss
NumericVector dZIPGauss(NumericVector x, NumericVector lam, NumericVector pi, NumericVector mug, NumericVector sigmag, bool log);
RcppExport SEXP _convReg_dZIPGauss(SEXP xSEXP, SEXP lamSEXP, SEXP piSEXP, SEXP mugSEXP, SEXP sigmagSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mug(mugSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmag(sigmagSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dZIPGauss(x, lam, pi, mug, sigmag, log));
    return rcpp_result_gen;
END_RCPP
}
// dHPGauss
NumericVector dHPGauss(NumericVector x, NumericVector lam, NumericVector pi, NumericVector mug, NumericVector sigmag, bool log);
RcppExport SEXP _convReg_dHPGauss(SEXP xSEXP, SEXP lamSEXP, SEXP piSEXP, SEXP mugSEXP, SEXP sigmagSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mug(mugSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmag(sigmagSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dHPGauss(x, lam, pi, mug, sigmag, log));
    return rcpp_result_gen;
END_RCPP
}
// dNbinomLnorm
NumericVector dNbinomLnorm(NumericVector x, NumericVector mu, NumericVector size, NumericVector mug, NumericVector sigmag, bool log);
RcppExport SEXP _convReg_dNbinomLnorm(SEXP xSEXP, SEXP muSEXP, SEXP sizeSEXP, SEXP mugSEXP, SEXP sigmagSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mug(mugSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmag(sigmagSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dNbinomLnorm(x, mu, size, mug, sigmag, log));
    return rcpp_result_gen;
END_RCPP
}
// dCoMPoissonGauss2
NumericVector dCoMPoissonGauss2(NumericVector x, NumericVector mu, NumericVector size, NumericVector mug, NumericVector sigmag, bool log);
RcppExport SEXP _convReg_dCoMPoissonGauss2(SEXP xSEXP, SEXP muSEXP, SEXP sizeSEXP, SEXP mugSEXP, SEXP sigmagSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mug(mugSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmag(sigmagSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dCoMPoissonGauss2(x, mu, size, mug, sigmag, log));
    return rcpp_result_gen;
END_RCPP
}
// dPoisGauss
NumericVector dPoisGauss(NumericVector x, NumericVector lam, NumericVector mug, NumericVector sigmag, bool log);
RcppExport SEXP _convReg_dPoisGauss(SEXP xSEXP, SEXP lamSEXP, SEXP mugSEXP, SEXP sigmagSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mug(mugSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmag(sigmagSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dPoisGauss(x, lam, mug, sigmag, log));
    return rcpp_result_gen;
END_RCPP
}
// dPoisLnorm
Rcpp::NumericVector dPoisLnorm(Rcpp::NumericVector x, Rcpp::NumericVector lam, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, bool log);
RcppExport SEXP _convReg_dPoisLnorm(SEXP xSEXP, SEXP lamSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dPoisLnorm(x, lam, mu, sigma, log));
    return rcpp_result_gen;
END_RCPP
}
// my_dnorm
NumericVector my_dnorm(NumericVector x, NumericVector means, NumericVector sds);
RcppExport SEXP _convReg_my_dnorm(SEXP xSEXP, SEXP meansSEXP, SEXP sdsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type means(meansSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sds(sdsSEXP);
    rcpp_result_gen = Rcpp::wrap(my_dnorm(x, means, sds));
    return rcpp_result_gen;
END_RCPP
}
// dkNbinomGauss
NumericVector dkNbinomGauss(NumericVector x, NumericVector mu, NumericVector size, NumericVector mug, NumericVector sigmag, int k);
RcppExport SEXP _convReg_dkNbinomGauss(SEXP xSEXP, SEXP muSEXP, SEXP sizeSEXP, SEXP mugSEXP, SEXP sigmagSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mug(mugSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmag(sigmagSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(dkNbinomGauss(x, mu, size, mug, sigmag, k));
    return rcpp_result_gen;
END_RCPP
}
// dkNbinomLnorm
NumericVector dkNbinomLnorm(NumericVector x, NumericVector mu, NumericVector size, NumericVector mug, NumericVector sigmag, int k);
RcppExport SEXP _convReg_dkNbinomLnorm(SEXP xSEXP, SEXP muSEXP, SEXP sizeSEXP, SEXP mugSEXP, SEXP sigmagSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mug(mugSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmag(sigmagSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(dkNbinomLnorm(x, mu, size, mug, sigmag, k));
    return rcpp_result_gen;
END_RCPP
}
// dkBinomGauss
NumericMatrix dkBinomGauss(NumericVector x, NumericVector size, NumericVector prob, NumericVector mug, NumericVector sigmag);
RcppExport SEXP _convReg_dkBinomGauss(SEXP xSEXP, SEXP sizeSEXP, SEXP probSEXP, SEXP mugSEXP, SEXP sigmagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mug(mugSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmag(sigmagSEXP);
    rcpp_result_gen = Rcpp::wrap(dkBinomGauss(x, size, prob, mug, sigmag));
    return rcpp_result_gen;
END_RCPP
}
// dkPoisGauss
NumericMatrix dkPoisGauss(NumericVector x, NumericVector lam, NumericVector mug, NumericVector sigmag);
RcppExport SEXP _convReg_dkPoisGauss(SEXP xSEXP, SEXP lamSEXP, SEXP mugSEXP, SEXP sigmagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mug(mugSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmag(sigmagSEXP);
    rcpp_result_gen = Rcpp::wrap(dkPoisGauss(x, lam, mug, sigmag));
    return rcpp_result_gen;
END_RCPP
}
// dkBinomLnorm
NumericMatrix dkBinomLnorm(NumericVector x, NumericVector size, NumericVector prob, NumericVector mug, NumericVector sigmag);
RcppExport SEXP _convReg_dkBinomLnorm(SEXP xSEXP, SEXP sizeSEXP, SEXP probSEXP, SEXP mugSEXP, SEXP sigmagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mug(mugSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmag(sigmagSEXP);
    rcpp_result_gen = Rcpp::wrap(dkBinomLnorm(x, size, prob, mug, sigmag));
    return rcpp_result_gen;
END_RCPP
}
// dkPoisLnorm
NumericMatrix dkPoisLnorm(NumericVector x, NumericVector lam, NumericVector mug, NumericVector sigmag);
RcppExport SEXP _convReg_dkPoisLnorm(SEXP xSEXP, SEXP lamSEXP, SEXP mugSEXP, SEXP sigmagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mug(mugSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmag(sigmagSEXP);
    rcpp_result_gen = Rcpp::wrap(dkPoisLnorm(x, lam, mug, sigmag));
    return rcpp_result_gen;
END_RCPP
}
// funoptim
NumericVector funoptim(NumericVector thetavar, NumericVector thetafixed, CharacterVector transforms, List Xmu1, List Xmu2, List Xsigma1, List Xsigma2, NumericVector idxvar, NumericVector idxfixed);
RcppExport SEXP _convReg_funoptim(SEXP thetavarSEXP, SEXP thetafixedSEXP, SEXP transformsSEXP, SEXP Xmu1SEXP, SEXP Xmu2SEXP, SEXP Xsigma1SEXP, SEXP Xsigma2SEXP, SEXP idxvarSEXP, SEXP idxfixedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type thetavar(thetavarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetafixed(thetafixedSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type transforms(transformsSEXP);
    Rcpp::traits::input_parameter< List >::type Xmu1(Xmu1SEXP);
    Rcpp::traits::input_parameter< List >::type Xmu2(Xmu2SEXP);
    Rcpp::traits::input_parameter< List >::type Xsigma1(Xsigma1SEXP);
    Rcpp::traits::input_parameter< List >::type Xsigma2(Xsigma2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type idxvar(idxvarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type idxfixed(idxfixedSEXP);
    rcpp_result_gen = Rcpp::wrap(funoptim(thetavar, thetafixed, transforms, Xmu1, Xmu2, Xsigma1, Xsigma2, idxvar, idxfixed));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _convReg_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// f01
NumericVector f01(NumericVector x);
RcppExport SEXP _convReg_f01(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(f01(x));
    return rcpp_result_gen;
END_RCPP
}
// fpos
NumericVector fpos(NumericVector x);
RcppExport SEXP _convReg_fpos(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(fpos(x));
    return rcpp_result_gen;
END_RCPP
}
// fint
NumericVector fint(NumericVector x);
RcppExport SEXP _convReg_fint(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(fint(x));
    return rcpp_result_gen;
END_RCPP
}
// fid
NumericVector fid(NumericVector x);
RcppExport SEXP _convReg_fid(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(fid(x));
    return rcpp_result_gen;
END_RCPP
}
// ftransform
NumericVector ftransform(NumericVector x, std::string name);
RcppExport SEXP _convReg_ftransform(SEXP xSEXP, SEXP nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type name(nameSEXP);
    rcpp_result_gen = Rcpp::wrap(ftransform(x, name));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_convReg_timesTwo", (DL_FUNC) &_convReg_timesTwo, 1},
    {"_convReg_indicator", (DL_FUNC) &_convReg_indicator, 1},
    {"_convReg_Indicator", (DL_FUNC) &_convReg_Indicator, 1},
    {"_convReg_W", (DL_FUNC) &_convReg_W, 3},
    {"_convReg_Y", (DL_FUNC) &_convReg_Y, 3},
    {"_convReg_YY", (DL_FUNC) &_convReg_YY, 3},
    {"_convReg_Z", (DL_FUNC) &_convReg_Z, 3},
    {"_convReg_dZIP", (DL_FUNC) &_convReg_dZIP, 5},
    {"_convReg_dHP", (DL_FUNC) &_convReg_dHP, 4},
    {"_convReg_dcomp", (DL_FUNC) &_convReg_dcomp, 5},
    {"_convReg_dtpois", (DL_FUNC) &_convReg_dtpois, 3},
    {"_convReg_dBinomGauss", (DL_FUNC) &_convReg_dBinomGauss, 6},
    {"_convReg_dBinomLnorm", (DL_FUNC) &_convReg_dBinomLnorm, 6},
    {"_convReg_dNbinomGauss", (DL_FUNC) &_convReg_dNbinomGauss, 6},
    {"_convReg_dZIPGauss", (DL_FUNC) &_convReg_dZIPGauss, 6},
    {"_convReg_dHPGauss", (DL_FUNC) &_convReg_dHPGauss, 6},
    {"_convReg_dNbinomLnorm", (DL_FUNC) &_convReg_dNbinomLnorm, 6},
    {"_convReg_dCoMPoissonGauss2", (DL_FUNC) &_convReg_dCoMPoissonGauss2, 6},
    {"_convReg_dPoisGauss", (DL_FUNC) &_convReg_dPoisGauss, 5},
    {"_convReg_dPoisLnorm", (DL_FUNC) &_convReg_dPoisLnorm, 5},
    {"_convReg_my_dnorm", (DL_FUNC) &_convReg_my_dnorm, 3},
    {"_convReg_dkNbinomGauss", (DL_FUNC) &_convReg_dkNbinomGauss, 6},
    {"_convReg_dkNbinomLnorm", (DL_FUNC) &_convReg_dkNbinomLnorm, 6},
    {"_convReg_dkBinomGauss", (DL_FUNC) &_convReg_dkBinomGauss, 5},
    {"_convReg_dkPoisGauss", (DL_FUNC) &_convReg_dkPoisGauss, 4},
    {"_convReg_dkBinomLnorm", (DL_FUNC) &_convReg_dkBinomLnorm, 5},
    {"_convReg_dkPoisLnorm", (DL_FUNC) &_convReg_dkPoisLnorm, 4},
    {"_convReg_funoptim", (DL_FUNC) &_convReg_funoptim, 9},
    {"_convReg_rcpp_hello_world", (DL_FUNC) &_convReg_rcpp_hello_world, 0},
    {"_convReg_f01", (DL_FUNC) &_convReg_f01, 1},
    {"_convReg_fpos", (DL_FUNC) &_convReg_fpos, 1},
    {"_convReg_fint", (DL_FUNC) &_convReg_fint, 1},
    {"_convReg_fid", (DL_FUNC) &_convReg_fid, 1},
    {"_convReg_ftransform", (DL_FUNC) &_convReg_ftransform, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_convReg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}