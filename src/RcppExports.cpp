// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// annealing
double annealing(NumericMatrix Corr, NumericVector S);
RcppExport SEXP evolqg_annealing(SEXP CorrSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Corr(CorrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type S(SSEXP);
    __result = Rcpp::wrap(annealing(Corr, S));
    return __result;
END_RCPP
}
// createRandomMatrix
Eigen::MatrixXd createRandomMatrix(int dimension, float ke);
RcppExport SEXP evolqg_createRandomMatrix(SEXP dimensionSEXP, SEXP keSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type dimension(dimensionSEXP);
    Rcpp::traits::input_parameter< float >::type ke(keSEXP);
    __result = Rcpp::wrap(createRandomMatrix(dimension, ke));
    return __result;
END_RCPP
}
