// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// LOD_fit
List LOD_fit(arma::vec y_data, arma::mat x_data, arma::vec mean_x_preds, arma::vec beta, double sigma_2_y, arma::mat sigma_x_preds, int no_of_samples, double threshold, int max_iterations, arma::mat LOD_u_l);
RcppExport SEXP _lodr_LOD_fit(SEXP y_dataSEXP, SEXP x_dataSEXP, SEXP mean_x_predsSEXP, SEXP betaSEXP, SEXP sigma_2_ySEXP, SEXP sigma_x_predsSEXP, SEXP no_of_samplesSEXP, SEXP thresholdSEXP, SEXP max_iterationsSEXP, SEXP LOD_u_lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y_data(y_dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_data(x_dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mean_x_preds(mean_x_predsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_2_y(sigma_2_ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_x_preds(sigma_x_predsSEXP);
    Rcpp::traits::input_parameter< int >::type no_of_samples(no_of_samplesSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type LOD_u_l(LOD_u_lSEXP);
    rcpp_result_gen = Rcpp::wrap(LOD_fit(y_data, x_data, mean_x_preds, beta, sigma_2_y, sigma_x_preds, no_of_samples, threshold, max_iterations, LOD_u_l));
    return rcpp_result_gen;
END_RCPP
}
// LOD_bootstrap_fit
List LOD_bootstrap_fit(int num_of_boots, arma::vec y_data, arma::mat x_data, int no_of_samples, double threshold, int max_iterations, arma::mat LOD_u_l);
RcppExport SEXP _lodr_LOD_bootstrap_fit(SEXP num_of_bootsSEXP, SEXP y_dataSEXP, SEXP x_dataSEXP, SEXP no_of_samplesSEXP, SEXP thresholdSEXP, SEXP max_iterationsSEXP, SEXP LOD_u_lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type num_of_boots(num_of_bootsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y_data(y_dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_data(x_dataSEXP);
    Rcpp::traits::input_parameter< int >::type no_of_samples(no_of_samplesSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type LOD_u_l(LOD_u_lSEXP);
    rcpp_result_gen = Rcpp::wrap(LOD_bootstrap_fit(num_of_boots, y_data, x_data, no_of_samples, threshold, max_iterations, LOD_u_l));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_lodr_LOD_fit", (DL_FUNC) &_lodr_LOD_fit, 10},
    {"_lodr_LOD_bootstrap_fit", (DL_FUNC) &_lodr_LOD_bootstrap_fit, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_lodr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
