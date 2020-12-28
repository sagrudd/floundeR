// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// sequence_quality
double sequence_quality(std::string qstring);
RcppExport SEXP _floundeR_sequence_quality(SEXP qstringSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type qstring(qstringSEXP);
    rcpp_result_gen = Rcpp::wrap(sequence_quality(qstring));
    return rcpp_result_gen;
END_RCPP
}
// fishy_fastq
int fishy_fastq(std::string fastq);
RcppExport SEXP _floundeR_fishy_fastq(SEXP fastqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fastq(fastqSEXP);
    rcpp_result_gen = Rcpp::wrap(fishy_fastq(fastq));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_floundeR_sequence_quality", (DL_FUNC) &_floundeR_sequence_quality, 1},
    {"_floundeR_fishy_fastq", (DL_FUNC) &_floundeR_fishy_fastq, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_floundeR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
