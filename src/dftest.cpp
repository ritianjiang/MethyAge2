#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
RcppExport SEXP dftest(SEXP df, SEXP vname)
{
  DataFrame DF = as<DataFrame>(df);
  std::string var = as<std::string>(vname);
  NumericVector v = DF[var];
  return(wrap(v));
}
