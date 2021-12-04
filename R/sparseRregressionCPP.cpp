#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sparseRegressionCPP(NumericMatrix cMatrix,
                                  NumericVector constant, NumericVector bvec, NumericVector xyvec,
                                  double penalty, int p, double tol) {
  double maxdiff = 1;
  while (maxdiff > tol) {
    NumericVector oldbvec = bvec;
    for (int i = 0; i < p; i++) {
      double approxb = 0;
      for (int j = 0; j < p; j++) {
        if (i != j) {
          approxb = approxb + cMatrix(i,j)*bvec[j];
        }
      }
      double est = (xyvec[i]-approxb)/constant[i];
      double tempconstant = penalty/constant[i];
      if (est - tempconstant > 0) {
        bvec[i] = est - tempconstant;
      } else if (abs(est) <= tempconstant) {
        bvec[i] = 0;
      } else {
        bvec[i] = est + tempconstant;
      }
      maxdiff = max(bvec-oldbvec);
    }
  }
  return bvec;
}
