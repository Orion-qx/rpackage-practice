# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

sparseRegressionCPP <- function(cMatrix, constant, bvec, xyvec, penalty, p, tol) {
    .Call('_hw4_sparseRegressionCPP', PACKAGE = 'hw4', cMatrix, constant, bvec, xyvec, penalty, p, tol)
}

