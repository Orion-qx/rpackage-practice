#'mylm
#'
#'This function takes formula and return a fitted linear model.
#'It deals with both continuous random variable and variables can be treated as continuous random variable.
#'It can also help with parameter selection for sparse regression.
#'
#'@param formula an object of class 'formula': a symbolic description of the model to be fitted.
#'@param data an data frame which includes all x's and y.
#'@param penalty an optional double used for model selection for sparse regression. Default is 0 penalty.
#'@param tol an optional double used for approximating estimation of coefficient for sparse regression. Default is 1e-17.
#'@param useSparse an optional logical value. If TRUE, the model will use coordinate descent to approximate estimator of coefficients.
#'@param useCpp an optional logical value. If TRUE, the model will use the same method of approximation but written in C++.
#'
#'@return mylm returns an list contains estimations of coefficients, the corresponding result from t test, fitted values, and residuals.
#'
#'@examples
#'library(Rcpp)
#'library(plyr)
#'mylm(mpg~cyl+disp, mtcars) #sourceCpp("sparseRregressionCPP.cpp")
#'
## usethis namespace: start
#' @useDynLib hw4, .registration = TRUE
## usethis namespace: end
NULL
## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
#'
#'@export
#'
mylm = function(f, data, penalty=0, tol=1e-17, useSparse=FALSE, useCpp=FALSE) {
  # browser()
  if (!is.formula(f)) {
    stop("input must be a formula")
  }
  vars <- all.vars(f)# interpret input formula
  nvars <- length(vars)
  y <- as.matrix(data[vars[1]])
  X <- as.matrix(data[vars[-1]])
  tm <- terms(f)
  includeIntercept <- attr(tm,"intercept") == 1
  lmvalues <- mysimplelm.fit(X, y, includeIntercept, f, vars, tm, penalty, tol, useSparse, useCpp)
  return (lmvalues)
}

mysimplelm.fit = function(X, y, include.intercept, f, vars, tm, penalty, tol, useSparse, useCpp) {
  nrowX <- nrow(X)
  ncolX <- ncol(X)
  if (is.null(nrowX)) { # if X is not a matrix
    if (is.vector(X)) { # if X is a vector
      X <- as.matrix(X)
    } else {
      stop("X is not a matrix/vector")
    }
  }
  if (ncol(y) != 1) {
    stop("y must be a vector")
  }

  if (nrowX != length(y)) { # make sure X and y have the same row number
    stop("matrix/vector dimensions differ")
  }
  p <- ncolX
  smalldata <- nrowX < ncolX
  if (include.intercept) {
    X <- cbind(1, X)
    p <- p+1
  }
  # smalldata=TRUE
  # useCpp = FALSE
  # use normal regression method or sparse regression
  if (smalldata | useSparse | penalty != 0) {
    cMatrix <- t(X)%*%X
    const <- diag(cMatrix)
    bvec <- numeric(p)
    xyvec <- t(X)%*%y
    if(useCpp) {
      beta <- sparseRegression(cMatrix, const, bvec, xyvec, penalty, p, tol)
      # print(beta)
    } else {
      beta <- sparseRegression(cMatrix, const, bvec, xyvec, penalty, p, tol)
    }
    # print(beta)
  } else {
    beta <- solve(t(X) %*% X) %*% t(X) %*% y
  }


  # get coefficients
  betaVec <- as.vector(beta)
  variables <- vars[-1]
  if (include.intercept) {
    attr(betaVec, "names") <- c("(Intercept)", variables)
  } else {
    attr(betaVec, "names") <- c(variables)
  }
  names(betaVec) <- attr(betaVec, "names")
  if (!smalldata) {
    # get fitted values
    nobserves <- length(y)
    fittedy <- as.numeric(X%*%beta)
    idx <- 1:nobserves
    attr(fittedy, "names") <- idx

    # get residuals-related data
    res <- as.numeric(y)-fittedy
    attr(res, "names") <- idx
    rankX <- ncol(X+1) # get rank
    dfres <- nobserves - rankX # get df.residual
    ressd <- as.numeric(sqrt(t(res)%*%(res)/dfres)) # get residual standard error

    # get beta variance / df
    betaVarMatrix <- ressd^2 * solve(t(X) %*% X)
    betastd <- sqrt(diag(betaVarMatrix))
    betatvals <- betaVec/betastd
    tpvals <- 2*(1-pt(abs(betatvals), df=dfres))
    betadf <- data.frame(betaVec, betastd, betatvals, tpvals, sig.code(tpvals), row.names = attr(betaVec, "names"))
    colnames(betadf) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "")

    # construct resulting table
    resl <- list("call"=f, "coefficients"=betadf, "residuals"=res, "rank"=rankX,
                 "fitted.values"=fittedy, "terms"=terms(f),
                 "residual.sd" = ressd, "df.residual" = dfres)
  } else {
    resl <- list("call"=f, "coefficients"=betaVec)
  }
  return(resl)
}

sparseRegression <- function(cMatrix, const, bvec, xyvec, penalty, ncolX, tol) {
  maxdiff = 1
  while (maxdiff > tol) {
    for (m in 1) {
      oldbvec <- bvec
      for (i in 1:ncolX) {
        approxb <- 0
        for (j in 1:ncolX) {
          if (i != j) {
            approxb <- approxb + cMatrix[i,j]*bvec[j]
          }
        }
        est <- (xyvec[i]-approxb)/const[i]
        tempconst <- penalty/const[i]
        if (est - tempconst > 0) {
          bvec[i] <- est - tempconst
        } else if (abs(est) <= tempconst) {
          bvec[i] <- 0
        } else {
          bvec[i] <- est + tempconst
        }
        maxdiff <- max(abs(bvec-oldbvec))
      }
    }
  }
  return(bvec)
}

sig.code <- function(vals) {
  vl <- length(vals)
  otpt <- vector(length = vl)
  for (v in 1:vl) {
    if (vals[v] > 0.1) {
      otpt[v] <- "   "
    } else if (vals[v] > 0.05) {
      otpt[v] <- ".  "
    } else if (vals[v] > 0.01) {
      otpt[v] <- "*  "
    } else if (vals[v] > 0.001) {
      otpt[v] <- "** "
    } else {
      otpt[v] <- "***"
    }
  }
  return(otpt)
}

# mtcars_modified <- mtcars[1:3,]
# mylm(mpg~cyl+disp+wt+qsec, data = mtcars_modified, useSparse = TRUE)
# mylm(mpg~cyl+disp+wt+qsec, data = mtcars_modified, useSparse = TRUE, penalty = 20)
# betaTrue <- c(1.5:4.5)
# X <- matrix(rnorm(500),nrow=100, ncol=4)
# y <- X%*%betaTrue+rnorm(100)
# X <- cbind(y,X)
# colnames(X) <- c("y", "col1", "col2", "col3", "col4")
# X <- as.data.frame(X)
# lm.my2 <- mylm(y~col1+col2+col3+col4, data = X)
# lm.my2$coefficients
# lm(y~col1+col2+col3+col4, data = X)

# sourceCpp("sparseRregressionCPP.cpp")
# betaTrue <- c(1.5:4.5)
# X <- matrix(rnorm(40000),nrow=10000, ncol=4)
# y <- X%*%betaTrue+rnorm(10000)
# X <- cbind(y,X)
# colnames(X) <- c("y", "col1", "col2", "col3", "col4")
# X <- as.data.frame(X)
# lm.my2 <- mylm(y~col1+col2+col3+col4, data = X)
# lm.my2cpp <- mylm(y~col1+col2+col3+col4, data = X, useCpp = TRUE)
# lm.original2 <- lm(y~col1+col2+col3+col4, data = X)
# bench::mark(mylm(y~col1+col2+col3+col4, data = X)$coefficients$Estimate,
#             mylm(y~col1+col2+col3+col4, data = X, useCpp = TRUE)$coefficients$Estimate,
#             as.vector(lm(y~col1+col2+col3+col4, data = X)$coefficients))




