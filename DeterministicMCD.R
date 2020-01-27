# --------------------------------------------------------------------
# Author: Group 1:
# Willem Houck, 429434
# Niels Janssen, 450759
# Tim de Jonge van Ellemeet, 425288
# Dani?l de Bondt 416090
# --------------------------------------------------------------------

## Use this code skeleton to implement the deterministic MCD and the plug-in
## robust regression estimator.  Please submit your implementation together
## with your report.

## IMPORTANT: Please do not change any function names and make sure that you
##            use the correct input and output as specified for each function.
##            This will simplify grading because each group's code will be
##            similarly structured.  However, feel free to add other arguments
##            to the function definitions as needed.



## Functions for initial estimators

# Input: the standardized data matrix z
# Output: the estimated covariance or correlation matrix
# Please do not use any other input or output for the initial estimators

# correlation matrix based on the hyperbolic tangent
corHT <- function(z) {
  
  return(cor(tanh(z)))
  
}

# spearman correlation matrix
corSpearman <- function(z) {
  
  return(cor(z,method='spearman'))
}

# correlation matrix based on normal scores of the ranks
corNSR <- function(z) {
  R = apply(z,2,rank,ties.method = 'average')
  
  return(cor(apply((R-1/3)/(nrow(R) + 1/3), 2,qnorm)))
  
}

# modified spatial sign covariance matrix
covMSS <- function(z) {
  s = 0
  for (i in 1:nrow(z)){
    ki = z[i, ] / sqrt(sum(z[i, ]**2))
    s = s+ ki %*% t(ki)
  }
  return(s/nrow(z))
  
  }
# covariance matrix based on first step of BACON
covBACON1 <- function(z) {
  z.norm = as.matrix(apply(z,1,function(x) sqrt(sum(x**2))))
  indices = (apply(z.norm, 2, order)[ 1:round(length(z.norm)/2, digits= 0), ])
  
  return (cor(z[c(indices),]))
  
}

# raw OGK estimator of the covariance matrix with median and Qn
rawCovOGK <- function(z) {
  # *enter your code here*
  # Assume z is already scaled properly and columns have 
  # scale equal to one
  # Hint: have a look at function covOGK() in package robustbase
  s = Qn
  # m <- function(x) {apply(x,2,median)}
  # d.inv = solve(diag(apply(z,2,s)))
  # print(d.inv)
  # print(dim(z), dim(d.inv))
  # x = z%*%d.inv
  p = ncol(z)
  u = matrix(,nrow=p,ncol=p)
  
  for (j in 1:p){
    for (k in 1:p){
      u[j,k] = .25*(s(z[, j] + z[, k])**2 - s(z[, j] - z[, k])**2)
    }
  }
  E = eigen(u, symmetric = TRUE)$vectors
  V = z%*%E
  L = diag(apply(V,2,s))**2
  # mu = as.matrix(m(V))
  epsilon = E%*%L%*%t(E)
  # mu.raw = D%*%mu
  # epsilon.raw = d.inv%*%epsilon%*%t(d.inv)
  return(epsilon)
  
  }

# rawCovOGK(iris.scale)
Cstep <- function(H, indices, z, h){
  # H should be a list containing the a subset with size h
  T.hat = list()
  S.hat = list()
  indices.all = list()
  #initialise
  T.k = apply(H, 2, mean)
  S.k = cor(H)
  T.k1 = T.k
  S.k1 = S.k
  indices.1 = indices
  while (identical(T.k, T.k1) & identical(S.k, S.k1)){
    #set previous values to next values such that when while loop is exited, the old values are returned.
    T.k = T.k1
    S.k = S.k1
    indices = indices.1
    #calculate distance for all observations, and get h smallest
    d = as.matrix(sqrt(mahalanobis(z,T.k,S.k)))
    indices.1 = (apply(d, 2, order)[ 1:h, ])
    H = z[indices.1, ]
    S.k1 = cor(H)
    T.k1 = apply(H, 2, mean)
    
  }
  # indices.all = indices
  # T.hat = T.k
  # S.hat= S.k
  return(list(T.k, S.k, indices))
}
## Main function for deterministic MCD algorithm

# Input:
# x ....... data matrix
# alpha ... proportion of observations to be used for subset size
# anything else you need

# Output
# A list with the following components:
# center ....... mean of the reweighted estimator
# cov .......... covariance matrix of the reweighted estimator
# weights ...... binary weights of the observations used for the reweighted
#                estimator (0 if outlier, 1 if used for estimation)
# raw.center ... mean of the raw estimator
# raw.cov ...... covariance matrix of the raw estimator
# best ......... indices of the observations in the best subset found by the
#                raw estimator
# any other output you want to return

covDetMCD <- function(x, alpha, ...) {
  # *enter your code here*
  #
  # Please note that the subset sizes for the MCD are not simply fractions of 
  # the number of observations in the data set, as discussed in the lectures.
  # You can use function h.alpha.n() from package robustbase to compute the 
  # subset size.
  x = as.matrix(x)
  #scale x based on coordinate wise median and Qn from robustbase package
  z =  sweep(sweep(x,2, apply(x, 2,median),FUN='-'),2,apply(x,2,Qn),FUN = '/')
  S = list(corHT(z), corSpearman(z), corNSR(z), covMSS(z), covBACON1(z), rawCovOGK(z)   )
  # eps.hat = list()
  # mu.hat = list()
  subsets = list()
  covs = list()
  #save raw estimates as (parameters x amount of initial covariances) parameters =(location, scale, indices)
  estimates.raw = list()
  h = h.alpha.n(alpha, n = nrow(z) , p = ncol(z))
  print(h)
  
  for (k in 1:length(S)){
    # for each initial cov, calculate the mean and cov for H0
    S.k = S[[k]]
    E = eigen(S.k)$vectors
    V = z%*%E
    L = diag(apply(V,2,Qn)**2)
    eps.k = E%*%L%*%t(E)
    decomp = chol(eps.k)
    mu.k = decomp %*% apply(z%*%solve(decomp), 2, median)
    # eps.hat[[k]] = eps.k
    # mu.hat[[k]] = mu.k
    
    #construct distances based on estimate k, select smallest h0 = n/2 and based on these observations
    #calculate d.k*. Then select subsets size h and input in Cstep
    d.k = as.matrix(sqrt(mahalanobis(z, mu.k,eps.k)))
    indices = (apply(d.k, 2, order)[ 1:round(nrow(z)/2,digits = 0), ])
    H0 = z[indices, ]
    d.k.star = as.matrix(sqrt(mahalanobis(z, apply(H0, 2, mean),cor(H0))))
    indices.star =  (apply(d.k.star, 2, order)[ 1:h, ])
    
    estimates.raw[[k]] = Cstep(z[indices.star, ], indices.star, z, h)
    covs[[k]] = estimates.raw[[k]][[2]]
    
  }
  #select estimator for which the determinent of the raw covariance is smallest and obtain raw estimates
  bestk= which.min(lapply(covs, det))
  T.raw = estimates.raw[[bestk]][[1]]
  S.raw = estimates.raw[[bestk]][[2]]
  indices.raw = estimates.raw[[bestk]][[3]]
  #Use raw MCD to detect outliers, and compute reweighted MCD
  quant = qchisq(.975, df=ncol(z))
  d = as.matrix((mahalanobis(z,T.raw,S.raw)))
  w = (apply(d, 1, function(x) (if (x<=quant){1} else{0})))
  z.weight = diag(w)%*%z
  z.weight = z.weight[as.logical(rowSums(z.weight != 0)), ]
  T.MCD = apply(z.weight,2,mean)
  S.MCD = cor(z.weight)
  return(list(T.MCD, S.MCD, w, T.raw, S.raw, indices.raw))
}
covDetMCD(iris[-5], 0.975)


## Function for regression based on the deterministic MCD

# Input:
# x ........ matrix of explanatory variables
# y ........ response variable
# alpha .... proportion of observations to be used for the subset size in the 
#            MCD estimator
# anything else you need

# Output
# A list with the following components:
# coefficients .... estimated regression coefficients based on the reweighted
#                   deterministic MCD covariance matrix
# fitted.values ... fitted values for all observations in the data
# residuals ....... residuals for all observations in the data
# MCD ............. R object for the deterministic MCD (entire output from
#                   function covDetMCD())
# any other output you want to return

lmDetMCD <- function(x, y, alpha, ...) {
  # *enter your code here*
}
