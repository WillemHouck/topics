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

# install.packages("robustbase")
# install.packages("ggplot2")

library("MASS")
library("robustbase")
library("ggplot2")
library("datasets")
library("hrbrthemes")
# install.packages("mixtools")
library("mixtools")
library("rgl")
# install.packages("plotly")
library('plotly')

# Preparing data ----------------------------------------------------------
rm(list=ls())
data(iris)
load("Eredivisie28.RData")
rownames(Eredivisie28) <- 1:nrow(Eredivisie28)
Eredivisie28$MarketValue <- log(Eredivisie28$MarketValue)
options(scipen=999)

#normal plots
par(mfrow=c(1,2))
fit_age <- fitdistr(Eredivisie28$Age, densfun="normal")
hist(Eredivisie28$Age, prob=TRUE, main="", xlim = c(10,35), xlab = "Age")
curve(dnorm(x, fit_age$estimate[1], fit_age$estimate[2]), col="red", lwd=2, add=T)

fit_MV <- fitdistr(Eredivisie28$MarketValue, densfun="normal")
hist(Eredivisie28$MarketValue, prob=TRUE, main="", xlab = "log of MarketValue", xlim = c(10,20))
curve(dnorm(x, fit_MV$estimate[1], fit_MV$estimate[2]), col="red", lwd=2, add=T)

par(mfrow=c(1,1))
plot(x = Eredivisie28$Age, y = Eredivisie28$MarketValue, pch=20, cex=0.4,
     xlab="Age" , ylab="log of MarketValue", 
     xlim = c(10,35), #ylim = c(-1000000,10000000)
)
# ggplot(data = Eredivisie28, aes(x = Age, y = MarketValue)) + 
#   geom_point() +
#   theme_ipsum()


# Writing functions -------------------------------------------------------
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
    if(sum(is.na(ki)) > 0){ #if NaN, set ki to zero (or don't add to s)
      s = s
    }else{
      s = s+ ki %*% t(ki)
    }
    
  }
  return(s/nrow(z))
}
# covariance matrix based on first step of BACON
covBACON1 <- function(z) {
  z.norm = as.matrix(apply(z,1,function(x) sqrt(sum(x**2))))
  indices = (apply(z.norm, 2, order)[ 1:round(length(z.norm)/2, digits= 0), ])
  
  return (cov(z[c(indices),]))
  
}

# raw OGK estimator of the covariance matrix with median and Qn
rawCovOGK <- function(z) {
  # # *enter your code here*
  # # Assume z is already scaled properly and columns have 
  # # scale equal to one
  # # Hint: have a look at function covOGK() in package robustbase
  # s = Qn
  # # m <- function(x) {apply(x,2,median)}
  # # d = solve(diag(apply(z,2,s)))
  # # print(d.inv)
  # # print(dim(z), dim(d.inv))
  # # x = z%*%d.inv
  # p = ncol(z)
  # u = matrix(,nrow=p,ncol=p)
  # 
  # for (j in 1:p){
  #   for (k in 1:p){
  #     u[j,k] = .25*(s(z[, j] + z[, k])**2 - s(z[, j] - z[, k])**2)
  #   }
  # }
  # E = eigen(u, symmetric = TRUE)$vectors
  # V = z%*%E
  # L = diag(apply(V,2,s)**2)
  # # mu = as.matrix(m(V))
  # epsilon = E%*%L%*%t(E)
  # # mu.raw = D%*%mu
  # # epsilon.raw = d.inv%*%epsilon%*%t(d.inv)
  # return(epsilon)
  # 
  covOGK(z, sigmamu = s_Qn, n.iter = 2)$cov
}

# rawCovOGK(iris.scale)
Cstep <- function(H, indices, z, h){
  # H should be a list containing the a subset with size h
  T.hat = list()
  S.hat = list()
  indices.all = list()
  #initialise
  T.k = NA
  S.k = NA
  T.k1 = apply(H, 2, mean)
  S.k1 = cov(H)
  indices.1 = indices
  count = 0
  while (!identical(S.k, S.k1)){
    #set previous values to next values such that when while loop is exited, the old values are returned.
    T.k = T.k1
    S.k = S.k1
    # print(S.k)
    indices = indices.1
    #calculate distance for all observations, and get h smallest
    d = as.matrix(sqrt(mahalanobis(z,T.k,S.k)))
    indices.1 = (apply(d, 2, order)[ 1:h, ])
    H = z[indices.1, ]
    S.k1 = cov(H)
    T.k1 = apply(H, 2, mean)
    count = count + 1
  }
  print(S.k)
  # indices.all = indices
  # T.hat = T.k
  # S.hat= S.k
  print(paste0("#iter of Cstep: ", count))
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
  z =  sweep(sweep(x,2, apply(x, 2,median),FUN='-'),2,apply(x,2,Qn),FUN = '/') #to make the estimators equivariant
  S = list(corHT(z), corSpearman(z), corNSR(z), covMSS(z), covBACON1(z), rawCovOGK(z)   )
  # eps.hat = list()
  # mu.hat = list()
  subsets = list()
  covs = list()
  #save raw estimates as (parameters x amount of initial covariances) parameters =(location, scale, indices)
  estimates.raw = list()
  h = h.alpha.n(alpha, n = nrow(z) , p = ncol(z))
  print(paste0("#obs for constructing raw estimates (h): ",h))
  
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
    
    d.0 = as.matrix(sqrt(mahalanobis(z, mu.k,eps.k)))
    indices = (apply(d.0, 2, order)[ 1:h, ])
    # H0 = z[indices, ]
    # d.k.star = as.matrix(sqrt(mahalanobis(z, apply(H0, 2, mean),cov(H0))))
    # indices.star =  (apply(d.k.star, 2, order)[ 1:h, ])
    
    estimates.raw[[k]] = Cstep(z[indices, ], indices, z, h)
    covs[[k]] = estimates.raw[[k]][[2]]
    
  }
  
  #compute fisher consistency correction for raw variance
  alpha_fisher_raw <- h/nrow(x)
  c_fisher_raw <- alpha_fisher_raw/(pgamma(qchisq(alpha_fisher_raw, df=ncol(z))/2,
                                           shape = ncol(z)/2 + 1,
                                           scale = 1))
  
  #select estimator for which the determinent of the raw covariance is smallest and obtain raw estimates
  bestk= which.min(lapply(covs, det))
  T.raw = estimates.raw[[bestk]][[1]]
  S.raw = estimates.raw[[bestk]][[2]] * c_fisher_raw
  indices.raw = estimates.raw[[bestk]][[3]]
  
  #Use raw MCD to detect outliers, and compute reweighted MCD
  quant = qchisq(.975, df=ncol(z))
  d = as.matrix((mahalanobis(z,T.raw,S.raw)))
  w = (apply(d, 1, function(x) (if (x<=quant){1} else{0})))
  z.weight = diag(w)%*%z
  z.weight = z.weight[as.logical(rowSums(z.weight != 0)), ]
  T.MCD = apply(z.weight,2,mean)
  S.MCD = cov(z.weight)
  
  #compute fisher consistency correction for reweighted variance
  alpha_fisher_new <- sum(w)/nrow(x)
  c_fisher_new <- alpha_fisher_new/(pgamma(qchisq(alpha_fisher_new, df=ncol(z))/2,
                                           shape = ncol(z)/2 + 1,
                                           scale = 1))
  S.MCD <- S.MCD * c_fisher_new
  
  #transform from equivariant z back to x
  T.MCD.x <- T.MCD * apply(x, 2,Qn) + apply(x, 2,median)
  S.MCD.x <- sweep(sweep(S.MCD,2,t(apply(x,2,Qn)),FUN='*'),1,apply(x,2,Qn),FUN = '*')
  T.raw.x <- T.raw * apply(x, 2,Qn) + apply(x, 2,median)
  S.raw.x <- sweep(sweep(S.raw,2,t(apply(x,2,Qn)),FUN='*'),1,apply(x,2,Qn),FUN = '*')
  
  return(list(T.MCD.x, S.MCD.x, w, T.raw.x, S.raw.x, indices.raw))
}

output_E = covDetMCD(Eredivisie28, 0.975)

output_iris <- covDetMCD(iris[-5], 0.75)
output_erediv <- covDetMCD(Eredivisie28, 0.75)


output_covmcd_erediv <- covMcd(Eredivisie28, alpha = 0.75, nsamp = "deterministic")

ellipse(mu = output_erediv[[4]], sigma = output_erediv[[5]], alpha = 0.025, newplot = F, col = "red",lty  = "dashed")
ellipse(mu = output_erediv[[1]], sigma = output_erediv[[2]], alpha = 0.025, newplot = F)

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

lmDetMCD <- function(x, y, alpha) {
  
  
  #run covDetMCD and compute regression coefficients
  y_and_x <- cbind(y,x)
  MCD_estimates <- covDetMCD(x = y_and_x, alpha = alpha)
  
  sigma_xy <- MCD_estimates[[2]][1,2:ncol(y_and_x)]
  sigma_xx <- MCD_estimates[[2]][2:ncol(y_and_x),2:ncol(y_and_x)]
  lm_beta <- solve(sigma_xx)%*%sigma_xy
  
  mu_y <- MCD_estimates[[1]][1]
  mu_x <- MCD_estimates[[1]][2:ncol(y_and_x)]
  lm_alpha <- mu_y - t(mu_x)%*%lm_beta
  
  coefficients <- list(lm_alpha, lm_beta)
  
  fitted_values <- as.numeric(lm_alpha) + x %*% lm_beta  #geen transform x?
  residuals <- y - fitted_values
  
  
  return(list(coefficients,fitted_values,residuals,MCD_estimates))
}

#Plug in
output_erediv_lm <- lmDetMCD(x = as.matrix(Eredivisie28$Age), y = as.matrix(Eredivisie28$MarketValue), alpha = 0.75)
abline(coef = output_erediv_lm[[1]], col = "orange")
# lmDetMCD(x = as.matrix(iris[,1:3]), y = as.matrix(iris[,4]), alpha = 0.75)

#LTS
lts_regression <- ltsReg(x = as.matrix(Eredivisie28$Age), y = as.matrix(Eredivisie28$MarketValue), alpha = 0.75)
abline(coef = lts_regression$coefficients, col = "blue")

#OLS
ols_regression <- lm(Eredivisie28$MarketValue ~ Eredivisie28$Age)
abline(coef = ols_regression$coefficients, col = "green")

legend("topleft", legend = c("plug-in","lts","ols","reweighted","raw"), col = c("orange", "blue", "green","black","red"), lty = c(1,1,1,1,2), cex = 0.8)


# Empirical Influence Function --------------------------------------------

 calculate_EIF <- function(data){ #not a general function
  set.seed(121134)
  observation_to_change <- sample(1:nrow(data), size = 1) #value to be changed
  
  x <- seq(from = 10,to = 35,by = 1)
  y <- seq(from = 10,to = 20,by = 1)
  
  plug_in_EIF_intercept <- as.data.frame(matrix(data = NA, nrow = length(x), ncol = length(y), dimnames = list(x,y)))
  plug_in_EIF_slope <- as.data.frame(matrix(data = NA, nrow = length(x), ncol = length(y), dimnames = list(x,y)))
  
  lts_EIF_intercept <- as.data.frame(matrix(data = NA, nrow = length(x), ncol = length(y), dimnames = list(x,y)))
  lts_EIF_slope <- as.data.frame(matrix(data = NA, nrow = length(x), ncol = length(y), dimnames = list(x,y)))
  
  ols_EIF_intercept <- as.data.frame(matrix(data = NA, nrow = length(x), ncol = length(y), dimnames = list(x,y)))
  ols_EIF_slope <- as.data.frame(matrix(data = NA, nrow = length(x), ncol = length(y), dimnames = list(x,y)))
  
  
  for(i in 1:length(x)){
    adjusted_age <- Eredivisie28$Age
    adjusted_age[observation_to_change] <- x[i]
    for(j in 1:length(y)){
      
      adjusted_value <- Eredivisie28$MarketValue
      adjusted_value[observation_to_change] <- y[j]
      
      plug_in_temporary_output <- lmDetMCD(x = adjusted_age,y = adjusted_value,alpha = 0.75)[[1]]
      plug_in_EIF_intercept[i,j] <- nrow(Eredivisie28) * (plug_in_temporary_output[[1]] - output_erediv_lm[[1]][[1]])
      plug_in_EIF_slope[i,j] <- nrow(Eredivisie28) * (plug_in_temporary_output[[2]] - output_erediv_lm[[1]][[2]])
      
      lts_temporary_output <- ltsReg(x = adjusted_age, y = adjusted_value, alpha = 0.75)
      lts_EIF_intercept[i,j] <- nrow(Eredivisie28) * (lts_temporary_output$coefficients[[1]] -lts_regression$coefficients[[1]])
      lts_EIF_slope[i,j] <- nrow(Eredivisie28) * (lts_temporary_output$coefficients[[2]] - lts_regression$coefficients[[2]])
      
      ols_temporary_output <- lm(adjusted_value ~ adjusted_age)
      ols_EIF_intercept[i,j] <- nrow(Eredivisie28) * (ols_temporary_output$coefficients[[1]] - ols_regression$coefficients[[1]])
      ols_EIF_slope[i,j] <- nrow(Eredivisie28) * (ols_temporary_output$coefficients[[2]] - ols_regression$coefficients[[2]])
    }
  }
  return(list(plug_in_EIF_intercept, plug_in_EIF_slope, lts_EIF_intercept, lts_EIF_slope, ols_EIF_intercept, ols_EIF_slope))
} 

EIF_results <- calculate_EIF(Eredivisie28)



p1 <- plot_ly(x = 10:20, y = 10:35, z = as.matrix(EIF_results[[1]])) %>% 
  add_surface(contours = list(z = list(show=TRUE, usecolormap=TRUE, highlightcolor="#ff0000", project=list(z=TRUE)))) %>%
  layout(scene = list(camera=list(eye = list(x=2, y=1, z=0.3)),
                      xaxis = list(title = "log(MarketValue)"),
                      yaxis = list(title = "Age"),
                      zaxis = list(title = "EIF")))

p2 <- plot_ly(x = 10:20, y = 10:35, z = as.matrix(EIF_results[[2]])) %>% 
  add_surface(contours = list(z = list(show=TRUE, usecolormap=TRUE, highlightcolor="#ff0000", project=list(z=TRUE)))) %>%
  layout(scene = list(camera=list(eye = list(x=2, y=1, z=0.3)),
                      xaxis = list(title = "log(MarketValue)"),
                      yaxis = list(title = "Age"),
                      zaxis = list(title = "EIF")))

p3 <- plot_ly(x = 10:20, y = 10:35, z = as.matrix(EIF_results[[3]]))%>% 
  add_surface(contours = list(z = list(show=TRUE, usecolormap=TRUE, highlightcolor="#ff0000", project=list(z=TRUE)))) %>%
  layout(scene = list(camera=list(eye = list(x=2, y=1, z=0.3)),
                     xaxis = list(title = "log(MarketValue)"),
                     yaxis = list(title = "Age"),
                     zaxis = list(title = "EIF")))

p4 <- plot_ly(x = 10:20, y = 10:35, z = as.matrix(EIF_results[[4]]))%>% 
  add_surface(contours = list(z = list(show=TRUE, usecolormap=TRUE, highlightcolor="#ff0000", project=list(z=TRUE)))) %>%
  layout(scene = list(camera=list(eye = list(x=2, y=1, z=0.3)),
                        xaxis = list(title = "log(MarketValue)"),
                        yaxis = list(title = "Age"),
                        zaxis = list(title = "EIF")))

p5 <- plot_ly(x = 10:20, y = 10:35, z = as.matrix(EIF_results[[5]])) %>% 
  add_surface(contours = list(z = list(show=TRUE, usecolormap=TRUE, highlightcolor="#ff0000", project=list(z=TRUE)))) %>%
  layout(scene = list(camera=list(eye = list(x=2, y=1, z=0.3)),
                      xaxis = list(title = "log(MarketValue)"),
                      yaxis = list(title = "Age"),
                      zaxis = list(title = "EIF")))

p6 <- plot_ly(x = 10:20, y = 10:35, z = as.matrix(EIF_results[[6]])) %>% 
  add_surface(contours = list(z = list(show=TRUE, usecolormap=TRUE, highlightcolor="#ff0000", project=list(z=TRUE)))) %>%
  layout(scene = list(camera=list(eye = list(x=2, y=1, z=0.3)),
                      xaxis = list(title = "log(MarketValue)"),
                      yaxis = list(title = "Age"),
                      zaxis = list(title = "EIF")))















