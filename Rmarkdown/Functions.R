#-------- Returns rows for which every column matches the given regular expression --------#
# Reference: https://web.stanford.edu/group/candes/knockoffs/software/knockoffs/tutorial-4-r.html
grepl_rows <- function(pattern, df) {
  cell_matches = apply(df, c(1,2), function(x) grepl(pattern, x))
  apply(cell_matches, 1, all)
}


#-------- Flatten a matrix to a vector with names from concatenating row/column names --------#
# Reference: https://web.stanford.edu/group/candes/knockoffs/software/knockoffs/tutorial-4-r.html
flatten_matrix <- function(M, sep='.') {
  x <- c(M)
  names(x) <- c(outer(rownames(M), colnames(M),
                      function(...) paste(..., sep=sep)))
  x
}


#-------- Extract numeric values and convert them to numeric type --------#
# Reference: https://web.stanford.edu/group/candes/knockoffs/software/knockoffs/tutorial-4-r.html
get_position <- function(x){
  sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)
}


#-------- Computes the number of true and false discoveries --------#
# Reference: https://web.stanford.edu/group/candes/knockoffs/software/knockoffs/tutorial-4-r.html
dr = function(selected) {
  positions = unique(get_position(selected)) # remove possible duplicates
  discoveries = length(positions)
  false_discoveries = length(setdiff(positions, tsm_df$Position))
  list(true_discoveries = discoveries - false_discoveries,
       false_discoveries = false_discoveries)
}


#-------- Computes the precision matrix --------#
#
# Parameters:
#   data: A design matrix where rows represent observations and columns represent variables.
#
# Returns: 
#   The precision matrix.
#
# Reference: https://github.com/Jeremy690/-A-Scale-free-Approach-for-False-Discovery-Rate-Control-in-Generalized-Linear-Models/blob/main/code/simulation/figure5/DS_MDS.R
get_M = function(data){
  ## data: design matrix
  n = dim(data)[1]
  p = dim(data)[2]
  C = matrix(0, nrow = p, ncol = p)
  tau_square = numeric(p)
  ## Nodewise Lasso regression
  for(j in 1:p){
    y = data[,j]
    X <- data[, -j]
    ## To save computation time, we only do the cv once
    if(j == 1){
      cvfit <- cv.glmnet(X, y, nfolds = 10, nlambda = 200, intercept = F, standardize = F)
      ## We will use this same lambda for the following Lasso regression
      lambda <- cvfit$lambda.min
    }
    beta1 <- as.vector(glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = F, standardize = F)$beta)
    C[j, -j] = -beta1
    tau_square[j] = mean((y-X%*%beta1)*y)
  }
  diag(C) = 1
  T_mat = diag(tau_square)
  M = solve(T_mat)%*%C
  M
}


#-------- Generates multivariate normal random samples --------#
#
# Parameters:
#   n: Number of samples to generate (default is 1).
#   mu: Mean vector.
#   Sigma: A positive-definite symmetric matrix specifying the covariance matrix of the variables.
#   tol: Tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma.
#
# Returns: 
#   A matrix where each row is a sample from the multivariate normal distribution.
mvrn = function (n = 1, mu, Sigma, tol = 1e-06) {
  p <- length(mu)
  eS1 = eS =  eigen(Sigma, symmetric = TRUE)
  for (i in 1:ncol(eS$vectors)) {
    if (sign(eS1$vectors[1,i]) != 0){
      eS$vectors[,i] =  eS1$vectors[,i] * sign(eS1$vectors[1,i])  
    }
  }
  ev <- eS$values
  X <- matrix(rnorm(p * n), n)
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  return(t(X))
}


#-------- Computes the conventional threshold for mirror statistics --------#
#
# Parameters:
#   S: A numeric vector of mirror statistics.
#   q_fdr: Target FDR level.
#   constant: A constant for FDP estimate.
#
# Returns: 
#   The rejection threshold for mirror statistics.
get_conventional_threshold <- function(S, q_fdr, constant){
  t = sort(abs(S))
  Ta = sapply(t,function(x){(constant+length(S[S<(-x)]))/max(1,length(S[S>x]))})
  sh = min(t[which(Ta<=q_fdr)])
  return(sh)
}


#-------- Computes the local FDR threshold for grouped mirror statistics --------#
#
# Parameters:
#   lfdr13: A numeric vector of local FDR values for points in 1st and 3rd quadrants.
#   lfdr2: A numeric vector of local FDR values for points in 2nd quadrant.
#   lfdr4: A numeric vector of local FDR values for points in 4th quadrant.
#   q_fdr: Target FDR level.
#   constant: A constant for FDP estimation.
#
# Returns: 
#   The rejection threshold for the symmetric statistics.
get_localfdr_threshold_half <- function(lfdr13, lfdr2, lfdr4,q_fdr,constant){
  t = sort(c(lfdr13,lfdr2,lfdr4))
  t[t==1]=0.9999
  Ta = sapply(t,function(x){(constant+0.5*(length(lfdr4[lfdr4<=x]) + length(lfdr2[lfdr2<=x])))/max(1,length(lfdr13[lfdr13<=x]))})
  sh = max(t[which(Ta<=q_fdr)]) 
  return(sh)
}


#-------- Selection based on the Symmetric Adaptive Selection procedure --------#
#
# Parameters:
#   X1: A numeric vector of the 1st dimension of symmetric statistics.
#   X2: A numeric vector of the 2nd dimension of symmetric statistics.
#   label_X1X2: Labels for null and non-null features.
#   c: A constant for FDP estimation.
#   c0: A scaling factor for the bandwidth of \hat_{f}_{0}(t).
#   c1: A scaling factor for the non-null bandwidth \hat_{f}(t).
#   q_fdr: Target FDR level.
#
# Returns: 
#   A vector of indices that pass the local FDR threshold based on the SAS procedure.
SAS = function(X1,X2,label_X1X2,c,c0=1,c1=1.4,q_fdr = 0.1){
  index1 = intersect(which(X1>0),which(X2>0))
  index2 = intersect(which(X1<0),which(X2>0))
  index3 = intersect(which(X1<0),which(X2<0))
  index4 = intersect(which(X1>0),which(X2<0))
  
  X1_4444 = c(X1[c(index2,index4)])
  X2_4444 = c(X2[c(index2,index4)])
  
  h0 = c0*c(bandwidth.nrd(X1_4444), bandwidth.nrd(X2_4444))
  h1 = c1*c(bandwidth.nrd(X1), bandwidth.nrd(X2))
  
  Kernel_h0 = function(x){0.3*sum(dnorm((x[1] - X1_4444)/h0[1])*dnorm((x[2] - X2_4444)/h0[2]))/(h0[1]*h0[2]*length(X1_4444))}
  Kernel_h1 = function(x){sum(dnorm((x[1] - X1)/h1[1])*dnorm((x[2] - X2)/h1[2]))/(h1[1]*h1[2]*length(X1))}
  f0_hat0 = function(x){
    return(0.25*(Kernel_h0(x) + Kernel_h0(c(-x[1],x[2])) + Kernel_h0(c(x[1],-x[2])) + Kernel_h0(-x)))
  }
  f0_hat = function(x){
    return(0.5*(f0_hat0(x) + f0_hat0(c(x[2],x[1]))))
  }
  f_hat = function(x){
    return(0.5*(Kernel_h1(x) + Kernel_h1(c(x[2],x[1]))))
  }
  Lfdr_hat13 = function(x){min(1,(f0_hat(x)/f_hat(x)))}
  
  lfdr_points_2 = rep(1,2*length(index2))
  lfdr_points_4 = rep(1,2*length(index4))
  lfdr_points13 = lfdr_points13 = rep(1,length(c(index1,index3)))
  lfdr_points13 = apply(cbind(X1[c(index1,index3)],X2[c(index1,index3)]), 1, Lfdr_hat13)
  lfdr_points_2= apply(cbind(c(X1[index2],-X1[index2]),c(-X2[index2], X2[index2])), 1, function(x){min(1,f0_hat(x)/f_hat(x))})
  lfdr_points_4= apply(cbind(c(-X1[index4],X1[index4]),c(X2[index4], -X2[index4])), 1, function(x){min(1,f0_hat(x)/f_hat(x))})
  
  
  fdr_q = get_localfdr_threshold_half(lfdr_points13, lfdr_points_2, lfdr_points_4, q_fdr, c)
  index13 = c(index1,index3)
  
  return(index13[which(lfdr_points13<=fdr_q)])
}


#-------- Selection based on the Derandomized SAS procedure --------#
#
# Parameters:
#   I_vec: A numeric vector of Derandomized statistics.
#   label_X1X2: Labels for null and non-null features.
#   q_fdr: Target FDR level.
#
# Returns: 
#   A vector of selected indices.
MSAS = function(I_vec,label_X1X2,q_fdr){
  I_sort = I_vec[order(I_vec)]
  label_sort = label_X1X2[order(I_vec)]
  I_sum = cumsum(I_sort)
  number_d = length(which(I_sum>q_fdr))
  number_td = length(intersect(which(I_sum>q_fdr), which(label_sort != 'Null')))
  fdr = (number_d - number_td)/ max(1,number_d)
  power = (number_td)/(length(which(label_X1X2 != 'Null')))
  return(c(fdr,power))
}


#-------- True discoveries and false discoveries based on the Derandomized SAS procedure for HIV data --------#
#
# Parameters:
#   I_vec: A numeric vector of Derandomized statistics.
#   label_X1X2: Labels for null and non-null features.
#   q_fdr: Target FDR level.
#
# Returns: 
#   The proxy of true discoveries and false discoveries.
MSAS_HIV = function(I_vec,label_X1X2,q_fdr){
  I_sort = I_vec[order(I_vec)]
  label_sort = label_X1X2[order(I_vec)]
  I_sum = cumsum(I_sort)
  discoveries = label_sort[which(I_sum>q_fdr)]
  result = dr(discoveries)
  return(result)
}


#-------- Selection based on the Derandomized SAS procedure for HIV data --------#
#
# Parameters:
#   I_vec: A numeric vector of Derandomized statistics.
#   label_X1X2: Labels for null and non-null features.
#   q_fdr: Target FDR level.
#
# Returns: 
#   Discovered features
MSAS_HIV_stability = function(I_vec,label_X1X2,q_fdr){
  I_sort = I_vec[order(I_vec)]
  label_sort = label_X1X2[order(I_vec)]
  I_sum = cumsum(I_sort)
  discoveries = label_sort[which(I_sum>q_fdr)]
  return(discoveries)
}


#-------- Selection based on the Derandomized SAS procedure for Supermarket data --------#
#
# Parameters:
#   I_vec: A numeric vector of Derandomized statistics.
#   q_fdr: Target FDR level.
#
# Returns: 
#   The number of discoveries.
MSAS_Sup = function(I_vec,q_fdr){
  label_X1X2 = c(1:length(I_vec))
  I_sort = I_vec[order(I_vec)]
  label_sort = label_X1X2[order(I_vec)]
  I_sum = cumsum(I_sort)
  discoveries = label_sort[which(I_sum>q_fdr)]
  return(length(discoveries))
}


