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


grepl_rows <- function(pattern, df) {
  cell_matches = apply(df, c(1,2), function(x) grepl(pattern, x))
  apply(cell_matches, 1, all)
}


flatten_matrix <- function(M, sep='.') {
  x <- c(M)
  names(x) <- c(outer(rownames(M), colnames(M),
                      function(...) paste(..., sep=sep)))
  x
}



get_position <- function(x)
  sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)

dr = function(selected) {
  positions = unique(get_position(selected)) # remove possible duplicates
  discoveries = length(positions)
  false_discoveries = length(setdiff(positions, tsm_df$Position))
  list(true_discoveries = discoveries - false_discoveries,
       false_discoveries = false_discoveries)
}



get_conventional_threshold <- function(S, q_loc_fdr){
  t = sort(abs(S))
  Ta = sapply(t,function(x){(1+length(S[S<(-x)]))/max(1,length(S[S>x]))})
  sh = min(t[which(Ta<=q_loc_fdr)]) 
  return(sh)
}


get_localfdr_threshold_half <- function(lfdr13, lfdr2, lfdr4,q_loc_fdr){
  t = sort(c(lfdr13,lfdr2,lfdr4))
  t[t==1]=0.9999
  Ta = sapply(t,function(x){(1+0.5*(length(lfdr4[lfdr4<=x]) + length(lfdr2[lfdr2<=x])))/max(1,length(lfdr13[lfdr13<=x]))})
  sh = max(t[which(Ta<=q_loc_fdr)]) 
  return(sh)
}



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

