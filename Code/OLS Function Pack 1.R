### OLS Function Pack 1.
## Juan David Rincón Mora, 2023.

# OLS Estimation Function ----
Beta_OLS.f <- function(X,Y){
  
  # Formula.
  beta_ht <- solve(t(X)%*%X)%*%t(X)%*%Y
  
  # Result.
  return(beta_ht)
}

# Estimated Variance Function ----
sigma2_ht.f <- function(errors, T_obs, K_var){
  
  # Formula.
  sigma2_ht <- as.numeric((t(errors)%*%errors)/(T_obs-K_var))
  
  # Result.
  return(sigma2_ht)
}

# R2 function ----
R2.f <- function(Y, Y_ht){
  
  # Sum of squares explained.
  ssr <- sum((Y_ht-mean(Y))^2)
  
  # Total sum of squares.
  sst <- sum((Y-mean(Y))^2)
  
  # Formula.
  R2 <- as.numeric(ssr/sst)
  
  # Result.
  return(R2)
}

# Adjusted R2 function ----
R2adj.f <- function(R2, T_obs, K_var){
  
  # Formula.
  R2adj <- 1-((1-R2)*((T_obs-1)/(T_obs-(K_var-1)-1)))
  
  # Result.
  return(R2adj)
}

