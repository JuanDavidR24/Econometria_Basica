### OLS Function Pack 1.
## Juan David Rincón Mora, 2023.

# OLS Estimation Function ----
Beta_OLS.f <- function(X,Y){
  
  # Formula.
  beta_gr <- solve(t(X)%*%X)%*%t(X)%*%Y
  
  # Result.
  return(beta_gr)
}

# Estimated Variance Function ----
sigma2_gr.f <- function(errors, T_obs, K_var){
  
  # Formula.
  sigma2_gr <- as.numeric((t(errors)%*%errors)/(T_obs-K_var))
  
  # Result.
  return(sigma2_gr)
}

# R2 function ----
R2.f <- function(Y, Y_gr){
  
  # Sum of squares explained.
  ssr <- sum((Y_gr-mean(Y))^2)
  
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

