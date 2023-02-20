##### EJEMPLO AUTOCORRELACIÓN (JUDGE ET AL, 1982).

### 1. IMPORTAR DATOS. ----
# install.packages("readxl")
library(readxl)

DATA <- as.data.frame(read_excel(file.choose()))
head(DATA)

### 2. MATRICES DE DISEÑO. ----
Y <- as.matrix(DATA$y)
head(Y)

X <- cbind(1, DATA$x2, DATA$x3)
head(X)

N <- length(Y)
K <- dim(X)[2]

### 3. ESTIMACIÓN INICIAL. ----

## Estimación Betas.
OLS.f <- function(X,Y){
  B_gr <- solve(t(X)%*%X)%*%t(X)%*%Y
  return(B_gr)
}

B_gr <- OLS.f(X,Y)
B_gr

## Estimación de la Variable Dependiente.
Y_gr <- X%*%B_gr
head(Y_gr)

## Estimación de los errores.
e_gr <- Y - Y_gr
head(e_gr)

## Estimación de la varianza.
sigma2_gr.f <- function(residuales, T_obs, K_var){
  sigma2_gr <- as.numeric((t(residuales)%*%residuales)/(T_obs-K_var))
  return(sigma2_gr)
}

sigma2_gr <- sigma2_gr.f(e_gr, N, K)
sigma2_gr

sd_gr <- sqrt(sigma2_gr)
sd_gr

## Estimación Matriz Var-Cov de los Betas.
varcov_betas_gr <- sigma2_gr*solve(t(X)%*%X)
varcov_betas_gr

sd_betas_gr <- sqrt(diag(varcov_betas_gr))
sd_betas_gr

## Coeficiente de Determinación.
R2.f <- function(Y, Y_gorro){
  sec <- sum((Y_gorro-mean(Y))^2)
  stc <- sum((Y-mean(Y))^2)
  
  R2 <- as.numeric(sec/stc)
  return(R2)
}

R2 <- R2.f(Y, Y_gr)
R2

### 4. IDENTIFICACIÓN AUTOCORRELACIÓN. ----

## Test de Durbin Watson (Método Largo).
dwL.f <- function(residuales){
  N <- length(residuales)
  
  dw_num <- c()
  for(i in 2:N){
    dw_num[i-1] <- (residuales[i]-residuales[i-1])^2
  }
  
  dwL <- sum(dw_num)/sum(residuales^2)
  
  return(dwL)
}

DWL <- dwL.f(e_gr)
DWL

# Los valores críticos con 5% de significancia son: 
# 0 (AC Pos) 1.100 (Inc) 1.537 (No AC) 2.463 (Inc) 2.900 (AC Neg) 4

# Al caer el estadístico en la zona de Autocorrelación Positiva se rechaza 
# la H0, por lo tanto, hay autocorrelación.

## Estimación de Rho (Método Corto).
rho_gr <- 1 - (0.5*DWL)
rho_gr

### 5. CORREGIR AUTOCORRELACIÓN. ----

## Estimación Matriz de Choleski para Autocorrelación.
choleskyAC.f <- function(rho, N){
  cholesky <- matrix(data=0, nrow=N, ncol=N)
  
  diag(cholesky) <- 1
  cholesky[1,1] <- sqrt(1-(rho^2))
  
  for(i in 2:N){
    cholesky[i,(i-1)] <- -rho
  }
  
  return(cholesky)
}

P_gr <- choleskyAC.f(rho_gr, N)
View(P_gr)

## Estimación Matriz Psi Inversa para Autocorrelación.
psiInvAC.f <- function(rho, N){
  psiM <- matrix(data=0, nrow=N, ncol=N)
  
  diag(psiM) <- 1+(rho^2)
  psiM[1,1] <- 1
  psiM[N,N] <- 1
  
  for(i in 2:N){
    psiM[i,(i-1)] <- -rho
    psiM[(i-1),i] <- -rho
  }
  
  return(psiM)
}

psiInv_gr <- psiInvAC.f(rho_gr, N)
View(psiInv_gr)

## Variables Modificadas.
Y_star <- as.matrix(P_gr%*%Y)
head(Y_star)

X_star <- as.matrix(P_gr%*%X)
head(X_star)

## Estimación Beta Star.
B_star <- OLS.f(X_star, Y_star)
B_star

## Estimación FGSL.
FGLS.f <- function(X, Y, psiInv){
  B_fgls <- solve(t(X)%*%psiInv%*%X)%*%t(X)%*%psiInv%*%Y
  return(B_fgls)
}

B_fgls <- FGLS.f(X, Y, psiInv_gr)
B_fgls

## Variables Modificadas (Cuchilla de Cochrane-Orcutt).
Y_star0 <- as.matrix(Y_star[-1,])
head(Y_star0)

X_star0 <- X_star[-1,]
head(X_star0)

## Estimación Beta Star 0.
B_star0 <- OLS.f(X_star0, Y_star0)
B_star0

## Estimación FGSL 0.
B_fgls0 <- FGLS.f(X, Y, t(P_gr[-1,])%*%P_gr[-1,])
B_fgls0

### 6. NUEVAS ESTIMACIONES. ----

## Estimaciones con info completa.
Y_star_gr <- X_star%*%B_star
head(Y_star_gr)

e_star_gr <- Y_star - Y_star_gr
head(e_star_gr)

sigma2_2gr <- sigma2_gr.f(e_star_gr, N, K)
sigma2_2gr

sd_2gr <- sqrt(sigma2_2gr)
sd_2gr

varcov_betas_2gr <- sigma2_2gr*solve(t(X_star)%*%X_star)
varcov_betas_2gr

sd_betas_2gr <- sqrt(diag(varcov_betas_2gr))
sd_betas_2gr

## Estimaciones con Cuchilla de Cochrane-Orcutt.
Y_star_0gr <- X_star0%*%B_star0
head(Y_star_0gr)

e_star_0gr <- Y_star0 - Y_star_0gr
head(e_star_0gr)

sigma2_20gr <- sigma2_gr.f(e_star_0gr, N-1, K)
sigma2_20gr

sd_20gr <- sqrt(sigma2_20gr)
sd_20gr

varcov_betas_20gr <- sigma2_20gr*solve(t(X_star0)%*%X_star0)
varcov_betas_20gr

sd_betas_20gr <- sqrt(diag(varcov_betas_20gr))
sd_betas_20gr

### 7. PRONÓSTICO. ----

## Dentro de muestra (Para obs. 17).
Y[17,]
X[17,]%*%B_star + rho_gr*(Y[16,]-X[16,]%*%B_star)
X[17,]%*%B_star0 + rho_gr*(Y[16,]-X[16,]%*%B_star0)

## Fuera de muestra (Para obs. 21).
# Suponiendo X[21] = (1, 20, 20)
c(1, 20, 20)%*%B_star + rho_gr*(Y[20,]-X[20,]%*%B_star)
c(1, 20, 20)%*%B_star0 + rho_gr*(Y[20,]-X[20,]%*%B_star0)
