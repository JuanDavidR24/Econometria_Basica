#### Mínimos Cuadrados Ordinarios.
#### Juan David Rincón, 2023

# 2. MODELO: CAPM AMAZON ----
# r_i = r_f + B(r_mkt - r_f)
# r_i - r_f = alpha + beta*(r_mkt - r_f) + e

## 2.1. Variables ---- 

#install.packages("readxl")
library(readxl)

# Importar datos: Datos CAPM Amazon.
archivo <- file.choose()
DATA <- as.data.frame(read_excel(archivo, sheet = 1))

# Datos (Primeros Valores).
head(DATA)

### 2.1.1. Representación Gráfica ----

#install.packages("ggplot2")
library(ggplot2)

# Gráfica Amazon Precios.
ggplot(data=DATA, aes(x = Fecha, y = AMZN)) + 
  geom_line(color = "blue") + 
  labs(x = "Año/Mes", y = "Precio $USD", 
       title = "Último precio mensual Amazon (AMZN)", 
       subtitle = "Enero 2015 - Diciembre 2022")

# Gráfica de Amazon Retornos.
ggplot(data=DATA, aes(x = Fecha, y = `AMZN Rtrn`)) + 
  geom_line(color = "blue") + 
  labs(x = "Año/Mes", y = "Retorno Mensual", 
       title = "Rentabilidad Mensual Amazon (AMZN)", 
       subtitle = "Enero 2015 - Diciembre 2022") +
  scale_y_continuous(labels = scales::percent)

# Gráfica NASDAQ Precios.
ggplot(data=DATA, aes(x = Fecha, y = NASDAQ)) + 
  geom_line(color = "red") + 
  labs(x = "Año/Mes", y = "Precio $USD", 
       title = "Último precio mensual Nasdaq Composite Index", 
       subtitle = "Enero 2015 - Diciembre 2022")

# Gráfica NASDAQ Retornos.
ggplot(data=DATA, aes(x = Fecha, y = `NASDAQ Rtrn`)) + 
  geom_line(color = "red") + 
  labs(x = "Año/Mes", y = "Retorno Mensual", 
       title = "Rentabilidad Mensual Nasdaq Composite Index", 
       subtitle = "Enero 2015 - Diciembre 2022") +
  scale_y_continuous(labels = scales::percent)

# Gráfica Bonos del Tesoro 10 años.
ggplot(data=DATA, aes(x = Fecha, y = `UST10 mv`)) + 
  geom_line(color = "black") + 
  labs(x = "Año/Mes", y = "Retorno Mensual", 
       title = "Tasa mensual Bonos del Tesoro de Estados Unidos a 10 años",
       subtitle = "Enero 2015 - Diciembre 2022") +
  scale_y_continuous(labels = scales::percent)

### 2.1.2. Construcción variables del modelo ----

# Amazon, Ri.
ri <- DATA$`AMZN Rtrn`

# Nasdaq, Rmkt.
rmkt <- DATA$`NASDAQ Rtrn`

# Bonos US 10y, Rf.
rf <- DATA$`UST10 mv`

# Prima del activo, Ri-Rf.
ri_rf <- ri - rf

# Prima del mercado, Rmkt-Rf.
rmkt_rf <- rmkt - rf

# 3. ESTIMACIÓN OLS, MCO ----
## 3.1. Matrices de Diseño ----

# Matriz Y (Variable dependiente).
Y <- as.matrix(ri_rf)

# Matriz X (Variables independientes).
X <- cbind(1, rmkt_rf)

# Número de datos y variables.
N <- dim(X)[1]
K <- dim(X)[2]

## 3.2. Mínimos Cuadrados Ordinarios ----

# Función OLS.
OLS.f <- function(X,Y){
  
  # Formula.
  B_gr <- solve(t(X)%*%X)%*%t(X)%*%Y
  
  # Resultado.
  return(B_gr)
}

# Betas Estimados.
B_gr <- OLS.f(X,Y)
B_gr

## 3.3. Estimación variable dependiente ----

# Variable dependiente estimada.
Y_gr <- X%*%B_gr
head(Y_gr)

## 3.4. Estimación de errores ----

# Errores estimados.
e_gr <- Y - Y_gr
head(e_gr)

## 3.5. Varianza estimada de los errores ----

# Función de varianza estimada.
sigma2_gr.f <- function(residuales, T_obs, K_var){
  
  # Formula.
  sigma2_gr <- as.numeric((t(residuales)%*%residuales)/(T_obs-K_var))
  
  # Resultado.
  return(sigma2_gr)
}

# Varianza estimada.
sigma2_gr <- sigma2_gr.f(e_gr, N, K)
sigma2_gr

# Desviación Estándar estimada.
sigma_gr <- sqrt(sigma2_gr)
sigma_gr

## 3.6. Matriz de Varianza-Covarianza de los betas ----

# Matriz de Varianza Covarianza.
varcov_beta_gr <- sigma2_gr*solve(t(X)%*%X)
varcov_beta_gr

# Desviaciones estándar de los betas.
sd_beta_gr <- sqrt(diag(varcov_beta_gr))
sd_beta_gr

## 3.7. Coeficiente de determinación ----

# Función Coeficiente de determinación R2.
R2.f <- function(Y, Y_gr){
  
  # Suma explicada de cuadrados.
  sec <- sum((Y_gr-mean(Y))^2)
  
  # Suma total de cuadrados.
  stc <- sum((Y-mean(Y))^2)
  
  # Formula.
  R2 <- as.numeric(sec/stc)
  
  # Resultado.
  return(R2)
}

# Coeficiente de determinación R2.
R2 <- R2.f(Y, Y_gr)
R2

### 3.7.1. Coeficiente de determinación ajustado ----

# Función R2 ajustado.
R2adj.f <- function(R2, T_obs, K_var){
  
  # Formula.
  R2adj <- 1-((1-R2)*((T_obs-1)/(T_obs-(K_var-1)-1)))
  
  # Resultado.
  return(R2adj)
}

# Coeficiente de determinación R2 ajustado.
R2adj <- R2adj.f(R2, N, K)
R2adj

# 4. USO DEL MODELO ----

# Data frame provisional.
df1 <- as.data.frame(cbind(ri_rf, rmkt_rf))

# Gráfica.
ggplot(df1, aes(x=rmkt_rf, y=ri_rf)) + 
  geom_point(color="blue", pch=18) + 
  geom_abline(color="red", intercept = B_gr[1], slope = B_gr[2]) +
  labs(x = "Prima de riesgo del mercado, (Rmkt - Rf)",
       y = "Prima de riesgo del activo, (Ri - Rf)",
       title = "CAPM Amazon",
       subtitle = paste("Y =", round(B_gr[1],6), "+", round(B_gr[2],6), "X")) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent)

## 4.1. Pronóstico ----

### 4.1.1. Pronóstico dentro de muestra ----

# Data frame provisional.
df2 <- data.frame(Fecha = DATA$Fecha, Y_gr = Y_gr, ri_rf = ri_rf)

# Gráfica pronóstico dentro de muestra.
ggplot(df2, aes(x=Fecha)) + 
  geom_line(aes(y=ri_rf, linetype="Efectiva", color="Efectiva")) +
  geom_line(aes(y=Y_gr, linetype="Estimación", color="Estimación")) + 
  labs(title = "Prima de riesgo Amazon: Efectiva vs Estimada",
       subtitle = "Enero 2015 - Diciembre 2022",
       x = "Año/Mes", y = "Prima de Riesgo del activo",
       linetype = "", color = "") +
  scale_linetype_manual(values = c("Efectiva"="solid", "Estimación"="dashed")) +
  scale_color_manual(values = c("Efectiva"="black", "Estimación"="red")) +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom")

# Promedio valor absoluto de los errores.
mean(abs(e_gr))

# Vector de coincidencias.
dCoinc <- (Y > 0) == (Y_gr > 0)

# Porcentaje de aciertos.
sum(dCoinc)/N

### 4.1.2. Pronóstico fuera de muestra ----

# Importar datos fuera de muestra.
DATA0 <- as.data.frame(read_excel(archivo, sheet = 2))

# Datos (Primeros Valores).
head(DATA0)

# Variables iniciales.
ri_0 <- DATA0$`AMZN Rtrn`
rmkt_0 <- DATA0$`NASDAQ Rtrn`
rf_0 <- DATA$`UST10 mv`

# Variables del modelo.
ri_rf_0 <- ri_0 - rf_0
rmkt_rf_0 <- rmkt_0 - rf_0

# Matrices de diseño.
Y0 <- as.matrix(ri_rf_0)
X0 <- cbind(1, rmkt_rf_0)

# Pronóstico.
Y0_gr <- X0%*%B_gr

# Data frame provisional.
df3 <- data.frame(Fecha = DATA0$Fecha, Y0_gr = Y0_gr, ri_rf_0 = ri_rf_0)

# Gráfica pronóstico fuera de muestra.
ggplot(df3, aes(x=Fecha)) + 
  geom_line(aes(y=ri_rf_0, linetype="Efectiva", color="Efectiva")) +
  geom_line(aes(y=Y0_gr, linetype="Estimación", color="Estimación")) + 
  labs(title = "Prima de riesgo Amazon: Efectiva vs Estimada",
       subtitle = "Enero 2013 - Diciembre 2014",
       x = "Año/Mes", y = "Prima de Riesgo del activo",
       linetype = "", color = "") +
  scale_linetype_manual(values = c("Efectiva"="solid", "Estimación"="dashed")) +
  scale_color_manual(values = c("Efectiva"="black", "Estimación"="red")) +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom")

# Promedio valor absoluto de los errores.
e0_gr <- Y0 - Y0_gr
mean(abs(e0_gr))

# Vector de coincidencias.
dCoinc0 <- (Y0 > 0) == (Y0_gr > 0)

# Porcentaje de aciertos.
sum(dCoinc0)/length(Y0)

# Vectores de información.
xene <- c(1, 0.1068 - 0.002894)
xfeb <- c(1, -0.0111 - 0.003206)
xmar <- c(1, 0.0669 - 0.002866)
xabr <- c(1, 0.0004 - 0.002832)

# Pronóstico.
pronM <- matrix(NA, nrow = 4, ncol = 2)
pronM[1,] <- c("Enero", round(xene%*%B_gr, 6))
pronM[2,] <- c("Febrero", round(xfeb%*%B_gr, 6))
pronM[3,] <- c("Marzo", round(xmar%*%B_gr, 6))
pronM[4,] <- c("Abril", round(xabr%*%B_gr, 6))
pronM
