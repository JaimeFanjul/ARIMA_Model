
## Jaime Fanjul García
# Time Series
# ARIMA Modeling

packages <- c("forecast","tseries","TTR")
new <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new)) install.packages(new)
a=lapply(packages, require, character.only=TRUE)

# clear workspace

rm(list=ls())

# Data weather from 1749-1983

datos = read.csv("C:/Users/YOUR PATH/data.csv")
datos

# Select the interest variable and transform it as time series
# indicate start date and frequency (in this case years).

datos <- ts(datos[, 2], start = c(1749,1), frequency = 12)

# clear values, as we will see to correct the variance to get homoscedasticity

datos <- replace(datos, datos == 0, 001)

# First of all graph the series

plot(datos, main = "Temperaturas min ( Grados Fahrenheit)",
     xlab = "Tiempo", ylab = "Temperaturas mín", col = "blue")

# it is possible to check that the time series is very closed to a stacionary series, 
# but the variance doesnt change constantly

# Try to find auto suggested differentiation 
ndif <- nsdiffs(datos)
ndif


# the variance is corrected taking log
logdatos=log10(datos)
plot(logdatos, main = "Temperaturas min ( Grados Fahrenheit)",
     xlab = "Tiempo", ylab = "Temperaturas min ( Grados Fahrenheit)", col = "blue")

ndif <- nsdiffs(logdatos)
ndif
# differentiation number 0, the time series looks like stacionary

# corroborate the deductions above with testing 
# ADF
adf.test(datos)
adf.test(logdatos)

# as a result of the test, the null hypothesis is rejected and the alternative hypothesis is not rejected.
# Supports and corroborates the graphs represented above.

# KPSS
.rs.restartR()
library(tseries)
kpss.test(datos, null = c("Level", "Trend"), lshort = TRUE)
kpss.test(logdatos, null = c("Level", "Trend"), lshort = TRUE)

# null hypothesis is rejected, both series are stationary
# by means of the two previous tests it is determined that the series are stationary,
# therefore these two tests would be enough.
# It would not be necessary to differentiate the series
# ERS

library(urca)
ers = ur.ers(datos, type = c("DF-GLS", "P-test"), model = c("constant", "trend"), lag.max = 4)
summary(ers)
ers = ur.ers(logdatos, type = c("DF-GLS", "P-test"), model = c("constant", "trend"), lag.max = 4)
summary(ers)

# PP
library(aTSA)
pp.test(datos, type = c("Z_rho", "Z_tau"), lag.short = TRUE, output = TRUE)
pp.test(logdatos, type = c("Z_rho", "Z_tau"), lag.short = TRUE, output = TRUE)


# Once we have reached this point it is time to pass the modeling
# The first step is to decide which transformation of the data we are going to model.
# We have seen that the best option is to take logarithms.

# To determine the "p" and "q" orders of the model:
# We can use the correlogram ("q" or MA part) and the partial correlogram (AR part)

par(mfrow = c(1,2))
Acf(datos, main = "Correlograma  (MA)")
Pacf(datos, main = "Correlograma Parcial (AR)" )

par(mfrow = c(1,2))
Acf(logdatos, main = "Correlograma  (MA)")
Pacf(logdatos, main = "Correlograma Parcial (AR)" )

### ARIMA Modeling

# The value where the partial correlogram is truncated shows the AR order.
# The value where the correlogram is truncated shows the MA order.

# In this way we can select the best model (or several and check which is the best)


train <- window(logdatos, start = c(1749,1), end = c(1980,12))
test <- window(logdatos, start = c(1981,1))
ARIMAfit <- auto.arima(train, approximation = FALSE, trace = FALSE)


# approximation = TRUE only for very long series or wide seasonal period
# trace = FALSE to not report all the considered models

summary(ARIMAfit)

# Recommend an ARIMA (2,0,0) x (0,0,3)
# There is no differentiation
# We make the prediction.

.rs.restartR()
library(forecast)
par(mfrow = c(1,1))
pred.log <- forecast(ARIMAfit, h = 12*3) # 3 years
pred.log

# We remove the logs from the prediction.

pred <- pred.log
for (columna in c("x","mean","lower","upper","fitted","residuals")) {
  if (is.ts(pred.log[[columna]]) | is.numeric(pred.log[[columna]]))
    pred[[columna]] = 10^(pred.log[[columna]])
}

pred[[4]]

# Graph

plot(pred[[4]], col = "red", lwd = 2, main = "Prediccion ARIMA")

lines(datos, col = "green", lwd = 2, lty = 2)


# Residuals analysing.

plot(ARIMAfit$residuals)

# We seek that the residuals follow the normal distribution, mean 0, homcedasticity and are independent

par(mfrow = c(1,2))
Acf(ts(ARIMAfit$residuals), main = "ACF Residual", col = "red")
Pacf(ts(ARIMAfit$residuals), main = "PACF Residual", col = "red")

# Most of the residues contained in the acceptance interval, close to 0

# We calculate the R ^ 2
R2=1-(sum((datos-pred[[4]])^2)/sum((datos-mean(datos))^2))
R2


# we obtain an r ^ 2 of 0.98, therefore our model is capable of
# explain 98% of the temperature variation


##################################################
# Model no log

train2 <- window(datos, start = c(1749,1), end = c(1980,12))
test2 <- window(datos, start = c(1981,1))
ARIMAfit2 <- auto.arima(train2, approximation = FALSE, trace = FALSE)

summary(ARIMAfit2)

# it is observed that the penalties increase
# a differentiation is made
# Recommend an ARIMA (2,1,0) x (0,1,2)
# We make the prediction.
.rs.restartR()
par(mfrow = c(1,1))
pred <- forecast(ARIMAfit2, h = 12*3) # 3 ainos
pred



plot(pred[[4]], col = "red", lwd = 2, main = "Prediccion ARIMA")

lines(datos, col = "green", lwd = 2, lty = 2)

# Lines does not work properly 


# Residuals

plot(ARIMAfit$residuals)

par(mfrow = c(1,2))
Acf(ts(ARIMAfit2$residuals), main = "ACF Residual", col = "red")
Pacf(ts(ARIMAfit2$residuals), main = "PACF Residual", col = "red")
# The residues get worse we have more significant residues not contained within the acceptance interval.


# R^2
R2=1-(sum((datos-pred[[4]])^2)/sum((datos-mean(datos))^2))
R2


# We get r ^ 2 of 0.97, it gets worse so that we can conclude:
# To improve the model it is decided to take
# log correcting the problem of non-homecedasticity detected in the representation of the series.
# It has been proven that the model taking log is more reliable, it is able to explain a greater variability of the
# the temperature over the years.
# On the other hand, it has been observed how the residuals improve for the model by taking log, being contained in
# in the assumption interval as 0.

# Regarding the decomposition of the time series in trend, seasonality and randomness, there is no
# It makes sense to do it in this exercise as it will not be used in modeling.
# We have obtained a more than acceptable model by taking log and correcting the variance that helps us to
# make predictions.