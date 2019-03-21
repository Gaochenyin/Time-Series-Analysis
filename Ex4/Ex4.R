setwd('C:\\Users\\Lenovo\\Desktop\\大四上\\时间序列\\4')
# Loading required packages
library(astsa)
library(fracdiff)
library(xtable)
library(latex2exp)
library(stats)
library(fGarch)
library(marima)
library(vars)
# Plot the time-line of arf
plot.ts(arf,
        col='steelblue')
par(mfrow=c(1,2))
# ACF and PACF analysis of arf
acf(arf,100,
    col='steelblue')
pacf(arf,100,
    col='steelblue')
# Begian ARFIMA model fit 
arf.mean <- arf-mean(arf)
arf.fd <- fracdiff(arf.mean,nar=1,
                   ar=1,h=0,M=500)
arf.fd.sum <- summary(arf.fd)
# Output LaTex code
xtable(arf.fd.sum$coef,
       digits = 4)
# Display model coefficients and Plot it
p <- rep(1,50)
for (k in 1:50){ p[k+1]=(k-arf.fd$d)*p[k]/(k+1) }
plot(1:50, p[-1],
     ylab=expression(pi(d)),
     xlab="Index", 
     type="h")
# Plot the first difference of arf
plot.ts(diff(arf),
        col='steelblue',
        ylab=TeX('First Difference of x_t'))
par(mfrow=c(1,2))
# ACF and PACF of first difference of arf
acf(diff(arf),100,
    col='steelblue')
pacf(diff(arf),100,
    col='steelblue')
# ARMA model fit to arf
arf.arma.110 <- arima(arf,order=c(1,1,1))
arf.arma.predicted <- arf.arma.110$residuals + arf
# Compared the true and the prediction of time series
plot(arf,col='steelblue',lwd=1)
lines(arf.arma.predicted,col='orange',lwd=1)
legend('topleft',
       legend = c('Real',
                  'Predicted'),
       col=c('steelblue',
             'orange'),
       lwd=c(2,2),
       lty=c(1,1))
# Output LaTex Code
arf.arma.ps <- rbind(arf.arma.110$coef,
      arf.arma.110$sigma2)
xtable(arf.arma.110$var.coef,digits = 4)

###################
# Plot of the European Stock Markets
plot.ts(EuStockMarkets,
        col='steelblue')
# For example
DAX <- EuStockMarkets[,1]
# Compute the return of the time series
DAX.r <- diff(log(DAX))
plot(DAX.r,type='l')
# fit GARCH(1,1) model for the return
DAX.garch <- garchFit(~garch(1,1),
         data=DAX.r,
         cond.dist = 'norm')
# Output the LaTex Code
xtable(DAX.garch@fit$matcoef,
       digits=4)
# Plot the estimated results
plot.ts(DAX.r,
     type='l',
     ylab='',
     col='darkgrey',
     col.axis='darkgrey')
par(new=T)
DAX.h <- ts(DAX.garch@fit$series$h,
   start = start(DAX),
   end = end(DAX),
   frequency = frequency(DAX))
plot.ts(DAX.h,
     type='l',
     col='steelblue',
     lwd=1,axes=F,
     ylab='')
axis(4,col="steelblue",
     col.ticks="steelblue",
     col.axis="steelblue")

########
# Plot all the climate time series
plot.ts(climhyd,
        col='steelblue')
# Plot the Precipitation and the prediction by SARIMA model
Precip.season <- sarima(sqrt(climhyd$Precip),
       p=0,d=0,q=0,
       P=0,D=1,Q=1,S=12)
Precip.season.predict <- sqrt(climhyd$Precip)+Precip.season$fit$residuals
plot.ts(sqrt(climhyd$Precip),
        col='steelblue',
        ylab='Precip',ylim=c(-10,50))
lines(Precip.season.predict,
      col='orange')
# Plot the inflow series and prediction by SARIMA model
inflow.season <- sarima(log(climhyd$Inflow),
                        p=0,d=0,q=0,
                        P=0,D=1,Q=1,S=12)
inflow.season.predict <- log(climhyd$Inflow)+inflow.season$fit$residuals
plot.ts(log(climhyd$Inflow),
        col='steelblue',
        ylab='Inflow',
        ylim=c(3,9))
lines(Precip.season.inflow.predict,
      col='orange')

# Prewhittened version of inflow and Plot
inflow.pw <- inflow.season$fit$residuals
plot(inflow.pw ,
     col='steelblue',
     ylab='Prewhitened flow residuals')
sma1 <- Precip.season.inflow$ttable[1]
# Construct the filter parameters
season_AR <- sapply(0:10,function(x)0.9815^x)*-0.0185
int <- rep(0,11)
AR_par <- c(1,int)
for(i in season_AR)
{
  AR_par <- c(AR_par,i,int)
}
pre.fil <- filter(sqrt(climhyd$Precip),
                  filter = AR_par,
                  sides=1)
# Plot the filtered transformed precipitation series
plot(pre.fil,
     col='steelblue',
     ylab='Filtered Precipitation series')
par(fig = c(0,0.35,.5,1),new = TRUE)
plot(AR_par,
     col='orange',
     type='l',xaxt='n',
     yaxt='n',
     ylab='',
     xlab='')
# Compute the CCF of the prewhitten version of inflow and filted transformed precipitation
ccf(inflow.pw,pre.fil,
    na.action=na.omit,
    panel.first=grid())

# 5.12

# Plot the original time series
matplot(econ5[,1:3],
     type='l',
     col= c('steelblue','orange','green'),
     ylab='Unemploy&GNP&Consum',
     lwd=2)
# Log transformation of econ5
econ5.log <- log(econ5)
m <- dim(econ5.log)[1]
time <- 1:m
# Least-square detrended time series
econ5.log.detrend <- apply(econ5.log,2,function(x)
  {
  fit <- lm(x~time)
  return(fit$residuals)
})
# Plot the detrended time series
matplot(econ5.log.detrend[,1:3],
     type='l',
     col= c('steelblue','orange','green'),
     ylab='Detrended Unemploy&GNP&Consum',
     lwd=2)
econ5.log.detrend.3 <- econ5.log.detrend[,1:3]

# VAR(econ5.log.detrend.3,lag.max = 100,ic='AIC')
# Define VARMA(1,1) model
model <- define.model(kvar = 3,ar=c(1),ma=c(1))
arp <- model$ar.pattern
map <- model$ma.pattern
model.fit <- marima(econ5.log.detrend.3,
       ar.pattern = arp,
       ma.pattern = map)
# Obtain the fitted residuals
model.residual <- t(resid(model.fit))
colnames(model.residual) <- c('Unemploy',
                              'GNP',
                              'Consum')
# Plot the model Residuals
plot.ts(model.residual,
        col='steelblue',
        main='Model Residual')
# Output the LaTex Code
xtable(model.fit$ar.estimates[,,2],digits = 4)
xtable(model.fit$ma.estimates[,,2],digits=4)
model.fit$Constant
# Run ACF and PACF test of model residuals
par(mfrow=c(3,2))
acf(model.residual[,1][-1],
    main='Unemploy')
pacf(model.residual[,1][-1],
    main='Unemploy')
acf(model.residual[,2][-1],
    main='GNP')
pacf(model.residual[,2][-1],
     main='GNP')
acf(model.residual[,3][-1],
    main='Consum')
pacf(model.residual[,3][-1],
     main='Consum')
