library(astsa)
library(tseries)
library(Rmisc)
library(ggplot2)
library(reshape2)
library(xtable)
my_theme <- theme_bw()+
theme(panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	panel.border = element_blank())
par(mai=c(0.5,0.5,0.5,1.5))
xy <- par('usr')
# Loading Original Data
sunspot <- read.table('sunspot.dat')
sun_ts <- ts(sunspot,
             start = 1700,end=1984,
             frequency =1)
len <- length(sun_ts)
plot.ts(sun_ts,
     col='steelblue',
     ylab='Sunspot',xlab='Date',
     lwd=2)
# compute ACF and PACF
acf(sun_ts,20,main='',lwd=2)
pacf(sun_ts,20,main='',lwd=2)
adf.test(sun_ts) # ADF test of original series

# Select p and q for seasonal ARMA(p,q)*(1,0) model
ps <- seq(4)
qs <- seq(4)
BIC <- matrix(0,nrow=4,ncol=4)
AIC <- matrix(0,nrow=4,ncol=4)
ESS <- matrix(0,nrow=4,ncol=4)
R2 <- matrix(0,nrow=4,ncol=4)
SS <- mean(var(sun_ts))
for(p in ps)
  for(q in qs)
  {
    mod <- sarima(sun_ts,p,0,q,1,0,0,11,details=FALSE)
    error <- sum(mod$fit$residuals^2,na.rm=TRUE)/(len-cad-1)
    BIC[p,q] <- mod$BIC
    AIC[p,q] <- mod$AIC
    ESS[p,q] <- error
    R2[p,q] <- 1-ESS[p,q]/SS
  }

# Define the plot_matrix function for selection
plot_matrix <- function(mat,title)
{
  frame <- melt(mat)
  p <- ggplot(frame,aes(x=Var1,y=Var2))+
    geom_tile(aes(fill=value),colour='white')+
    scale_fill_gradient(low = "white",high = "#1a476f")+
    ylab('p')+xlab('q')+ggtitle(title)+my_theme
  return(p)
}
p1 <- plot_matrix(BIC,'BIC(p,q)')
p2 <- plot_matrix(AIC,'AIC(p,q)')
p3 <- plot_matrix(ESS,'ESS(p,q)')
p4 <- plot_matrix(R2,'R2(p,q)')
multiplot(p1,p2,p3,p4,cols=2)
# Select the optimal lag p, q and Output LaTex code 
p_best <- 3
q_best <- 4
mod_best_sun <- sarima(sun_ts,p_best,0,q_best,1,0,0,11)
xtable(mod_best_sun$ttable,digits = 4)
# Obtain the Predict results
pre_ts_all <- sun_ts+mod_best$fit$residuals
# Compare the Prediction with the Real series
plot.ts(sun_ts,
        col='steelblue',
        ylab='Sunspot',
        xlab='Date',
        ylim=c(-20,250),
        lwd=2)
lines(pre_ts_all,
      col='orange')
legend(x=xy[2],
	y=xy[4]-yinch(0.8),
	legend=c('True',
        'Predicted'),
	xpd=TRUE,
	bty='n',
	col=c('steelblue','orange'),
	lty=1,
	lwd=2,
	title = '')
# Forecast the sunspots out for 4 years from 1985 to 1988
Sunspot <- sun_ts
pre_ts <- sarima.for(Sunspot ,4,p_best,0,q_best,1,0,0,11,
                     plot.all = FALSE)
# Loading real series from 1985 to 1988
sun_ts_all <- read.table('sunspot2.dat')
sun_ts_all <- ts(sun_ts_all,frequency=1,
                 start = 1700,end=1988)
real_four <- sun_ts_all[286:289]

# Compare the Forecast with Real series in plot with 95% CI
predict_four <- pre_ts$pred
predict_four_up <- pre_ts$pred+1.96*pre_ts$se
predict_four_low <- pre_ts$pred-1.96*pre_ts$se
real_predict_four <- cbind(real_four,predict_four)
plot.ts(real_predict_four[,1],
        type='l',
        col='steelblue',
        ylab='Sunspot',
        xlab='Date',
        ylim=c(-20,200),
        lwd=2,xaxt='n')
axis(side=1,at=1985:1988,labels=1985:1988)
lines(real_predict_four[,2],
      col='orange',
      lwd=2)
lines(predict_four_up,lty=2,
      col='grey',lwd=2)
lines(predict_four_low,lty=2,
      col='grey',lwd=2)
legend(x=xy[2],
	y=xy[4]-yinch(0.8),
	legend=c('True',
		'Predicted',
		'CI'),
	xpd=TRUE,
	bty='n',
    col=c('steelblue','orange','grey'),
    lty=c(1,1,3),
    lwd=c(1,1,2),
    title = '')

# Combine the forecast Details and Output LaTex code
pre_info_sun <- data.frame(cbind(real_four,
                             pre_ts$pred,
                             pre_ts$se,
                             predict_four_up,
                             predict_four_low))
colnames(pre_info_sun) <- c('True','Predicted',
                        'Standard Deviation','Up Bounds',
                        'Low Bounds')
xtable(pre_info_sun,digits = 3)

######3.35
lag_ADF = (length(sales)-1)^(1/3)
adf.test(sales) # ADF test of Original series
# Overall Plot 
plot.ts(sales,
        col='steelblue',
        ylab='Sales',xlab='Timeline',
        lwd=2)

# Fit linear regression model
fit <- lm(sales~time(sales),na.action = NULL)
# Plot residuals and its ACF,PACF at preset layout
plot_fit <- function(fit,ylab)
{
  layout(matrix(c(1,1,1,1,1,1,1,1,2,2,3,3,2,2,3,3),ncol=4,byrow=TRUE))
  plot.ts(fit$residuals,
          col='steelblue',
          ylab=ylab,xlab='Timeline',
          lwd=2)
  acf(fit$residuals,150,main=ylab)
  pacf(fit$residuals,150,main=ylab)
  
}
# ADF of model residuals
adf.test(fit$residuals)
# Plot differenced series with its ACF and PACF
plot.ts(diff(sales),
        col='steelblue',
        ylab='First Difference',xlab='Timeline',
        lwd=2)
acf(diff(sales),150,main='First Difference Time Series')
pacf(diff(sales),150,main='First Difference Time Series')
adf.test(diff(sales)) # ADF of the differenced series
sales_diff <- diff(sales)
# Define model selection function and ggplot 
statistic_test <- function(dat)
{
  ps <- seq(4)
  qs <- seq(4)
  BIC <- matrix(0,nrow=4,ncol=4)
  AIC <- matrix(0,nrow=4,ncol=4)
  ESS <- matrix(0,nrow=4,ncol=4)
  R2 <- matrix(0,nrow=4,ncol=4)
  SS <- mean(var(sales_diff))
  for(p in ps)
    for(q in qs)
    {
      mod <- sarima(dat,p,0,q,details=FALSE)
      error <- sum(mod$fit$residuals^2,na.rm=TRUE)/(len-cad-1)
      BIC[p,q] <- mod$BIC
      AIC[p,q] <- mod$AIC
      ESS[p,q] <- error
      R2[p,q] <- 1-ESS[p,q]/SS
    }
  p1 <- plot_matrix(BIC,'BIC(p,q)')
  p2 <- plot_matrix(AIC,'AIC(p,q)')
  p3 <- plot_matrix(ESS,'ESS(p,q)')
  p4 <- plot_matrix(R2,'R2(p,q)')
  multiplot(p1,p2,p3,p4,cols=2) 
}

# choose the best model and Output the LaTex code
mod_best <- sarima(sales_diff,4,0,3,
                   no.constant = TRUE)
xtable(mod_best$ttable,digits=4)
lead_diff <- diff(lead)
# CCF of differenced sales and lead
ccf_value <- ccf(sales_diff,lead_diff,main='Sales v.s. Lead',
    ylab='CCF')
# the maximum value of CCF
ccf_value$lag[which.max(ccf_value$acf)]
# lag plot of differenced lead and sales
lag2.plot(lead_diff,sales_diff,8)
# Fit differenced sales of differenced lead at lag 3
lag_reg <- ts.intersect(y=sales_diff,x=lag(lead_diff,3),dframe = T)
fit_lag <- lm(y~x,lag_reg)
# Residual Analysis of fitted model
plot_fit(fit_lag,'Residuals Analysis')
statistic_test(fit_lag$residuals)
# Select best lag of the residuals and Output LaTex code
mod_best <- sarima(fit_lag$residuals,
                   3,0,3)
mod_all_best <- sarima(lag_reg$y,3,0,3,xreg = lag_reg$x)
xtable(mod_all_best$ttable,digits = 3)
mod_pre <- mod_all_best$fit$residuals+lag_reg$y
mod_pre_limit <- mod_pre+1.96*mod_all_best$fit$sigma2
mod_pre_low <- mod_pre-1.96*mod_all_best$fit$sigma2
# plot predict with real data and 95% CI
plot.ts(lag_reg$y,
        col='steelblue',
        lwd=2,
        xlab='Time',ylab='Values',
        ylim=c(-15,15))
lines(mod_pre,
      col='orange',
      lwd=2)
lines(mod_pre_limit,
      col='grey',
      lwd=1,lty=3)
lines(mod_pre_low,
      col='grey',
      lwd=1,lty=3)
legend(x=xy[2],y=xy[4]-yinch(0.8),
       legend=c('True',
                'Predicted',
                '95%CI'),
       xpd=TRUE,
       bty='n',
       col=c('steelblue',
       	'orange',
       	'grey'),
       lty=c(1,1,1),
       lwd=c(2,2,1),
       title = '')
##3.36
# Plot overall Cpg
plot.ts(cpg,
        col='steelblue',
        lwd=2,
        ylab='Retail Price',
        xlab='Time')
# logged Transformation
logg_cpg <- log(cpg)
# Fit the logged series with time t
fit_logg <- lm(logg_cpg~time(logg_cpg))
fitted_cpg <- fit_logg$fitted.values
# Plot fitted values with real logged series
plot(logg_cpg,
        col='steelblue',
        lwd=2,
        ylab='Logged Retail Price',
        xlab='Time')
fitted_cpg <- ts(array(fitted_cpg),start = 1980,end=2008,
   frequency = 1)
lines(fitted_cpg,
      col='orange',
      lwd=2)
legend(x=xy[2],y=xy[4]-yinch(0.8),
       legend=c('True',
                'Predicted'),
       xpd=TRUE,
       bty='n',
       col=c('steelblue',
       	'orange'),
       lty=c(1,1),
       lwd=c(2,2),
       title = '')
# Residual Analysis of fit model
plot_fit(fit_logg,'Residuals of Linear Model')
fit_residuals <- fit_logg$residuals 
# Refit the model with autocorrelated error and Output the LaTex code
fit_all_arima <- sarima(logg_cpg,1,0,0,xreg = time(logg_cpg))
xtable(fit_all_arima$ttable,digits = 3)
# Plot the model fit with real series
plot.ts(logg_cpg,
        col='steelblue',
        lwd=2,
        ylab='Logged Retail Price',
        xlab='Time',
        ylim=c(-5,15))
lines(fit_all_arima$fit$residuals+logg_cpg,
      col='orange',
      lwd=2))
legend(x=xy[2],y=xy[4]-yinch(0.8),
       legend=c('Real',
                'Predicted \nwith autocorrelated'),
       xpd=TRUE,
       bty='n',
       col=c('steelblue',
       	'orange'),
       lty=c(1,1),
       lwd=c(2,2),
       title = '')
