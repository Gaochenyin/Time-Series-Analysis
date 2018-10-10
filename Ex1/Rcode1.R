install.packages('astsa')
library(astsa)
varve
len1 <- length(varve)

par(mfrow=c(2,1))
#1
plot(varve,main="varve", ylab="",xlab='')
abline(v=len1/2,col='red',lty=2)
v1 <- var(varve[1:len1/2])
v2 <- var(varve[len1/2:len1])
text(len1/4,150,paste('v1=',round(v1,3)))
text(3*len1/4,150,paste('v2=',round(v2,3)))
#2
plot(log(varve),main="log(varve)", ylab="")
abline(v=len1/2,col='red',lty=2)
v1_ln <- var(log(varve[1:len1/2]))
v2_ln <- var(log(varve[len1/2:len1]))
text(len1/4,4.8,paste('v1\'=',round(v1_ln,3)))
text(3*len1/4,4.8,paste('v2\'=',round(v2_ln,3)))

par(mfrow=c(2,1))
hist(varve,col='blue',border='white',xlab='')
hist(log(varve),col='blue',border='white',xlab='')
# correlation testing
#1
q1 <- qqnorm(varve)
#2
q2 <- qqnorm(log(varve))
cor(q1$x,q1$y)
cor(q2$x,q2$y)

y <- log(varve)

# log-transformed
plot(y)
w1 <- rep(1,100)/100
y_w1 <- filter(y,sides = 2,filter = w1)
# moving average
lines(y_w1,col='red',lwd=2)
par(fig = c(0, 0.35, .5, 1), new = TRUE) # the insert
nwgts = c(rep(0,20), w1, rep(0,20))
# kenerl function
plot(nwgts, type="l", xaxt='n', ylim=c(-0.01,0.03),yaxt='n', ann=FALSE)

acf(y,main='log(varve)')

y_d <- diff(y)
par(mfrow=c(2,1))
# diff of y
plot(y_d,ylab=expression(paste('diff( ',y[t],')')))
# acf of y
acf(y_d,-5,main='')


# 2.10
plot(oil,type='l',ylim=c(10,400),ylab='Oil&Gas')
lines(gas,col='blue')
#(c)
par(mai=c(0.5,0.5,0.5,1.5))
plot(diff(log(oil)),
     main=expression(paste(nabla,'Log of Oil&Gas')),
     ylab='Diff',
     lwd=2,
     ylim=c(-0.3,0.4))
lines(diff(log(gas)),
      col='blue',
      lwd=2)
xy <- par('usr')
legend(x=xy[2],y=xy[4]-yinch(0.8),legend=c('Oil',
                                           'Gas'),xpd=TRUE,
       bty='n',col=c('black','blue'),lty=1,lwd=2,title = 'Variables'
)
par(mfrow=c(1,2))
acf(diff(log(oil)),main=expression(paste('ACF of ',nabla,'Log Oil')))
acf(diff(log(gas)),main=expression(paste('ACF of ',nabla,'Log Gas')))    
ccf(diff(log(oil)),diff(log(gas)),
    main=expression(paste('CCF of ',nabla,' Log Oil&Gas')),
    ylab='CCF')

#

transformed_gas <- diff(log(gas))
transformed_oil <- diff(log(oil))
par(mfrow=c(1,3))
lag2.plot(transformed_gas,transformed_oil,3)

#


indi <- ifelse(transformed_oil < 0, 0, 1)
mess <- ts.intersect(transformed_gas, transformed_oil, 
                     poilL = lag(transformed_oil,-1), indi)
summary(fit <- lm(transformed_gas~ transformed_oil + poilL + indi, data=mess))
#
ts.plot(fit$residuals/sd(fit$residuals),
        ylab='Model residuals',
        xlab='')
abline(h=3,col='red',lty=2)
abline(h=-3,col='red',lty=2)
anova(fit)
