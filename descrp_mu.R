 
 
library(ggplot2)
library(xts)
library(reshape2) # this is required for the melt function below
library(lmtest)
library(base)
library(dplyr)
library(lubridate)
library(dlm)
library(plyr)
###########Plot of the data############
setwd("C:/Users/Muhannad/Desktop/multvariat/Time_series/Tim_fall_2018/project")
data1=read.csv("oil_prod.csv ",stringsAsFactors=FALSE)
data2=read.csv("oil_price.csv ",stringsAsFactors=FALSE)
data1=as.numeric(data1[-c(1,2),-c(1,2)])
data2=data2[-c(1:98),]
data2=head(data2,-5)
data=cbind(as.numeric(data1),as.numeric(data2[,3]))
data=data[seq(1,288),]
##
data01=read.csv("Brent_oil1.csv ",stringsAsFactors=FALSE)
data02=read.csv("Brent_oil2.csv ",stringsAsFactors=FALSE)
data01=data01[-c(1:233),]
data001=rbind(data01,data02)
#
#  Convert to date if not already
data001$Date=as.Date(data001$Date,"%d-%b-%y")

#  Get months
data001$Month <- months(data001$Date)

#  Get years
data001$Year <- format(data001$Date,format="%y")

#  Aggregate 'X2' on months and year and get mean
tt=aggregate( data001$Price ~ Month + Year , data001 , mean )
tt1=tt[order(tt$Month),]
tt$date1=paste(tt$Month, tt$Year, sep="-" ) 

tt$date1=parse_date_time(tt$date1, "my") 

#tt$date1=as.Date(tt$date1,"%B-%y")

tt= tt[order(tt$date1),,drop=FALSE]
data=cbind(data,tt[,3])
#taking the log of the amount of the oil produced.
data[,1]=log(data[,1])
##


data=as.data.frame(data)
colnames(data)=c("oil_prodiction", "oil_prices","Future_price")
 


## to convert to quarterly
monthly <- ts(data,start=c(1994,1),frequency=12)
quarterly <- aggregate(monthly, nfrequency=4, mean)
my_ts=quarterly

#
####################################################################
 

str(my_ts)
opar<-par(mfrow=c(2,1), mar = c(3, 4, 1, 2) + 0.1, oma = c(0, 0, 2, 0))

pdf("descriptive_techniques_plots.pdf", width = 11, height = 11, paper = "a4")
for(i in 1:3) { plot(my_ts[,i],type="l",ylab=colnames(my_ts)[i])
  mtext("A first look at the different variables", outer = TRUE, cex = 1.3) }
save(my_ts,file="my_ts.RData")


###########################################################################


  
lag.plot(my_ts , lags=4)
rev_myts <- ts(matrix(rev(my_ts[,seq(5,1)]),ncol = 5 , dimnames = dimnames(my_ts)), start = start(my_ts)
               , end = end(my_ts) , frequency = frequency(my_ts))
lag.plot(rev_myts , set.lags=4)
pairs(my_ts)
par(mfrow=c(3,2))
for(i in 1:3) monthplot(my_ts[,i],ylab = dimnames(my_ts)[[2]][i])

y <- my_ts[,c(1,2,3)]
par(mfrow = c(3,2) , mar = c(3 , 4.1 , 3 , 2.1))
for(i in 1:3) { acf(y[,i] , main = dimnames(y)[[2]][i])
  pacf(y[,i] , main = dimnames(y)[[2]][i])}
#y<-ts(matrix(as.numeric(y),ncol=3,dimnames=dimnames(y)),
#            start= start(y) , end=end(y), frequency=frequency(y))
#str(y)
acf(y)
x11()
pacf(y)
lag.plot(y,set.lags=4)
revy <- ts(matrix(rev(y[,c(3,2,1)]),ncol=3,dimnames=dimnames(y)),
           start= start(y) , end=end(y), frequency=frequency(y))
lag.plot(revy,set.lags=4)
cor(y)
pairs(y)

par(mfrow=c(3,1))
for(i in 1:3) monthplot(y[,i],ylab = dimnames(y)[[2]][i])

for(i in 1:3){ X11()
  plot(decompose(y[,i]))}
library(GGally)
library(Deducer)
dy <- data.frame(y )# , dimnames =  dimnames(data.frame(y))[[2]])
ggpairs(data.frame(y) ,  color = dimnames(data.frame(y))[[2]] )
ggcorplot(cor.matrix(dy) , dy , var_text_size = 3 , cor_text_limits = c(5,10) , line.method = "loess")


####### Let's specify a training and testing dataset for our 
##################future analysis############################
y_train <- window(y, end = end(y) - c(2, 0))
y_test <- window(y, start = end(y) - c(1, 3)) 

#we only leave two years for trianing data set because we want the model to learn
#from the postcrash behavior of the series.otherwise it would perform poorly

 

data <- list(data_set = y , training_set = y_train , testing_set = y_test)

save(data , file = "project_data.RDATA")

#####
###


###############################################
#SUTSE for 3 dimensional data
###############################################
 
load("project_data.RDATA")
y <-data$data_set
y_train <- data$training_set
y_test <- data$testing_set
str(y_train)
str(y_test)

#Building our three dimensional SUTSE model

mod_one <- dlmModPoly() + dlmModSeas(4) 
mod_sutseM <- dlm(FF = FF(mod_one) %x% diag(3),
                  GG = GG(mod_one) %x% diag(3),
                  V = diag(3), W = diag(x = 0, 15, 15),
                  m0 = rep(0, 15), C0 = diag(x = 1e8, 15, 15))
rm(mod_one)
build_sutseM <- function(psi) {
  U <- diag(exp(0.5 * psi[1:3]), 3, 3)
  U[lower.tri(U)] <- psi[4:6]
  V(mod_sutseM) <- tcrossprod(U) + diag(1e-5, 3, 3)
  diag(U) <- exp(0.5 * psi[7:9])
  U[lower.tri(U)] <- psi[10:12]
  W(mod_sutseM)[4:6, 4:6] <- tcrossprod(U)
  diag(U) <- exp(0.5 * psi[13:15])
  U[lower.tri(U)] <- psi[16:18]
  W(mod_sutseM)[7:9, 7:9] <- tcrossprod(U)
  mod_sutseM
}

### start optimizer from a number of values, randomly selected
set.seed(1252)
n <- 18 # number of parameters to estimate
l <- 18 # number of starting values (assuming l <= n)
Q <- qr.Q(qr(matrix(rnorm(n * n), n, n)))[, 1:l] * rep(runif(l, 1, 10), each = n)

## WARNING: THIS TAKES A FAIRLY LONG TIME
## skip to the next comment if you are short on time
out <- vector("list", l) # for the outputs
system.time(
  for (j in seq.int(l)) {
    cat(j, " of ", l, "\n")
    out[[j]] <- try(dlmMLE(y_train, Q[, j], build_sutseM, control = list(maxit = 500)))
    # if(out[[j]]$convergence==0) browser()
  }
)
smry <- t(sapply(out, function(x) if (!inherits(x, "try-error")) x[c("convergence", "value")]
                 else c(1, NaN)))
out1 <- out[smry[, "convergence"] == 0]
opt <- which.min(sapply(out1, function(x) x$value))
fit_sutse <- out1[[opt]]
#save(fit_sutse,file="Time serSUTSE_3d.RDATA")

load("Time serSUTSE_3d.RDATA")
### fitted model
mod_sutseM <- build_sutseM(fit_sutse$par)
V(mod_sutseM)
cov2cor(V(mod_sutseM))
W(mod_sutseM)[4:6, 4:6]
cov2cor(W(mod_sutseM)[4:6, 4:6])
W(mod_sutseM)[7:9, 7:9]
cov2cor(W(mod_sutseM)[7:9, 7:9])

############################################################
# SUTSE Diagnosis
############################################################
SUTSE_filt <- dlmFilter(y , mod_sutseM )
err <- residuals(SUTSE_filt , sd=F)
windows()
par(mfrow = c(3,2) , mar = c(1,2,1,1))
tsdiag(SUTSE_filt)
lag.plot(err,set.lag=4)

 
tsdiag(SUTSE_filt)
 


lag.plot(err , set.lags = c(1,4))
#############################################
invest=my_ts
investFilt=SUTSE_filt
par(mar=c(3.1,4.1,1.1,2.1), cex=0.5)
require(zoo)
plot(as.zoo(invest), main = "", mar = c(0, 2.1, 0, 1.1),
     oma = c(2.1,0,.1,.1), cex.axis = 0.5, type='o')

sdev <- residuals(investFilt)$sd
lwr <- investFilt$f + qnorm(0.25) * sdev
upr <- investFilt$f - qnorm(0.25) * sdev
 
par(mar=c(2, 3, 1, 0) + 0.1, cex=0.7)
plot(invest[,2], type='o', xlab="", ylab="", main="Kalman filter")
lines(window(investFilt$f[,2], start = start(invest) + c(2,0)), type='o', lty=2, pch=4)
lines(window(lwr[,2], start = start(invest) + c(2,0)), col="darkgrey")
lines(window(upr[,2], start = start(invest) + c(2,0)), col="darkgrey")
legend("topleft", inset = 0.05,
       legend=c("Observed", "One-step-ahead forecast", "50% prediction interval"),
       pch=c(1,4,-1), lty=c(1,2,1), col=c(rep("black", 2), "darkgrey"), bty='n') 


##############################################################################

# A multivariate Dynamic regression model

Oil_price <-y_train[,"oil_prices"]
covariate <- matrix(y[,c(1,3)] , ncol = 2)
dyn_reg <- dlm(FF = matrix(c(1,0,0),nrow=1) , 
               GG = diag(1,3) , V = diag(1) , 
               W = diag(1,3) , JFF = matrix(c(0,1,2) , nrow = 1) , 
               m0 = rep(0, 3), C0 = 1e5 * diag(3), 
               X = covariate )

build_dyn <-function(Psi) {
  V(dyn_reg) <-exp(Psi[1])
  W(dyn_reg) <-diag(exp(Psi[2:4]))
  dyn_reg   }
set.seed(100)
n <- 4 # number of parameters to estimate
l <- 4 # number of starting values (assuming l <= n)
Q <- qr.Q(qr(matrix(rnorm(n * n), n, n)))[, 1:l] * rep(runif(l, 1, 10), each = n)

 
out <- vector("list", l) # for the outputs
system.time(
  for (j in seq.int(l)) {
    cat(j, " of ", l, "\n")
    out[[j]] <- try(dlmMLE(Oil_price, Q[, j], build_dyn, control = list(maxit = 500)))
    # if(out[[j]]$convergence==0) browser()
  }
)
smry <- t(sapply(out, function(x) if (!inherits(x, "try-error")) x[c("convergence", "value")]
                 else c(1, NaN)))
out1 <- out[smry[, "convergence"] == 0]
opt <- which.min(sapply(out1, function(x) x$value))
fit_dyn <- out1[[opt]]
out1[[opt]]$value

dyn_reg <- build_dyn(fit_dyn$par)
 

 
V(dyn_reg)
W(dyn_reg)
 ########################################################################################
############################################
###########################################################################################
 dev.off()

 
dyn_filt <- dlmFilter(y[,"oil_prices"] , dyn_reg)
tsdiag(dyn_filt)
res <- residuals(dyn_filt,sd=F)
lag.plot(res,set.lags=c(1,4))
 
###############################################################
invest=my_ts
investFilt=dyn_filt
par(mar=c(3.1,4.1,1.1,2.1), cex=0.5)
require(zoo)
plot(as.zoo(invest), main = "", mar = c(0, 2.1, 0, 1.1),
     oma = c(2.1,0,.1,.1), cex.axis = 0.5, type='o')

sdev <- residuals(investFilt)$sd
lwr <- investFilt$f + qnorm(0.25) * sdev
upr <- investFilt$f - qnorm(0.25) * sdev
 

par(mar=c(2, 3, 1, 0) + 0.1, cex=0.7)
plot(invest[,2], type='o', xlab="", ylab="", main="Kalman filter")
lines(window(investFilt$f[,2], start = start(invest) + c(2,0)), type='o', lty=2, pch=4)
lines(window(lwr[,2], start = start(invest) + c(2,0)), col="darkgrey")
lines(window(upr[,2], start = start(invest) + c(2,0)), col="darkgrey")
legend("topleft", inset = 0.05,
       legend=c("Observed", "One-step-ahead forecast", "50% prediction interval"),
       pch=c(1,4,-1), lty=c(1,2,1), col=c(rep("black", 2), "darkgrey"), bty='n') 
  

 
reg_mod <- dlm(FF = matrix(c(1,0,0),nrow=1) , 
               GG = diag(1,3) , V = diag(1) , 
               W = diag(1,3) , JFF = matrix(c(0,1,2) , nrow = 1) , 
               m0 = rep(0, 3), C0 = 1e5 * diag(3), 
               X = covariate )
 

 
sreg_mod <- reg_mod + dlmModSeas(4)
 
build_sreg <-function(Psi) {
  V(sreg_mod) <-exp(Psi[1])
  W(sreg_mod) <-diag(c(exp(Psi[2:5]),0,0))
  sreg_mod   }

load("reg+sea.RDATA")
 

 
sreg_mod <- build_sreg(fit_sreg$par)
V(sreg_mod)
matrix(diag(W(sreg_mod)[1:4,1:4]),ncol=1)
 
 
#Filtering and diagnosis
sreg_filt <- dlmFilter(y[,2], sreg_mod)
tsdiag(sreg_filt)
 
 
## _Forecasting_  
 
x11()
 
nAhead <- nrow(y_test)
unempt <-ts(c(Oil_price , rep(NA , nAhead ) ) , start = start(Oil_price) , frequency = 4)
ar_filt <- dlmFilter(unempt ,sreg_mod )
 

 
temp <- window(y[,"oil_prices"] , start = start(y_test)-c(3,3))
tf <- window(ar_filt$f , start = start(y_test), end=end(y_test))
resid <- residuals(ar_filt )
sd <- window(resid$sd , start = start(y_test))
lower = tf - 0.67 * sd
lower_2 = tf - 1.64 * sd
upper = tf + 0.67 * sd
upper_2 = tf + 1.64 * sd
plot(temp , type = "o" , ylab = "oil_prices", main="Forecasting DLR with seas.")
lines(tf ,type = "o" ,pch = 19, col="green")
lines(lower , lty = 2 , col="blue" )
lines(lower_2  , col="red" )
lines(upper_2  , col="red" )
lines(upper , lty = 2 , col="blue" )
legend("bottomleft",inset=0.05,legend=c("observed","forecast","75 % prediction interval","90 % prediction interval"),
       pch=c(1,19,-1,-1),lty=c(1,1,2,1),col=c("black","green","blue","red"),bty='n')
  
