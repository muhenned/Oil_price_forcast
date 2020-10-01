#####################################
#Descriptive plots and analysis
#of the data
#####################################
library(ggplot2)
library(xts)
library(reshape2) # this is required for the melt function below
library(lmtest)
library(base)
library(dplyr)
library(lubridate)
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
 ##
  
data=as.data.frame(data)
colnames(data)=c("oil_prodiction", "oil_prices","Future_price")
my_ts<-ts(data,start=c(1994,1),freq=12)


## to convert to quarterly
monthly <- ts(data,start=c(1994,1),frequency=12)
quarterly <- aggregate(monthly, nfrequency=4, mean)
my_ts=quarterly

##
opar<-par(mfrow=c(2,1), mar = c(3, 4, 1, 2) + 0.1, oma = c(0, 0, 2, 0))
for(i in 1:3) { plot(my_ts[,i],type="l",ylab=colnames(my_ts)[i])
  mtext(" Oil Price and production", outer = TRUE, cex = 1.3) }
lines(lowess(my_ts[,2]), lwd=2)

plot(data, type = 'o', pch = 20, cex = 0.75, ylab = "ppm", main = "CO2 concentration at Mauna Loa")
plot(stl(my_ts[,2], s.window = "periodic"))
 acf(my_ts[,2])
 
 
opar<-par(mfrow=c(3,2), mar = c(3, 4, 1, 2) + 0.1, oma = c(0, 0, 2, 0))
my_ts<-ts(matrix(as.numeric(my_ts),ncol=ncol(my_ts),dimnames = dimnames(my_ts)),
          start= start(my_ts) , end = end(my_ts) , frequency = frequency(my_ts))


for(i in 1:5) { plot(my_ts[,i],type="l",ylab=colnames(my_ts)[i])
  mtext("A first look at the different variables", outer = TRUE, cex = 1.3) }

plot(globtemp, type = 'o', pch = 20, ylab = ~degree~C,
     main = "Global mean land-ocean temperature\nDeparture from 1951-1980 average")
lines(lowess(globtemp, f = 0.3), lwd = 2, col = "hotpink") ## add smoothed trend
########################################
###########################################
#########################################
########################################
###########################################
#########################################
 
#SUTSE for 3 dimensional data
###############################################
#setwd("C:/Users/admin/Desktop/STAN")
load("project_data.RDATA")
y <-data$data_set
y_train <- data$training_set
y_test <- data$testing_set
#str(y_train)
#str(y_test)
library(dlm)
```
\ 
&nbsp;
```{r SUTSE_build, echo=T  }

#Building our three dimensional SUTSE model

mod_one <- dlmModPoly() + dlmModSeas(4) 
mod_sutseM <- dlm(FF = FF(mod_one) %x% diag(3),
                  GG = GG(mod_one) %x% diag(3),
                  V = diag(3), W = diag(x = 0, 15, 15),
                  m0 = rep(0, 15), C0 = diag(x = 1e8, 15, 15))
```
```{r , echo=F}
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
#load("C:/Users/admin/Desktop/STAN/Time series project data and code/SUTSE_3d.RDATA")
#mod_sutseM <- build_sutseM(fit_sutse$par)
```
\ 
&nbsp;
After estimating the different parameters via _Maximum likelihoud_ method using the _training dataset_ for the estimation. To check the validity of the model, first we filter on the whole dataset then we run the diagnostics on the filtered object.

##Diagnostics
```{r figures , echo=T , results='hide' , fig.show='hide'}
load("SUTSE_3d.RDATA")
mod_sutseM <- build_sutseM(fit_sutse$par)

SUTSE_filt <- dlmFilter(y , mod_sutseM )
err <- residuals(SUTSE_filt , sd=F)
windows()
tsdiag(SUTSE_filt)
lag.plot(err,set.lags=4)
pacf(err)
```
<!--##Diagnostics
  ```{r ,fig.align='center' , fig.height=12}
#knitr::include_graphics("Rplot01.png")
```



