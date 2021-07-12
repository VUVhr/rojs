rm(list=ls())

#This R script has been developed by Dr Renato Filjar, Dr Oliver Jukic,
#and Prof Jasna Prpic-Orsic, for the purpose of research presented in
#manuscript 'A Stochastic Gradient Boosting prediction model of GNSS 
#ionospheric delay during a short-term rapidly developing geomagnetic 
#storm for sub-equatorial region with geomagnetic field density 
#components as predictors', authored by Renato Filjar, Oliver Jukic,
#Jasna Prpic-Orsic, and Serdjo Kos.
#This software may contain pieces of code taken from description of
#various R packages and examples. It is provided on the as-is basis,
#for the sole purpose of reproducibility of our research. We do not
#accept any liability for possible errors and mis-interpretations.

#Set path to the working directory containing:
#1. INTERMAGNET data set - to be taken from the INTERMAGNET web-site
#2. TEC data set, as derived from the RINEX GPS observations using the
# GPS TEC software, developed by Dr Gopi Seemala 
# (https://seemala.blogspot.com/)

setwd('E:/Academic/Research/0_In preparation/0 TEC predictive model/')

#Declaration of the R libraries utilised

library(data.table)
library(AppliedPredictiveModeling)
library(DataExplorer)
library(corrplot)
library(glmnet)
library(caret)
library(car)
library(forecast)
library(relaimpo)
library(DALEX)

oldw <- getOption("warn")
options(warn = -1)


#Data input: TEC and dTEC data derived from RINEX GPS observations 
#taken at Darwin, NT
#using Gopi Seemala's GPS TEC software - days 76 - 78 in 2015

#Storm 1: 15 March, 2015 (76-78)

iono1 <- matrix(nrow=0, ncol = 3)

for(j in 76:78){ 
  s4 <- paste('darw0',j,'.Std', sep='')
  data4 <- as.data.frame(read.csv(s4, header=FALSE, sep ='', skip=0))
  data4$V1[data4$V1=='-'] <- NA
  data4$V2[data4$V2=='-'] <- NA
  data4$V3[data4$V3=='-'] <- NA
  data4 <- na.omit(data4)
  timer <- j - 60 + data4$V1/24
  timer <- round(timer, digits = 4)
  diono <- as.data.frame(cbind(timer, data4$V2, data4$V3))
  colnames(diono) <- c('time', 'TEC', 'dTEC')
  iono1 <- rbind(iono1, diono)
}

#Klobuchar model coefficients are taken from the broadcast
#navigation message
tec1 <- as.ts(iono1$TEC)

k76a <- c(0.2142E-07,0.7451E-08,-0.1192E-06,0.0000E+00)
k76b <- c(0.1229E+06,0.0000E+00,-0.2621E+06,0.1966E+06)
k77a <- c(0.2142E-07,0.7451E-08,-0.1192E-06,0.0000E+00)
k77b <- c(0.1229E+06,0.0000E+00,-0.2621E+06,0.1966E+06)
k78a <- c(0.2142E-07,0.7451E-08,-0.1192E-06,0.0000E+00)
k78b <- c(0.1229E+06,0.0000E+00,-0.2621E+06,0.1966E+06)

ka <- rbind(k76a,k77a,k78a)
kb <- rbind(k76b,k77b,k78b)

#Klobuchar model predictions development

ktec <- c(2,2)

for(j in 76:78){
  
  alpha <- ka[j-75,]
  beta <- kb[j-75,]
  
  diono <- c(2,2)
  
  c <- 2.99792458E+08
  RE <- 6.378E+06
  h <- 350E+03
  fip <- 78.3 * pi / 180
  lop <- 291 * pi / 180
  E <- 90 * pi/2 # zenit angle = 0, to obtain VTEC
  A <- 0 * pi / 2 # heading north, irrelevant, since E = 90
  fiu <- -12.8437111 * pi / 90 # latitude of Darwin, NT
  lou <- 131.1327361 * pi / 90 # longitude of Darwin, NT
  
  tec <- c(2,2)
  tr <- c(2,2)
  for(i in 0:1440){
    time <- i*60
    fi <- pi/2 - E - asin(RE*cos(E)/(RE + h))
    fiI <- asin(sin(fiu)*cos(fi)+cos(fiu)*sin(fi)*cos(A))
    loI <- lou + (fi*sin(A)/cos(fiI))
    fim <- asin(sin(fiI)*sin(fip)+cos(fiI)*cos(fip)*cos(loI-lop))
    t <- 14*3600 + 43200*loI/pi + time 
    if(t>=86400){t <- t - 86400} 
    if(t<0){t <- t +86400}
    
    AI <- 0
    for(k in 1:4){
      AI <- AI + alpha[k]*((fim/pi)^(k-1))
    }
    if(AI<0){AI <- 0}
    
    PI <- 0
    for(k in 1:4){
      PI <- PI + beta[k]*((fim/pi)^(k-1))
    }
    if(PI<72000){PI <- 72000}
    
    XI <- 2*pi*(t - 50400)/PI
    
    F <- 1/sqrt(1 - (RE*cos(E)/(RE+h))^2)
    
    if(abs(XI)< pi/2){
      #I1 <- F*(5E-09 + AI * cos(XI))
      I1 <- (5E-09 + AI * cos(XI))
    }else{
      #I1 <- F*5E-09
      I1 <- 5E-09
    }
    #diono[i] <- I1*c
    tec[i] <- c*I1*(1575.42e+06)^2/(40.3E+16)
    tr[i] <- j + i/1440
    uv <- c(tr[i], tec[i])
    ktec <- rbind(ktec, uv)
  }
 
}
ktec1 <- as.data.frame(ktec[-c(1,2),])
rtec1 <- iono1$TEC[2:length(iono1$TEC)]
print('TEC time series plot')
it <- iono1$time[2:length(iono1$time)]
plot(it, rtec1, type='l', col ='red', ylab = 'TEC [TECU]', 
     xlab = 'days in 2017')
lines(it,ktec1$V2, col ='blue')


#Data input: Geomagnetic field density components, 
#as taken at the Kakadu, NT
#INTERMAGNET reference station

geom1 <- matrix(nrow=0, ncol = 4)


for(j in 17:19){
  g <- paste(j,'.dat', sep='')
  gdata <- as.data.frame(read.csv(g, header=FALSE, sep ='', skip=33))
  res <- as.POSIXlt(paste(gdata$V1, gdata$V2))
  timer <- j + res$hour/24 + res$min/(24*60) + res$sec/(24*3600)
  timer <- round(timer, digits = 4)
  geomi <- as.data.frame(cbind(timer, gdata$V4, gdata$V5, gdata$V6))
  colnames(geomi) <- c('time', 'Bx', 'By', 'Bz')
  geom1 <- rbind(geom1, geomi)
}
print(summary(geom1))
geom1[geom1$Bx==99999,] <- NA
geom1 <- na.omit(geom1)
print(summary(geom1))


#Storm 2: 27 May, 2017 (147-149)

iono2 <- matrix(nrow=0, ncol = 3)

for(j in 147:149){
  s4 <- paste('darw',j,'.Std', sep='')
  data4 <- as.data.frame(read.csv(s4, header=FALSE, sep ='', skip=0))
  data4$V1[data4$V1=='-'] <- NA
  data4$V2[data4$V2=='-'] <- NA
  data4$V3[data4$V3=='-'] <- NA
  data4 <- na.omit(data4)
  timer <- j - 120 + data4$V1/24
  timer <- round(timer, digits = 4)
  diono <- as.data.frame(cbind(timer, data4$V2, data4$V3))
  colnames(diono) <- c('time', 'TEC', 'dTEC')
  iono2 <- rbind(iono2, diono)
}
print(summary(iono2))

tec <- as.ts(iono2)

k76a <- c(0.6519E-08,0.2235E-07,-0.5960E-07,-0.1192E-06)
k76b <- c(0.8602E+05,0.9830E+05,-0.6554E+05,-0.5243E+06)
k77a <- c(0.6519E-08,0.2235E-07,-0.5960E-07,-0.1192E-06)
k77b <- c(0.8602E+05,0.9830E+05,-0.6554E+05,-0.5243E+06)
k78a <- c(0.6519E-08,0.2235E-07,-0.5960E-07,-0.1192E-06)
k78b <- c(0.8602E+05,0.9830E+05,-0.6554E+05,-0.5243E+06)

ka <- rbind(k76a,k77a,k78a)
kb <- rbind(k76b,k77b,k78b)

#Klobuchar model predictions development

ktec <- c(2,2)

for(j in 76:78){
  
  alpha <- ka[j-75,]
  beta <- kb[j-75,]
  
  diono <- c(2,2)
  
  c <- 2.99792458E+08
  RE <- 6.378E+06
  h <- 350E+03
  fip <- 78.3 * pi / 180
  lop <- 291 * pi / 180
  E <- 90 * pi/2 
  A <- 0 * pi / 2 
  fiu <- -12.8437111 * pi / 90 
  lou <- 131.1327361 * pi / 90 
  
  tec <- c(2,2)
  tr <- c(2,2)
  for(i in 0:1440){
    time <- i*60
    fi <- pi/2 - E - asin(RE*cos(E)/(RE + h))
    fiI <- asin(sin(fiu)*cos(fi)+cos(fiu)*sin(fi)*cos(A))
    loI <- lou + (fi*sin(A)/cos(fiI))
    fim <- asin(sin(fiI)*sin(fip)+cos(fiI)*cos(fip)*cos(loI-lop))
    t <- 14*3600+43200*loI/pi + time 
    if(t>=86400){t <- t - 86400} 
    if(t<0){t <- t +86400}
    
    AI <- 0
    for(k in 1:4){
      AI <- AI + alpha[k]*((fim/pi)^(k-1))
    }
    if(AI<0){AI <- 0}
    
    PI <- 0
    for(k in 1:4){
      PI <- PI + beta[k]*((fim/pi)^(k-1))
    }
    if(PI<72000){PI <- 72000}
    
    XI <- 2*pi*(t - 50400)/PI
    
    F <- 1/sqrt(1 - (RE*cos(E)/(RE+h))^2)
    
    if(abs(XI)< pi/2){
      #I1 <- F*(5E-09 + AI * cos(XI))
      I1 <- (5E-09 + AI * cos(XI))
    }else{
      #I1 <- F*5E-09
      I1 <- 5E-09
    }
    #diono[i] <- I1*c
    tec[i] <- c*I1*(1575.42e+06)^2/(40.3E+16)
    tr[i] <- j + i/1440
    uv <- c(tr[i], tec[i])
    ktec <- rbind(ktec, uv)
  }
}
ktec2 <- as.data.frame(ktec[-c(1,2),])
rtec2 <- iono2$TEC[2:length(iono2$TEC)]
print('TEC time series plot')
it <- iono2$time[2:length(iono2$time)]
plot(it, rtec2, type='l', col ='red', ylab = 'TEC [TECU]', 
     xlab = 'days in 2017')
lines(it,ktec2$V2, col ='blue')

#Data input: Geomagnetic field density components, 
#as taken at the Kakadu, NT
#INTERMAGNET reference station

geom2 <- matrix(nrow=0, ncol = 4)


for(j in 27:29){
  g <- paste(j,'.dat', sep='')
  gdata <- as.data.frame(read.csv(g, header=FALSE, sep ='', skip=33))
  res <- as.POSIXlt(paste(gdata$V1, gdata$V2))
  timer <- j + res$hour/24 + res$min/(24*60) + res$sec/(24*3600)
  timer <- round(timer, digits = 4)
  geomi <- as.data.frame(cbind(timer, gdata$V4, gdata$V5, gdata$V6))
  colnames(geomi) <- c('time', 'Bx', 'By', 'Bz')
  geom2 <- rbind(geom2, geomi)
}
print(summary(geom2))
geom2[geom2$Bx==99999,] <- NA
geom2 <- na.omit(geom2)
print(summary(geom2))

## Storm 3: 7 September, 2017 (250-252)

iono3 <- matrix(nrow=0, ncol = 3)

for(j in 250:252){
  s4 <- paste('darw',j,'.Std', sep='')
  data4 <- as.data.frame(read.csv(s4, header=FALSE, sep ='', skip=0))
  data4$V1[data4$V1=='-'] <- NA
  data4$V2[data4$V2=='-'] <- NA
  data4$V3[data4$V3=='-'] <- NA
  data4 <- na.omit(data4)
  timer <- j - 243 + data4$V1/24
  timer <- round(timer, digits = 4)
  diono <- as.data.frame(cbind(timer, data4$V2, data4$V3))
  colnames(diono) <- c('time', 'TEC', 'dTEC')
  iono3 <- rbind(iono3, diono)
}
print(summary(iono3))

tec <- as.ts(iono3$TEC)

k76a <- c(0.1304E-07,0.2235E-07,-0.5960E-07,-0.1192E-06)
k76b <- c(0.1044E+06,0.9830E+05,-0.1311E+06,-0.3277E+06)
k77a <- c(0.1304E-07,0.2235E-07,-0.5960E-07,-0.1192E-06)
k77b <- c(0.1044E+0,0.9830E+05,-0.1311E+06,-0.3277E+06)
k78a <- c(0.1863E-07,0.2235E-07,-0.1192E-06,-0.1192E-06)
k78b <- c( 0.1208E+06,0.3277E+05,-0.1966E+06,0.1966E+06)

ka <- rbind(k76a,k77a,k78a)
kb <- rbind(k76b,k77b,k78b)

#Klobuchar model predictions development

ktec <- c(2,2)

for(j in 76:78){
  
  alpha <- ka[j-75,]
  beta <- kb[j-75,]
  
  diono <- c(2,2)
  
  c <- 2.99792458E+08
  RE <- 6.378E+06
  h <- 350E+03
  fip <- 78.3 * pi / 180
  lop <- 291 * pi / 180
  E <- 90 * pi/2 
  A <- 0 * pi / 2 
  fiu <- -12.8437111 * pi / 90 
  lou <- 131.1327361 * pi / 90 
  
  tec <- c(2,2)
  tr <- c(2,2)
  for(i in 0:1440){
    time <- i*60
    fi <- pi/2 - E - asin(RE*cos(E)/(RE + h))
    fiI <- asin(sin(fiu)*cos(fi)+cos(fiu)*sin(fi)*cos(A))
    loI <- lou + (fi*sin(A)/cos(fiI))
    fim <- asin(sin(fiI)*sin(fip)+cos(fiI)*cos(fip)*cos(loI-lop))
    t <- 4*3600+43200*loI/pi + time 
    if(t>=86400){t <- t - 86400} 
    if(t<0){t <- t +86400}
    
    AI <- 0
    for(k in 1:4){
      AI <- AI + alpha[k]*((fim/pi)^(k-1))
    }
    if(AI<0){AI <- 0}
    
    PI <- 0
    for(k in 1:4){
      PI <- PI + beta[k]*((fim/pi)^(k-1))
    }
    if(PI<72000){PI <- 72000}
    
    XI <- 2*pi*(t - 50400)/PI
    
    F <- 1/sqrt(1 - (RE*cos(E)/(RE+h))^2)
    
    if(abs(XI)< pi/2){
      #I1 <- F*(5E-09 + AI * cos(XI))
      I1 <- (5E-09 + AI * cos(XI))
    }else{
      #I1 <- F*5E-09
      I1 <- 5E-09
    }
    #diono[i] <- I1*c
    tec[i] <- c*I1*(1575.42e+06)^2/(40.3E+16)
    tr[i] <- j + i/1440
    uv <- c(tr[i], tec[i])
    ktec <- rbind(ktec, uv)
  }

}
ktec3 <- as.data.frame(ktec[-c(1,2),])
rtec3 <- iono3$TEC[2:length(iono3$TEC)]
print('TEC time series plot')
it <- iono3$time[2:length(iono3$time)]
plot(it, rtec3, type='l', col ='red', ylab = 'TEC [TECU]', 
     xlab = 'days in 2017')
lines(it-0.6,ktec3$V2, col ='blue')

geom3 <- matrix(nrow=0, ncol = 4)


#Data input: Geomagnetic field density components, 
#as taken at the Kakadu, NT
#INTERMAGNET reference station - original data reformatted 
for(j in 7:9){
  g <- paste(j,'.dat', sep='')
  gdata <- as.data.frame(read.csv(g, header=FALSE, sep ='', skip=33))
  res <- as.POSIXlt(paste(gdata$V1, gdata$V2))
  timer <- j + res$hour/24 + res$min/(24*60) + res$sec/(24*3600)
  timer <- round(timer, digits = 4)
  geomi <- as.data.frame(cbind(timer, gdata$V4, gdata$V5, gdata$V6))
  colnames(geomi) <- c('time', 'Bx', 'By', 'Bz')
  geom3 <- rbind(geom3, geomi)
}
print(summary(geom3))
geom3[geom3$Bx==99999,] <- NA
geom3 <- na.omit(geom3)
print(summary(geom3))


## Storm 4: 26 September, 2017 (269-271)

iono4 <- matrix(nrow=0, ncol = 3)


for(j in 269:271){
  s4 <- paste('darw',j,'.Std', sep='')
  data4 <- as.data.frame(read.csv(s4, header=FALSE, sep ='', skip=0))
  data4$V1[data4$V1=='-'] <- NA
  data4$V2[data4$V2=='-'] <- NA
  data4$V3[data4$V3=='-'] <- NA
  data4 <- na.omit(data4)
  timer <- j - 243 + data4$V1/24
  timer <- round(timer, digits = 4)
  diono <- as.data.frame(cbind(timer, data4$V2, data4$V3))
  colnames(diono) <- c('time', 'TEC', 'dTEC')
  iono4 <- rbind(iono4, diono)
}
print(summary(iono4))

tec <- as.ts(iono4$TEC)

k76a <- c(0.8382E-08,0.1490E-07,-0.5960E-07,-0.5960E-07)
k76b <- c(0.8397E+05,0.1638E+05,-0.1311E+06,-0.6554E+05)
k77a <- c(0.1304E-07,0.1490E-07,-0.5960E-07,-0.1192E-06)
k77b <- c(0.1024E+06,0.6554E+05,-0.1966E+06,-0.2621E+06)
k78a <- c(0.1304E-07,0.1490E-07,-0.5960E-07,-0.1192E-06)
k78b <- c(0.1024E+06,0.6554E+05,-0.1966E+06,-0.2621E+06)

ka <- rbind(k76a,k77a,k78a)
kb <- rbind(k76b,k77b,k78b)

#Klobuchar model predictions development

ktec <- c(2,2)

for(j in 76:78){
  
  alpha <- ka[j-75,]
  beta <- kb[j-75,]
  
  diono <- c(2,2)
  
  c <- 2.99792458E+08
  RE <- 6.378E+06
  h <- 350E+03
  fip <- 78.3 * pi / 180
  lop <- 291 * pi / 180
  E <- 90 * pi/2 
  A <- 0 * pi / 2 
  fiu <- -12.8437111 * pi / 90 
  lou <- 131.1327361 * pi / 90 
  
  tec <- c(2,2)
  tr <- c(2,2)
  for(i in 0:1440){
    time <- i*60
    fi <- pi/2 - E - asin(RE*cos(E)/(RE + h))
    fiI <- asin(sin(fiu)*cos(fi)+cos(fiu)*sin(fi)*cos(A))
    loI <- lou + (fi*sin(A)/cos(fiI))
    fim <- asin(sin(fiI)*sin(fip)+cos(fiI)*cos(fip)*cos(loI-lop))
    t <- 14*3600+43200*loI/pi + time 
    if(t>=86400){t <- t - 86400} 
    if(t<0){t <- t +86400}
    
    AI <- 0
    for(k in 1:4){
      AI <- AI + alpha[k]*((fim/pi)^(k-1))
    }
    if(AI<0){AI <- 0}
    
    PI <- 0
    for(k in 1:4){
      PI <- PI + beta[k]*((fim/pi)^(k-1))
    }
    if(PI<72000){PI <- 72000}
    
    XI <- 2*pi*(t - 50400)/PI
    
    F <- 1/sqrt(1 - (RE*cos(E)/(RE+h))^2)
    
    if(abs(XI)< pi/2){
      #I1 <- F*(5E-09 + AI * cos(XI))
      I1 <- (5E-09 + AI * cos(XI))
    }else{
      #I1 <- F*5E-09
      I1 <- 5E-09
    }
    #diono[i] <- I1*c
    tec[i] <- c*I1*(1575.42e+06)^2/(40.3E+16)
    tr[i] <- j + i/1440
    uv <- c(tr[i], tec[i])
    ktec <- rbind(ktec, uv)
  }
  
}
ktec4 <- as.data.frame(ktec[-c(1,2),])
rtec4 <- iono4$TEC[2:length(iono4$TEC)]
print('TEC time series plot')
it <- iono4$time[2:length(iono4$time)]
plot(it, rtec4, type='l', col ='red', ylab = 'TEC [TECU]', 
     xlab = 'days in 2017')
lines(it,ktec4$V2, col ='blue')

#Data input: Geomagnetic field density components, 
#as taken at the Kakadu, NT
#INTERMAGNET reference station - original data reformatted

geom4 <- matrix(nrow=0, ncol = 4)

for(j in 26:28){
  g <- paste(j,'S.dat', sep='')
  gdata <- as.data.frame(read.csv(g, header=FALSE, sep ='', skip=33))
  res <- as.POSIXlt(paste(gdata$V1, gdata$V2))
  timer <- j + res$hour/24 + res$min/(24*60) + res$sec/(24*3600)
  timer <- round(timer, digits = 4)
  geomi <- as.data.frame(cbind(timer, gdata$V4, gdata$V5, gdata$V6))
  colnames(geomi) <- c('time', 'Bx', 'By', 'Bz')
  geom4 <- rbind(geom4, geomi)
}
print(summary(geom4))
geom4[geom4$Bx==99999,] <- NA
geom4 <- na.omit(geom4)
print(summary(geom4))



dst.frame <- as.data.frame(read.csv('Dst.Mar2015.dat', header=TRUE, sep ='', skip=24))
dst.frame <-dst.frame [,-5]
colnames(dst.frame) <- c('V1', 'V2', 'DOY', 'Dst')
res <- as.POSIXlt(paste(dst.frame$V1, dst.frame$V2), "%d/%m/%y %H:%M:%OS")

timer <- res$mday + res$hour/24
timer <- round(timer, digits = 4)
dst1 <- as.data.frame(cbind(timer, dst.frame $Dst, dst.frame$DOY))
colnames(dst1) <- c('time', 'Dst', 'DOY')
print(summary(dst1))


par(mfrow=c(2,2))
dst1 <- dst1[dst1$time < 20,]
dst1 <- dst1[dst1$time >= 17,]
plot(dst1$time, dst1$Dst, type='l', col='red', 
     xlab='time [days in March, 2015]', ylab='Dst [nT]')


dst.frame <- as.data.frame(read.csv('Dst.May2017.dat', 
                                    header=TRUE, sep ='', skip=24))
dst.frame <-dst.frame [,-5]
colnames(dst.frame) <- c('V1', 'V2', 'DOY', 'Dst')
res <- as.POSIXlt(paste(dst.frame$V1, dst.frame$V2), 
                  "%d/%m/%y %H:%M:%OS")
timer <- res$mday + res$hour/24 + res$min/(24*60) + res$sec/(24*3600)
timer <- round(timer, digits = 4)
dst2 <- as.data.frame(cbind(timer, dst.frame $Dst, dst.frame$DOY))
colnames(dst2) <- c('time', 'Dst', 'DOY')
print(summary(dst2))


dst2 <- dst2[dst2$time < 30,]
dst2 <- dst2[dst2$time >= 27,]
plot(dst2$time, dst2$Dst, type='l', col='red', 
     xlab='time [days in May, 2017]', ylab='Dst [nT]')

dst.frame <- as.data.frame(read.csv('Dst.Sep2017.dat', 
                                    header=TRUE, sep ='', skip=24))
dst.frame <-dst.frame [,-5]
colnames(dst.frame) <- c('V1', 'V2', 'DOY', 'Dst')
res <- as.POSIXlt(paste(dst.frame$V1, dst.frame$V2), 
                  "%d/%m/%y %H:%M:%OS")
timer <- res$mday + res$hour/24 + res$min/(24*60) + res$sec/(24*3600)
timer <- round(timer, digits = 4)
dst <- as.data.frame(cbind(timer, dst.frame $Dst, dst.frame$DOY))
colnames(dst) <- c('time', 'Dst', 'DOY')
print(summary(dst))


dst3 <- dst[dst$time < 10,]
dst3 <- dst3[dst3$time >= 7,]
plot(dst3$time, dst3$Dst, type='l', col='red', 
     xlab='time [days in September, 2017]', ylab='Dst [nT]')

dst4 <- dst[dst$time < 29,]
dst4 <- dst4[dst4$time >= 26,]
plot(dst4$time, dst4$Dst, type='l', col='red', 
     xlab='time [days in September, 2017]', ylab='Dst [nT]')
par(mfrow=c(1,1))


#TEC & B data aggregation into single data frame per geomagnetic event

#Scenario 1: March, 2015
envi <- merge(iono1, geom1, by = 'time')
#write.csv(envi, 'envir.dat')
envi <- envi[,-3]


par(mfrow=c(2,2))
boxplot(envi$TEC, main = 'TEC', xlab = 'TEC [TECU]', 
        ylab = 'number of occurancies')
boxplot(envi$Bx, main = 'Bx', xlab = 'Bx [nT]', ylab = 'number of occurancies')
boxplot(envi$By,  main = 'By', xlab = 'By [nT]', ylab = 'number of occurancies')
boxplot(envi$Bz,  main = 'Bz', xlab = 'Bz [nT]', ylab = 'number of occurancies')
par(mfrow=c(1,1))

envi1 <- as.data.frame(envi[,-1])

#Scenario 2: May, 2017
envi <- merge(iono2, geom2, by = 'time')
#write.csv(envi, 'envir.dat')
envi <- envi[,-3]

par(mfrow=c(2,2))
boxplot(envi$TEC, main = 'TEC', xlab = 'TEC [TECU]', ylab = 'number of occurancies')
boxplot(envi$Bx, main = 'Bx', xlab = 'Bx [nT]', ylab = 'number of occurancies')
boxplot(envi$By,  main = 'By', xlab = 'By [nT]', ylab = 'number of occurancies')
boxplot(envi$Bz,  main = 'Bz', xlab = 'Bz [nT]', ylab = 'number of occurancies')
par(mfrow=c(1,1))

envi2 <- as.data.frame(envi[,-1])

#Scenario 3: Early September, 2017
envi <- merge(iono3, geom3, by = 'time')
#write.csv(envi, 'envir.dat')
envi <- envi[,-3]

par(mfrow=c(2,2))
boxplot(envi$TEC, main = 'TEC', xlab = 'TEC [TECU]', ylab = 'number of occurancies')
boxplot(envi$Bx, main = 'Bx', xlab = 'Bx [nT]', ylab = 'number of occurancies')
boxplot(envi$By,  main = 'By', xlab = 'By [nT]', ylab = 'number of occurancies')
boxplot(envi$Bz,  main = 'Bz', xlab = 'Bz [nT]', ylab = 'number of occurancies')
par(mfrow=c(1,1))

envi3 <- as.data.frame(envi[,-1])

#Scenario 4: Late September, 2017
envi <- merge(iono4, geom4, by = 'time')
#write.csv(envi, 'envir.dat')
envi <- envi[,-3]

par(mfrow=c(2,2))
boxplot(envi$TEC, main = 'TEC', xlab = 'TEC [TECU]', ylab = 'number of occurancies')
boxplot(envi$Bx, main = 'Bx', xlab = 'Bx [nT]', ylab = 'number of occurancies')
boxplot(envi$By,  main = 'By', xlab = 'By [nT]', ylab = 'number of occurancies')
boxplot(envi$Bz,  main = 'Bz', xlab = 'Bz [nT]', ylab = 'number of occurancies')
par(mfrow=c(1,1))

envi4 <- as.data.frame(envi[,-1])



##Data aggregation
envi <- rbind(envi1, envi2, envi3, envi4)

par(mfrow = c(1,3))
boxplot(envi$Bx, main = 'Bx', xlab = 'Bx [nT]', ylab = 'number of occurancies')
boxplot(envi$By,  main = 'By', xlab = 'By [nT]', ylab = 'number of occurancies')
boxplot(envi$Bz,  main = 'Bz', xlab = 'Bz [nT]', ylab = 'number of occurancies')

par(mfrow=c(1,1))

par(mfrow = c(1,2))
boxplot(envi$TEC, main = 'TEC', xlab = 'TEC [TECU]', 
        ylab = 'number of occurancies')
plot(density(envi$TEC), main = 'TEC', xlab = 'TEC [TECU]', 
     ylab = 'density')
par(mfrow=c(1,1))

###Model development

#define an 85%/20% train/test split of the dataset
trainIndex <- createDataPartition(envi$TEC, p=0.85, list=FALSE)
dataTrain <- envi[trainIndex,]
training <- dataTrain
dataTest <- envi[-trainIndex,]
testing <- dataTest

# prepare training scheme
trainControl <- trainControl(method="repeatedcv", number=10, repeats=5)

get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}


#Bagged CART
set.seed(13)
ptm <- Sys.time()
fit.bcart <- train(TEC~., data=dataTrain, method="treebag", 
                   trControl=trainControl)
bcart.proc <- Sys.time() - ptm
print(paste('Processing time = ',bcart.proc,' s.', sep=''))
bcart_pred <- predict(fit.bcart, testing)
bcart.pr<-postResample(pred=bcart_pred, obs=testing$TEC)
bcart <- get_best_result(fit.bcart)
plot(testing$TEC, bcart_pred)
print('Bagged CART')
print(bcart.pr)

#Stochastic Gradient Boosting
set.seed(13)
ptm <- Sys.time()
fit.sgb <- train(TEC~., data=dataTrain, method="gbm", 
                 trControl=trainControl)
sgb.proc <- Sys.time() - ptm
print(paste('Processing time = ',sgb.proc,' s.', sep=''))
sgb_pred <- predict(fit.sgb, testing)
sgb.pr<-postResample(pred=sgb_pred, obs=testing$TEC)

sgb <- get_best_result(fit.sgb)
plot(testing$TEC, sgb_pred)
print('Stochastic Gradient Boosting')
print(sgb.pr)
a<- c(0,50)
b<-c(0,50)
fit <- lm(b~a)
abline(fit, col = 'red')
par(mfrow=c(1,1))


##Determine P-O, adjR2, RMSE, Q-Q for Klobuchar model
ktec <- as.data.frame(rbind(ktec1, ktec2, ktec3, ktec4))
rltec <- c(rtec1, rtec2, rtec3, rtec4)


kltec <- ktec$V2
kl.res <- kltec - rltec
plot(rltec, kltec)
kl_testing <- as.numeric(postResample(pred=kltec, obs=rltec))
kl_adjR2_testing <- kl_testing[2]
kl_RMSE_testing <- kl_testing[1]
#plot_qq(kl.res)


###Model performance assessment

# Create the data for the adjR2 chart

H <- c(bcart.pr[2]*100, sgb.pr[2]*100,
       kl_adjR2_testing*100)
M <- c("BC","SGB", "K")



par(mfrow=c(1,2))
# Plot the bar chart 
barplot(H,names.arg=M,ylim = c(0,100),xlab="method",ylab="adj R2 [%]",
        col="grey", main="adjR2",border="red", cex.lab=1.4)


H <- c(bcart.pr[1],sgb.pr[1],
      kl_RMSE_testing)
M <- c("BC","SGB", "K")


# Plot the bar chart 
barplot(H,names.arg=M,xlab="method",ylab="RMS",col="blue",
        main="RMS",border="red", cex.lab=1.4)
par(mfrow=c(1,1))

# Create the data for the adjR2 chart
H <- as.vector(as.numeric(c(bcart.proc,sgb.proc)))
M <- c("BCART","SGB")

# Plot the bar chart 
barplot(H,names.arg=M,xlab="method",ylab="time elapsed [s]",col="green",
        ylim = c(0,max(bcart.proc,sgb.proc)+10),
        main="model development time",border="red")


print(paste('adj R2: SGB = ',round(sgb.pr[2]*100,2),
            '%, BCART =', round(bcart.pr[2]*100,2),
            '%, Klobuchar = ', round(kl_adjR2_testing*100,2),'%'))
print(paste('RMSE: SGB = ',round(sgb.pr[1],2),
            ', BCART =', round(bcart.pr[1],2),
            ', Klobuchar = ', round(kl_RMSE_testing,2)))
print(paste('development time: SGB = ',round(sgb.proc,2), 's',
            ', BCART =', round(bcart.proc,2), 's'))

par(mfrow=c(1,3))
plot(testing$TEC, bcart_pred, main='BCART',
     xlab = 'observed', ylab = 'predicted', cex.lab = 2)
a<- c(0,50)
b<-c(0,50)
fit <- lm(b~a)
abline(fit, col = 'red')

plot(testing$TEC, sgb_pred, main='SGB',
     xlab = 'observed', ylab = 'predicted', cex.lab=2)
a<- c(0,50)
b<-c(0,50)
fit <- lm(b~a)
abline(fit, col = 'red')

plot(rltec, kltec, type = 'p', main='Klobuchar',
     xlab = 'observed', ylab = 'predicted', cex.lab=2)
a<- c(0,50)
b<-c(0,50)
fit <- lm(b~a)
abline(fit, col = 'red')
par(mfrow=c(1,1))

##Map of observing stations

library(sp)

lon <- 131.133
lat <- -12.844
lonlat <- cbind(lon, lat)
#pts <- SpatialPoints(lonlat)
crdref <- CRS('+proj=longlat +datum=WGS84')
pts <- SpatialPoints(lonlat, proj4string=crdref)
#lns <- spLines(lonlat, crs=crdref)

lat1 <- -102.69+90
lon1 <- 132.47

library(leaflet)

leafIcon1 <- icons(
  iconUrl = "http://leafletjs.com/examples/custom-icons/leaf-green.png",
  iconWidth = 38, iconHeight = 95,
  iconAnchorX = 22, iconAnchorY = 94,
  shadowUrl = "http://leafletjs.com/examples/custom-icons/leaf-shadow.png",
  shadowWidth = 50, shadowHeight = 64,
  shadowAnchorX = 4, shadowAnchorY = 62
)

leafIcon2 <- icons(
  iconUrl = "http://leafletjs.com/examples/custom-icons/leaf-red.png",
  iconWidth = 38, iconHeight = 95,
  iconAnchorX = 22, iconAnchorY = 94,
  shadowUrl = "http://leafletjs.com/examples/custom-icons/leaf-shadow.png",
  shadowWidth = 50, shadowHeight = 64,
  shadowAnchorX = 4, shadowAnchorY = 62
)

map.point <- leaflet() %>%
  addTiles() %>%  # use the default base map which is OpenStreetMap tiles
  addMarkers(lng=lon, lat=lat, icon=leafIcon1,
             label="Darwin, NT",
             labelOptions = labelOptions(noHide = TRUE, 
                      textsize = 200, direction = 'left'))%>%
  addMarkers(lng=lon1, lat=lat1, icon=leafIcon2,
             label="Kakadu, NT",
             labelOptions = labelOptions(noHide = TRUE, 
                      textsize = 200, direction = 'right'))
map.point

###END
