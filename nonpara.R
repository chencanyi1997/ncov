library(tidyverse)
library(np)
library(KernSmooth)
library(locfit)
set.seed(7)

data <- read_csv('./Wuhan-2019-nCoV.csv')
N_data <- read_csv('./常住人口.csv',skip = 3)
N_data <- na.omit(N_data)
# data <- read_csv('https://raw.githubusercontent.com/canghailan/Wuhan-2019-nCoV/master/Wuhan-2019-nCoV.csv')
pName <- '浙江省'
data_p <- data[!is.na(data$province) & data$province == pName & is.na(data$city),]
date <- data_p$date
idx <- 1:length(date)
R <- data_p$cured
# determine h
# optimal for normal kernel
# med <- median(idx)
# n <- length(idx)
# sig <- median(abs(idx-med)) / 0.6745
# h <- sig * (4/(3*n))^(1/5)

## fit R
R.fit <-npreg(R~idx, regtype = "ll",
            bwmethod = "cv.aic",
            gradient = TRUE)
plot(idx,R,
     cex=0.3)
lines(idx,fitted(R.fit),col = 'blue')
lines(idx,gradients(R.fit),col='red')

## fit I
I <- data_p$confirmed
I.fit <- npreg(I~idx, regtype = "ll",
               bwmethod = "cv.aic",
               gradient = TRUE)
plot(idx,I, xlab = 'date',
     ylab = 'Infected Number',
     cex = .2)
lines(idx,fitted(I.fit),col = 'blue')
lines(idx,gradients(I.fit),col = 'red')

## fit S
N <- N_data[N_data$地区==pName,]$`2018年`*1e4 # 常住人口数
S <- N - R - I
S.fit <- npreg(S~idx, regtype = "ll",
               bwmethod = "cv.aic",
               gradient = TRUE)
plot(idx,S, xlab = 'date',
     ylab = 'Susceptible Number',
     cex = .2)
lines(idx,fitted(S.fit),col = 'blue')
plot(idx,dS,col = 'red')

## estimate parameters
dI <- gradients(I.fit)
dR <- gradients(R.fit)
dS <- gradients(S.fit)
plot(idx,(dS+dR+dI),cex = .3,xaxt = 'n')
n <- length(idx)
axis(1,at = idx[seq(n) %% 20 ==1], labels = date[seq(n) %% 20 ==1])
Ihat <- fitted(I.fit)
Rhat <- fitted(R.fit)
Shat <- N - Ihat - Rhat
plot(I,dR,main = 'I vs dR/dt')
gamma <-dR/I
plot(idx,gamma,main = 'gamma')

r_beta <- (dI + dR)*N/(Shat*Ihat)

plot(idx,r_beta)
R0 <- r_beta/gamma

plot(tail(idx,5),tail(R0,5),cex = .3)
mean(R0,trim = 0.2)
mean(tail(R0,10))
