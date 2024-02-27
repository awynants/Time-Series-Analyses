library(readxl)
data <- read_excel(file.choose())
data <- data[,2:3]

#1
library(ggplot2)
library(ggfortify)
library(aTSA)
library(forecast)
N = 300
Ts = 0.2
Fs = 5
y <- ts(data$SystolicBP, frequency=Fs)
ggp1 <- autoplot(y, geom = "line", main = "Blood pressure signal", xlab = "time (s)", ylab = "systolic blood pressure (mm Hg)") + theme_bw()
plot(ggp1)


#1

autoplot(y)
autoplot(diff(y), geom = "line", main = "Differenced blood pressure signal", xlab = "time (s)", ylab = "systolic blood pressure (mm Hg") + theme_bw()
#Mean stationarity does look much better after differentiating, more mean reverting
#Maybe still some periodic volatility clustering left. On the whole somewhat homoscedastic.
#For completeness ADF test even though it's not very powerful
adf.test(y) #ADF test suggests stationary after eliminating drift -> agrees that the series needs differentiation
#So we work with the differenced series
ydiff <- diff(y)


#2

N <- length(y); 
acf(y, lag.max=N/10) #Periodicity with period length like 5 so 1 second. Patients were not measured at rest 
#so related to fast breath maybe? We might see breath as an AR generator close to the unit circle
pacf(y, lag.max=N/10)
#Maybe try seasonal differentiation? 

N <- length(ydiff)
acf(ydiff, lag.max=N/5) #Periodicity with period length like 5 so 1 second. Patients were not measured at rest 
#so related to fast breath maybe? We might see breath as an AR generator close to the unit circle
pacf(ydiff, lag.max=N/5)
#Maybe try seasonal differentiation? 


ysdiff <- diff(ydiff, lag = 5)
N <- length(ysdiff)
acf(ysdiff, lag.max = N/5) #Looks better, autocorrelation drops off now. Can't clearly see seasonality anymore.
#NVM if you plot it with longer lag you see it doesn't drop off. Ok maybe the acf looks better without seasonal diff
pacf(ysdiff, lag.max = N/5) #Some AR periodicity though

#If I use seasonal difference it looks like ARIMA(4, 1, 2/6)(2, 1, 1)[5]
#Without seasonal differences it looks like ARIMA(6, 1, 4/8)(0, 0, ?)[5]

#Seasonal differences is more logical

N2 <- length(ysdiff)
par(mfrow = c(1,1))
acf(ysdiff, lag.max = N/5)
pacf(ysdiff, lag.max = N/5)

N1 <- length(ydiff)
N2 <- length(ysdiff)
par(mfrow = c(2,2))
acf(ydiff, lag.max = N1/5)
acf(ysdiff, lag.max = N2/5)
pacf(ydiff, lag.max = N1/5)
pacf(ysdiff, lag.max = N2/5)

#Look at spectrum, just AR components for now
N <- length(ysdiff)
s_period <- spectrum(ysdiff,method="pgram", fast = FALSE, log = "no")
s_ar <- spectrum(ysdiff,method="ar", type = "b", log = "yes", n.freq = (N/2)+1, order = 17)
Spectral <- data.frame(freq = s_period$freq, Period <- s_period$spec, AR <- s_ar$spec[2:(N/2+1)])
plot(ggplot(Spectral, aes(x=freq)) +  geom_line(aes(y = Period, colour = "periodogram")) + geom_line(aes(y=AR, colour = "AR")))
#AIC selects order 17 AR
#There's a stop band as well around 1 so there must be MA components

AR1 <- ar(ysdiff,aic = FALSE, order = 17, method = c("ols"))
A1 <- c(1,-AR1$ar)
AR_Model1 <-data.frame(Roots_fp = 1/polyroot(A1), Radius = abs(1/polyroot(A1)), Angle = Arg(1/polyroot(A1))/2/pi*Fs)
#17 is too large bc it's giving roots that aren't on the periodogram, overfit while the periodogram is a nonparametric estimate
#The higher the radius the stronger the frequency should be, frequency 1.3 is a bit low down in that sense,
#Radius lower than about 0.8 is hard to identify visually

#So we lower the order, basically trial and error
#Say 8
N <- length(ysdiff)
s_period <- spectrum(ysdiff,method="pgram", fast = FALSE, log = "no")
s_ar <- spectrum(ysdiff,method="ar", type = "b", log = "yes", n.freq = (N/2)+1, order = 8)
Spectral <- data.frame(freq = s_period$freq, Period <- s_period$spec, AR <- s_ar$spec[2:(N/2+1)])
plot(ggplot(Spectral, aes(x=freq)) +  geom_line(aes(y = Period, colour = "periodogram")) + geom_line(aes(y=AR, colour = "AR")))

AR1 <- ar(ysdiff,aic = FALSE, order = 8, method = c("ols"))
A1 <- c(1,-AR1$ar)
AR_Model1 <-data.frame(Roots_fp = 1/polyroot(A1), Radius = abs(1/polyroot(A1)), Angle = Arg(1/polyroot(A1))/2/pi*Fs)
#Still a bit odd that 1.17 is in there, even lower then

#order 4 seems to match the maxima in the periodogram
N <- length(ysdiff)
s_period <- spectrum(ysdiff,method="pgram", fast = FALSE, log = "no")
s_ar <- spectrum(ysdiff,method="ar", type = "b", log = "yes", n.freq = (N/2)+1, order = 4)
Spectral <- data.frame(freq = s_period$freq, Period <- s_period$spec, AR <- s_ar$spec[2:(N/2+1)])
plot(ggplot(Spectral, aes(x=freq)) +  geom_line(aes(y = Period, colour = "periodogram")) + geom_line(aes(y=AR, colour = "AR")))

AR1 <- ar(ysdiff,aic = FALSE, order = 4, method = c("ols"))
A1 <- c(1,-AR1$ar)
AR_Model1 <-data.frame(Roots_fp = 1/polyroot(A1), Radius = abs(1/polyroot(A1)), Angle = Arg(1/polyroot(A1))/2/pi*Fs)
#Three stop bands so 6 MA terms? Or does only the 1 in the middle count idfk

#->ARMA(4, 0, 2)? Or (4, 0, 6) depending on whether you count the edges of the graph


#Let's try without seasonal differencing
#Look at spectrum, just AR components for now
N <- length(ydiff)
s_period <- spectrum(ydiff,method="pgram", fast = FALSE, log = "no")
s_ar <- spectrum(ydiff,method="ar", type = "b", log = "yes", n.freq = (N/2)+1, order = 17)
Spectral <- data.frame(freq = s_period$freq, Period <- s_period$spec, AR <- s_ar$spec[2:(N/2+1)])
plot(ggplot(Spectral, aes(x=freq)) +  geom_line(aes(y = Period, colour = "periodogram")) + geom_line(aes(y=AR, colour = "AR")))
#AIC selects order 17 AR as well, not sure if signal is attenuated

AR1 <- ar(ydiff,aic = FALSE, order = 17, method = c("ols"))
A1 <- c(1,-AR1$ar)
AR_Model1 <-data.frame(Roots_fp = 1/polyroot(A1), Radius = abs(1/polyroot(A1)), Angle = Arg(1/polyroot(A1))/2/pi*Fs)
#Overfits as well bc the angles don't match the periodogram


#try 10
#Look at spectrum, just AR components for now
N <- length(ydiff)
s_period <- spectrum(ydiff,method="pgram", fast = FALSE, log = "no")
s_ar <- spectrum(ydiff,method="ar", type = "b", log = "yes", n.freq = (N/2)+1, order = 10)
Spectral <- data.frame(freq = s_period$freq, Period <- s_period$spec, AR <- s_ar$spec[2:(N/2+1)])
plot(ggplot(Spectral, aes(x=freq)) +  geom_line(aes(y = Period, colour = "periodogram")) + geom_line(aes(y=AR, colour = "AR")))
#AIC selects order 17 AR as well, not sure if signal is attenuated

AR1 <- ar(ydiff,aic = FALSE, order = 10, method = c("ols"))
A1 <- c(1,-AR1$ar)
AR_Model1 <-data.frame(Roots_fp = 1/polyroot(A1), Radius = abs(1/polyroot(A1)), Angle = Arg(1/polyroot(A1))/2/pi*Fs)
#Radius of the peak in the middle is too small + some nonexistent peaks, should go down even further

#order 8, 
N <- length(ydiff)
s_period <- spectrum(ydiff,method="pgram", fast = FALSE, log = "no")
s_ar <- spectrum(ydiff,method="ar", type = "b", log = "yes", n.freq = (N/2)+1, order = 6)
Spectral <- data.frame(freq = s_period$freq, Period <- s_period$spec, AR <- s_ar$spec[2:(N/2+1)])
plot(ggplot(Spectral, aes(x=freq)) +  geom_line(aes(y = Period, colour = "periodogram")) + geom_line(aes(y=AR, colour = "AR")))
#AIC selects order 17 AR as well, not sure if signal is attenuated

AR1 <- ar(ydiff,aic = FALSE, order = 6, method = c("ols"))
A1 <- c(1,-AR1$ar)
AR_Model1 <-data.frame(Roots_fp = 1/polyroot(A1), Radius = abs(1/polyroot(A1)), Angle = Arg(1/polyroot(A1))/2/pi*Fs)
#Cannot get it to have the peak in the middle as the highest radius. At 6 it only models real peaks but the AR graph 
#Looks completely wrong

#Try with the original series bc that's what it says in the assignment


N <- length(y)
s_period <- spectrum(y,method="pgram", fast = FALSE, log = "no")
s_ar <- spectrum(y,method="ar", type = "b", log = "yes", n.freq = (N/2)+1)
Spectral <- data.frame(freq = s_period$freq, Period <- s_period$spec, AR <- s_ar$spec[2:(N/2+1)])
plot(ggplot(Spectral, aes(x=freq)) +  geom_line(aes(y = Period, colour = "periodogram")) + geom_line(aes(y=AR, colour = "AR")))
#Low frequency stuff is apparently due to baroreflex and Mayer waves, AR order 11


AR1 <- ar(y,aic = FALSE, order = 11, method = c("ols"))
A1 <- c(1,-AR1$ar)
AR_Model1 <-data.frame(Roots_fp = 1/polyroot(A1), Radius = abs(1/polyroot(A1)), Angle = Arg(1/polyroot(A1))/2/pi*Fs)
#Radii are still not in the right order, maybe a bit lower

N <- length(y)
s_period <- spectrum(y,method="pgram", fast = FALSE, log = "no")
s_ar <- spectrum(y,method="ar", type = "b", log = "yes", n.freq = (N/2)+1, order = 9)
Spectral <- data.frame(freq = s_period$freq, Period <- s_period$spec, AR <- s_ar$spec[2:(N/2+1)])
plot(ggplot(Spectral, aes(x=freq)) +  geom_line(aes(y = Period, colour = "periodogram")) + geom_line(aes(y=AR, colour = "AR")))


AR1 <- ar(y,aic = FALSE, order = 9, method = c("ols"))
A1 <- c(1,-AR1$ar)
AR_Model1 <-data.frame(Roots_fp = 1/polyroot(A1), Radius = abs(1/polyroot(A1)), Angle = Arg(1/polyroot(A1))/2/pi*Fs)
#10 or 9 look better than 11, but not lower

#So educated guess based on order 10: 8AR terms, 6MA terms 
#Cholesky decomposition could also be useful here


#3
library(forecast)
automodel <- auto.arima(ysdiff)
#wtf it chooses just 1 ma and 1 seasonal ma component for ysdiff, but spectrum has 2 clear peaks which have to be AR

#Try with just the differenced series
automodel <- auto.arima(ydiff)
#Ok autoplot suggests 2 maxima and 2 minima, seems to correspond ok to periodogram

#Think I should have used the original series from the start, it decides differentiation for itself
automodel <- auto.arima((y))


#4
#This might be about generators, the plot with the roots on the unit circle

PredictTS(y,automodel,100,N,5)
