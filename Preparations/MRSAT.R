library(mrsat)

#loading the demo data
data(Auditory_demo)
View(Auditory_demo)

## define how conditions are grouped
my.signal <- list(noint = c(1,3), 
                  obrel = c(5,8), 
                  obrelsub = c(11, 14))
my.noise <- list(noint = c(2,4), 
                 obrel = c(6, 7, 9, 10), 
                 obrelsub = c(12, 13, 15, 16))

## get bins
d.bins <- get.bins(Auditory_demo, auditory=TRUE)
View(d.bins$bins)
#check how the RT are binned.
plot(d.bins$opt.bins)

## obtain dprime values
d.dprime <- get.dprime(d.bins$bins, signal.list=my.signal, noise.list=my.noise, 
                       is.index=TRUE, binmax=14)

View(d.dprime)
attributes(d.dprime)
?get.dprime

# x <- data.frame(0)
# attr(x, "class") <- c("mrsat.data", "data.frame")
# class(x)
# row

## define structure of parameters
## in this case, different parameters for each condition
pc333 <- list(asym=c(1, 2, 3), rate=c(1, 2, 3), incp=c(1, 2, 3))

## fit the curve assuming, and plot
fit.333 <- fit.SATcurve(d.dprime, par.cond = pc333)
plot(fit.333, main="333")

## compare that to the curves with different asymptotes but the same rate and intercept

pc311 <- list(asym=c(1, 2, 3), rate=c(1, 1, 1), incp=c(1, 1, 1))
fit.311 <- fit.SATcurve(d.dprime, par.cond = pc311)
plot(fit.311, main="311")

#compare outputs of the two models side by side
SATsummary.list(list(fit.333, fit.311))

#or just compare AIC of the two models
AIC(fit.333, fit.311)

#fitting a 311 model with fixed asymptote
fit.311fa <- fit.SATcurve(d.dprime, fix.asym=TRUE, par.cond = pc311, quiet = FALSE, trace = TRUE)
#and compare with the other two models
SATsummary.list(list(fit.333, fit.311, fit.311fa))


?`mrsat-package`
?Auditory_demo
?SATsummary.list
?acp
?average.dprime

#need to define the scale.list
exp1.signal <- list(noint = c(1,3), 
                    obrel = c(5,8), 
                    obrelsub = c(11, 14))

exp1.noise <- list(noint = c(2,4), 
                   obrel = c(6, 7, 9, 10), 
                   obrelsub = c(12, 13, 15, 16))

#load data

data(Auditory_demo)

### this isn't really meaningful, but for the purpose of demo,
### tag bins in two diffent ways (with "fixed" vs. "free" window)
### and then obtain two data frames containing slightly different
### dprime and lag values.

s01.bins.fixed <- get.bins(Auditory_demo, auditory=TRUE, window = "fixed")
s01.bins.free <- get.bins(Auditory_demo, auditory=TRUE,  window = "free")


s01.dp.fixed <- get.dprime(s01.bins.fixed$bins, 
                           signal.list = exp1.signal, noise.list = exp1.noise, is.index=TRUE,
                           binmax=14)

s01.dp.free <- get.dprime(s01.bins.free$bins, 
                          signal.list = exp1.signal, noise.list = exp1.noise, is.index=TRUE,
                          binmax=14)

### average the two data frames

mean.dp <- average.dprime(list(s01.dp.fixed, s01.dp.free))

plot(mean.dp)

?extractAIC.SATcurve
?fit.SATcurve

#load data
data(Auditory_demo)

#need to define the scale.list
exp1.signal <- list(noint = c(1,3), 
                    obrel = c(5,8), 
                    obrelsub = c(11, 14))

exp1.noise <- list(noint = c(2,4), 
                   obrel = c(6, 7, 9, 10), 
                   obrelsub = c(12, 13, 15, 16))

# tag bins
s01.bins <- get.bins(Auditory_demo, auditory=TRUE, window = "fixed")

#calculate dprimes
s01.dp<- get.dprime(s01.bins$bins, 
                    signal.list = exp1.signal, noise.list = exp1.noise, is.index = TRUE,
                    binmax=14)

# specify par.cond
my.pc <- list(asym=c(1,2,3), rate=c(1,2,2), incp=c(1,1,1))

# fit the data
s01.fit <- fit.SATcurve(s01.dp, par.cond = my.pc)

#plot the data and fitter curves
plot(s01.fit)

#to calculate AIC and BIC, do
AIC(s01.fit)
BIC(s01.fit)

#you may also calculate AIC analogue using MSE
# as suggested by Burnham and Anderson (1998)
AIC(s01.fit$MSE)

?get.bins

data(Auditory_demo)
s01.bins <- get.bins(Auditory_demo, auditory = TRUE)
View(s01.bins)
#this would provide a slightly different bins 
s01.bins.free <- get.bins(Auditory_demo, auditory = TRUE, window="free")
View(s01.bins.free)

?get.dprime

#need to define the scale.list
exp1.signal <- list(noint = c(1,3), 
                    obrel = c(5,8), 
                    obrelsub = c(11, 14))

exp1.noise <- list(noint = c(2,4), 
                   obrel = c(6, 7, 9, 10), 
                   obrelsub = c(12, 13, 15, 16))

#then load data and tag bins

data(Auditory_demo)
s01.bins <- get.bins(Auditory_demo, auditory=TRUE)

# finally, obtain dprime
s01.dp <- get.dprime(s01.bins$bins, 
                     signal.list = exp1.signal, noise.list = exp1.noise, 
                     is.index = TRUE, binmax=14)


plot(s01.dp)

?get.param

pc1 <- list(asym=c(1,2,3), rate=c(1,2,2), incp=c(1,1,1))
get.param(pc1)

pc2 <- list(asym=c(1,2,3,4,5), rate=c(1,2,2,3,3), incp=c(1,1,1,2,2))
get.param(pc2)

?logLik.SATcurve
?opt.bins

data(Auditory_demo)
rt_col <- grep("rt[0-3][0-9]",colnames(Auditory_demo), value=TRUE)
RT <- Auditory_demo[,rt_col]
d.opt <- opt.bins(unlist(na.omit(RT)))
plot(d.opt)

?plot.SATcurve
?plot.mrsat.data
?read.dat
?read.edat.txt
?read.param
?summary.SATcurve
?vareq


