x <- -3
y <- 0
z <- 3

my_f <- function (x) {if(x < 0) {
  print("x is a negative number")
} else if(x == 0) {
  print("x is zero")
} else {
  print("x is a positive number")
}}

my_f(x)
my_f(y)
my_f(z)

install.packages("Z:/Downloads/satf_1.0.tar.gz", repos=NULL, dependencies=TRUE)

# SATF <- function(t, asymptote, invrate, intercept)
# (t >= intercept)*(asymptote*(1-exp(-1/invrate*(t-intercept))))

#asymptote*   ifelse(t > 0, (1-exp(-(slope)^(t-intercept))), 0)


satf <- function(t, intercept, slope, asymptote) {
  ifelse(t > intercept, asymptote*(1-exp(-slope*(t-intercept))), 0)
}

plot(function(t) {satf(t, 0.5, 1, 5)}, xlim=c(-2, 5))

?pnorm

p_yes_noise <- function(dprime, crit) {
  1 - pnorm(crit, mean = 0, sd = 1)
}

p_yes_signal <- function(dprime, crit) {
  1 - pnorm(crit, mean = dprime, sd = 1)
}

p_yes_noise(1, 0)
p_yes_signal(1, 0)

p_yes_dprime <- function(dprime, crit, is_signal) {
  1 - pnorm(crit, mean = is_signal*dprime, sd = 1)
}

dprime <- satf(seq(0, 4.5, 0.5), 0.5, 1, 5)
dprime

p_yes_dprime(dprime, 0.5*dprime, 0)

r_yes <- function(n, p_yes) {
  runif(n) <= p_yes
} 

x <- r_yes(20000, 0.3)
mean(x)

# random number generation uniform distribution
?runif