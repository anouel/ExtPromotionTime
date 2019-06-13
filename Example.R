
#----------------------------------------------------------------
# Simulating data.
#----------------------------------------------------------------

set.seed(123)
n <- 200
# covariate 1 and 2
X1 <- (runif(n) - .5) / sqrt(1 / 12)
X2 <- as.numeric(runif(n) < 0.5)
# censoring times: truncated exponential distribution
C <- round(rexp(n, 1 / 5), 5)
Cbin <- (C > 30)
while (sum(Cbin) > 0) {
  C[Cbin] <- round(rexp(sum(Cbin), 1 / 5), 5)
  Cbin <- (C > 30)
}
expb <- exp(0.5 + X1 - 0.5 * X2)
cure <- exp(-expb) # cure probabilities
# event times with baseline cdf of a truncated exponential
U <- runif(n)
T <- round(-6 * log(1 + (1 - exp(-20 / 6)) * log(1 - (1 - cure) * U) / expb), 5)
T[(runif(n) < cure)] <- 99999 # cured subjects
OT <- pmin(C, T) # observed times
delta <- (OT == T) # censoring indicator
Dat <- data.frame(OT, delta, X1, X2)
#----------------------------------------------------------------


#----------------------------------------------------------------
# Desired link function and its first and second derivatives.
#----------------------------------------------------------------

eta <- exp
deta <- exp
ddeta <- exp
# eta=function(...){exp(sin(...))};deta=function(...){eta(...)*cos(...)};ddeta=function(...){deta(...)*cos(...)-eta(...)*sin(...)}
g <- function(x, betA, ETA = eta) {
  as.vector(ETA(t(betA) %*% x))
}
dg <- function(x, betA) x * g(x, betA, ETA = deta)
ddg <- function(x, betA) (x %*% t(x)) * g(x, betA, ETA = ddeta)
#----------------------------------------------------------------

#----------------------------------------------------------------
# Getting the estimated parameters and there variances.
#----------------------------------------------------------------

source("Prog.R")
estF(Dat)
