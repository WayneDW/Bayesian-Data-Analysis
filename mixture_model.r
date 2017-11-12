set.seed(666)
setwd("C:/Users/Wei/Documents/Purdue STAT 695 Bayesian Data Analysis/HW5")

data = read.csv(file="mix_reg.txt", header=TRUE)

X = data$x
y = data$y[order(X)]
X = sort(X)

n = nrow(data)

X = cbind(1, X, X^2)
### Part 1

shape = 1
scale = 1000
std = 1000

library(invgamma)

sigmoid = function(x) 1 / (1 + exp(-x))

log_priors = function(pars) {
    sum(dnorm(pars[1:5], mean=0, sd=std, log=T)) + 
    sum(dinvgamma(exp(pars[6:7]), shape, rate=1, scale=scale, log=TRUE)) +
    dbeta(sigmoid(pars[8]), 1 , 1, log=TRUE)
}

# joint likelihood
log_likelihood = function(pars) {
    mu1 = c(pars[1:2], 0)
    mu2 = pars[3:5]
    sd1 = exp(pars[6])
    sd2 = exp(pars[7])
    lambda = sigmoid(pars[8])
    sum((dnorm(y, mean=X %*% mu1, sd=sd1, log=T) + log(lambda)) * Z) + 
    sum((dnorm(y, mean=X %*% mu2, sd=sd2, log=T) + log(1 - lambda)) * (1 - Z))
}

log_posterior = function(pars) log_priors(pars) + log_likelihood(pars)

pars = rep(0.5, 8)
Z = rbinom(n, 1, 0.5)

burnIn = 1000
iterations = 2 * burnIn

log_posterior(pars)
for (i in 1: 10) {
    optimal = optim(pars, log_posterior, control=list(fnscale=-1), hessian=TRUE)
    pars = optimal$par
    log_posterior_raw = log_posterior(pars)
    chains = array(dim=c(iterations + 1, 8))
    chains[1, ] = pars
    for (j in 1: iterations) {
        # better avoid saving the inverse of a matrix, compute them instead
        proposal = chains[j, ] + rnorm(8, sd=0.1)

        # write exp(num) as num to avoid overflow; symmetric proposal
        log_acceptance_prob = log_posterior(proposal)  -  log_posterior(chains[j, ])
        chains[j + 1, ] = chains[j, ]
        if (log(runif(1)) < log_acceptance_prob)
            chains[j + 1, ] = proposal
    }
    pars_draws = chains[-(1: burnIn), ]
    print(paste(i, "th round: ", "Acceptance rate", round(nrow(unique(chains)) / nrow(chains), 4)))
    pars = tail(pars_draws, 1)


    Z = dnorm(y, X %*% c(pars[1:2], 0), sd=exp(pars[6])) > dnorm(y, X %*% pars[3:5], sd=exp(pars[7])) # if the probability one belongs to one group is larger than another
    Z = as.numeric(Z)
    log_posterior_update = log_posterior(pars)
    print(c(log_posterior_raw, log_posterior_update))
}

plot(X[, 2] * Z, y * Z, ylim=c(-10, 60), col="red", pch=19)
points(X[, 2] * (1 - Z), y * (1 - Z), col="black", pch=15)


qt = array(NA, c(200, 2, 3))
pars_draws[, 6] = exp(pars_draws[, 6])
pars_draws[, 7] = exp(pars_draws[, 7])

for (i in 1:200) {
    beta = cbind(pars_draws[, 1:2], 0)
    std = sqrt(pars_draws[, 6])
    y_samples = beta %*% X[i, ] + rnorm(burnIn+1, sd=std)
    qt[i, 1, ] = quantile(y_samples, c(0.05, 0.5, 0.95))

    beta = pars_draws[, 3:5]
    std = sqrt(pars_draws[, 6])
    y_samples = beta %*% X[i, ] + rnorm(burnIn+1, sd=std)
    qt[i, 2, ] = quantile(y_samples, c(0.05, 0.5, 0.95))
}

plot(X[, 2], y, xlab='X', ylab='y')
lines(X[, 2], qt[, 1, 1], col='red')
lines(X[, 2], qt[, 1, 2], col='red')
lines(X[, 2], qt[, 1, 3], col='red')
lines(X[, 2], qt[, 2, 1], col='blue')
lines(X[, 2], qt[, 2, 2], col='blue')
lines(X[, 2], qt[, 2, 3], col='blue')
