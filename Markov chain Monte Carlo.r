set.seed(666)
setwd("C:/Users/Wei/Documents/Purdue STAT 695 Bayesian Data Analysis/HW3")

data_covariates = read.csv("ny_doe_covariates.csv", header=TRUE)
data_assignments = read.csv("ny_doe_treatment_assignments.csv", header=TRUE)
data_grades = read.csv("ny_doe_cumulative_grades.csv", header=TRUE)


library(standardize)
colNames = colnames(data_assignments)
data_assignments = scale(data_assignments)

X = cbind(rep(1, length(data_grades)), data_assignments)
colnames(X)[1] = "Intercept"


past_grades = data_covariates[11: 13]
past_grades = rowSums(past_grades)
Y = as.vector(data_grades > past_grades)
Y[Y==TRUE] = 1
Y[Y==FALSE] = 0

log_posterior = function(corr) {
  return(sum(dcauchy(corr, scale=2.5, log=TRUE)) +
         sum(Y * pnorm(X %*% corr, log.p=TRUE) + 
         (1 - Y) * pnorm(X %*% corr, lower.tail=FALSE, log.p=TRUE)))
}

optim_pars = optim(rep(0, ncol(X)), log_posterior, method="BFGS", control=list(fnscale=-1), hessian=TRUE)

initialValue = optim_pars$par
Hessian = optim_pars$hessian / (2.4 / sqrt(ncol(X))) # BDA pg. 290
burnIn = 10000
iterations = 20000

MCMC_MH = function(X, initialValue, Hessian) {
  chains = array(dim=c(iterations + 1, ncol(X)))
  chains[1, ] = initialValue
  Rchol = chol(-Hessian)
  BarProgress = txtProgressBar(min = 1, max = nrow(chains), style = 3)
  for (i in 1: iterations) {
    # never directly save the inverse of a matrix, compute them instead
    proposal = chains[i, ] + backsolve(Rchol, rnorm(ncol(X)))

    # better write exp(num) as num to avoid overflow; symmetric proposal
    log_acceptance_prob = log_posterior(proposal) - log_posterior(chains[i, ]) 

    chains[i + 1, ] = chains[i, ]
    if (log(runif(1)) < log_acceptance_prob)
      chains[i + 1, ] = proposal
    setTxtProgressBar(BarProgress, i)
  }
  close(BarProgress)
  return(chains[-(1: burnIn), ])  
}

pars_draws = MCMC_MH(X, initialValue, Hessian)
par(mfrow=c(1, 3))
for (i in 1: ncol(X)) {
  interval = quantile(pars_draws[ , i], probs=c(0.025, 0.975))
  if (interval[1] * interval[2] < 0) { # coefficient not significantly differs from 0
    hist(pars_draws[ , i], nclass=30, , main=colnames(X)[i], xlab="" )
    abline(v = 0, col="red")
  }
  else {
    print(colnames(X)[i])
  }
}


# Evaluate the main effects of QR, PA, SPBP with all two-factor interactions containing one of these factors. 
# Moreover, we exclude the intercept item since it is not significant.

QR = X[, 2]; PA = X[, 3]; IT = X[, 4]; SPBP = X[, 5]; ARIS = X[, 6]
XX = cbind(QR, PA, SPBP, QR * PA, QR * SPBP, PA * SPBP, 
          QR * IT, QR * ARIS, PA * IT, PA * ARIS, SPBP * IT, SPBP * ARIS)
colnames(XX) = c("QR", "PA", "SPBP", "QR*PA", "QR*SPBP", "PA*SPBP", 
                 "QR*IT", "QR*ARIS", "PA*IT", "PA*ARIS", "SPBP*IT", "SPBP*ARIS")

log_posterior = function(corr) {
  return(sum(dcauchy(corr, scale=2.5, log=TRUE)) +
         sum(Y * pnorm(XX %*% corr, log.p=TRUE) + 
         (1 - Y) * pnorm(XX %*% corr, lower.tail=FALSE, log.p=TRUE)))
}

optim_pars = optim(rep(0, ncol(XX)), log_posterior, method="BFGS", control=list(fnscale=-1), hessian=TRUE)

initialValue = optim_pars$par
Hessian = optim_pars$hessian / (2.4 / sqrt(ncol(XX))) # BDA pg. 290


pars_draws = MCMC_MH(XX, initialValue, Hessian)
par(mfrow=c(3, 4))
for (i in 1: ncol(XX)) {
  interval = quantile(pars_draws[ , i], probs=c(0.025, 0.975))
  if (interval[1] * interval[2] < 0) { # coefficient not significantly differs from 0
    hist(pars_draws[ , i], nclass=30, , main=colnames(XX)[i], xlab="" )
    abline(v = 0, col="red")
  }
  else {
    print(colnames(XX)[i])
  }
}



# Evaluate the main effects of QR, PA, SPBP, PA*IT and all three-factor interactions involving one of these factors.
# "QR"   "PA"   "IT"   "SPBP" "ARIS"

XX = cbind(QR, PA, SPBP, PA * IT, 
           QR*PA*IT, QR*PA*SPBP, QR*PA*ARIS, QR*IT*SPBP, QR*IT*ARIS, QR*SPBP*ARIS,
           PA*IT*SPBP, PA*IT*ARIS, PA*SPBP*ARIS,
           IT*SPBP*ARIS)
colnames(XX) = c("QR", "PA", "SPBP", "PA*IT", "QR*PA*IT", "QR*PA*SPBP", "QR*PA*ARIS", "QR*IT*SPBP", "QR*IT*ARIS", "QR*SPBP*ARIS",
    "PA*IT*SPBP", "PA*IT*ARIS", "PA*SPBP*ARIS", "IT*SPBP*ARIS")

log_posterior = function(corr) {
  return(sum(dcauchy(corr, scale=2.5, log=TRUE)) +
         sum(Y * pnorm(XX %*% corr, log.p=TRUE) + 
         (1 - Y) * pnorm(XX %*% corr, lower.tail=FALSE, log.p=TRUE)))
}

optim_pars = optim(rep(0, ncol(XX)), log_posterior, method="BFGS", control=list(fnscale=-1), hessian=TRUE)

initialValue = optim_pars$par
Hessian = optim_pars$hessian / (2.4 / sqrt(ncol(XX))) # BDA pg. 290


pars_draws = MCMC_MH(XX, initialValue, Hessian)
colnames(pars_draws) = colnames(XX)
par(mfrow=c(3, 4))
for (i in 1: ncol(XX)) {
  interval = quantile(pars_draws[ , i], probs=c(0.025, 0.975))
  if (interval[1] * interval[2] < 0) { # coefficient not significantly differs from 0
    hist(pars_draws[ , i], nclass=30, , main=colnames(XX)[i], xlab="" )
    abline(v = 0, col="red")
  }
  else {
    print(colnames(XX)[i])
  }
}




# Evaluate the effects of QR, PA, SPBP, PA*IT
XX = cbind(QR, PA, SPBP, PA * IT)
colnames(XX) = c("QR", "PA", "SPBP", "PA*IT")

#XX = data.frame(QR=QR, PA=PA, SPBP=SPBP, PA.IT=PA*IT)

log_posterior = function(corr) {
  return(sum(dcauchy(corr, scale=2.5, log=TRUE)) +
         sum(Y * pnorm(XX %*% corr, log.p=TRUE) + 
         (1 - Y) * pnorm(XX %*% corr, lower.tail=FALSE, log.p=TRUE)))
}

optim_pars = optim(rep(0, ncol(XX)), log_posterior, method="BFGS", control=list(fnscale=-1), hessian=TRUE)

initialValue = optim_pars$par
Hessian = optim_pars$hessian / (2.4 / sqrt(ncol(XX))) # BDA pg. 290


pars_draws = MCMC_MH(XX, initialValue, Hessian)
colnames(pars_draws) = colnames(XX)
pars_draws = data.frame(pars_draws)

quantile(pars_draws$QR, probs=c(0.025, 0.975))
quantile(pars_draws$PA, probs=c(0.025, 0.975))
quantile(pars_draws$SPBP, probs=c(0.025, 0.975))
quantile(pars_draws$PA.IT, probs=c(0.025, 0.975))



## posterior checks
XX = data.frame(XX)
group = list(QR=XX$QR>0, PA=XX$PA>0, SPBP=XX$SPBP>0, PA.IT=XX$PA.IT>0)
mat = matrix(NA, nrow=nrow(pars_draws), ncol=2)
posterior_check = list(QR=mat, PA=mat, SPBP=mat, PA.IT=mat)

BarProgress = txtProgressBar(min=1, max=nrow(pars_draws), style=3)
for (i in 1: nrow(pars_draws)) {
    diff = rep(NA, nrow(XX)) 
    for (j in 1: nrow(XX)) {
        prediction = pnorm(sum(XX[j, ] * pars_draws[i, ]))
        diff[j] = Y[j] - prediction
    }
    # posterior_check$element = [treatment group, control group]
    posterior_check$QR[i, ] = c(mean(diff[group$QR]), mean(diff[!group$QR])) 
    posterior_check$PA[i, ] = c(mean(diff[group$PA]), mean(diff[!group$PA])) 
    posterior_check$SPBP[i, ] = c(mean(diff[group$SPBP]), mean(diff[!group$SPBP])) 
    posterior_check$PA.IT[i, ] = c(mean(diff[group$PA.IT]), mean(diff[!group$PA.IT]))
    #setTxtProgressBar(BarProgress, i)
}
close(BarProgress)


mt = cbind(posterior_check$QR, posterior_check$PA, posterior_check$SPBP, posterior_check$PA.IT)

up_quantile = function(vec) quantile(vec, prob=0.975) 
down_quantile = function(vec) quantile(vec, prob=0.025)

evals = data.frame(mean=apply(mt, 2, mean),
                   up=apply(mt, 2, up_quantile),
                   down=apply(mt, 2, down_quantile))

row.names(evals) = c('QR_trt', 'QR_ctr', 'PA_trt', 'PA_ctr', 'SPBP_trt', 'SPBP_ctr', 'PA.IT_trt', 'PA.IT_ctr')

plot(evals$mean, ylim=c(min(evals$down), max(evals$up)), xaxt="n")
polygon(c(1: 8, rev(1: 8)),c(evals$down, rev(evals$up)),col=gray(0.8))
lines(evals$mean, pch=15, panel.first = grid())
points(evals$up)
points(evals$down)
abline(h=0, col="blue")
axis(1, at=1:8, labels=rownames(evals))



##　analyze the magnitude of the effect of interventions
XX = cbind(rep(1, nrow(XX)), QR, PA, SPBP, PA * IT)
colnames(XX) = c("Intercept", "QR", "PA", "SPBP", "PA*IT")
Y = (data_grades - past_grades)[1: 232, ]

# we have done the variable selection, we proceed by using Bayesian linear regression to estimate the coefficient.

# correlation is small, we are safe to use ridge regression
cor(XX[, 2:5]) 


coef = function(k) {
    R = t(XX) %*% XX
    inv = solve(R + k * diag(5)) # when k is large, inverse suffers less from numerical problem
    
    df = data.frame(mean=inv %*% t(XX) %*% Y,
                    sd=sqrt(diag(inv %*% R %*% inv)))
    round(df, 4)
}

coef(0.01)
coef(1)
coef(100)


