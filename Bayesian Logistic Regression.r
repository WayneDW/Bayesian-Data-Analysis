
######################################################################################################################################
### Computational problem
######################################################################################################################################

# set.seed(8) # can warm up but the result is bad
set.seed(26)
#setwd("//myhome.itap.purdue.edu/puhome/pu.data/Desktop/STAT 695/HW1")
setwd("C:/Users/Wei/Desktop/STAT 695/HW1")

###################################################################
# Q1, build the prior
###################################################################

# Suppose we have normal priors for beta_0 and beta_1
set_sol = seq(0, 3, 0.001)
thres_beta_0 = sapply(set_sol, function(sigma) abs(quantile(rnorm(1000, log(0.3), sigma), 0.9) - 0.4))
sigma_0 = set_sol[which.min(thres_beta_0)]

thres_beta_1 = sapply(set_sol, function(sigma) abs(quantile(rnorm(1000, log(1.5), sigma), 0.95) - 2))
sigma_1 = set_sol[which.min(thres_beta_1)]


prior_beta_0 = function(beta_0) dnorm(beta_0, log(0.3), sigma_0)
prior_beta_1 = function(beta_1) dnorm(beta_1, log(1.5), sigma_1)

x_grid = seq(-2.5, -0.5, 0.03)
y_grid = seq(-1.5, 1.5, 0.03)

points_beta_0 = sapply(x_grid, prior_beta_0)
points_beta_1 = sapply(y_grid, prior_beta_1)

par(mfrow=c(1,2))
plot(x_grid, points_beta_0, type="l", col="blue", main="density for beta_0")
plot(y_grid, points_beta_1, type="l", col="blue", main="density for beta_1")

###################################################################
# Q2, build the likelihood
###################################################################
# reference: https://www.mathworks.com/help/stats/examples/bayesian-analysis-for-a-logistic-regression-model.html?requestedDomain=www.mathworks.com
data = read.csv(file="computation_data.csv", header=TRUE)
x = data$x
y = data$y

p_ratio = function(beta_0, beta_1) 1 / (1 + exp(- beta_0 - beta_1 * x))

# The following is a bad way, due to the numerical problem when there are too many points
# likelihood = function(beta_0, beta_1) prod(dbinom(y, 1, p_ratio(beta_0, beta_1)))
log_likelihood = function(beta_0, beta_1) sum(log(dbinom(y, 1, p_ratio(beta_0, beta_1)))) # solve log likelihood first


#install.packages("scatterplot3d")
library(scatterplot3d)

# beta_0 ~ norm(0.3, 0.077), beta_1 ~ norm(1.5, 0.3)
x_grid = seq(-2.5, -0.5, 0.03)
y_grid = seq(-1.5, 0, 0.03)
grids = expand.grid(x_grid,y_grid) # create a grid for the 2-d tuples
val_log_likelihood = mapply(log_likelihood, grids$Var1, grids$Var2) # mapply is used to deal with multiple parameters
val_log_likelihood = val_log_likelihood - max(val_log_likelihood) # !!! very important step, scale the likelihood

par(mfrow=c(1,2))
scatterplot3d(grids$Var1, grids$Var2, val_log_likelihood, color="lightblue",pch=21,main="Log Likelihood")
scatterplot3d(grids$Var1, grids$Var2, exp(val_log_likelihood), color="lightblue",pch=21,main="Likelihood")

###################################################################
# Q3, build posterior, and maximize the likelihood
###################################################################

# example: https://www.r-bloggers.com/how-to-use-optim-in-r/


log_posterior = function(beta) {
  return (log(prior_beta_0(beta[1])) 
        + log(prior_beta_1(beta[2])) 
        + log_likelihood(beta[1], beta[2]))
}

optimal = optim(c(0, 0), log_posterior, control=list(fnscale=-1), hessian=TRUE)

u = optimal$par
hessian = optimal$hessian
cov = -solve(hessian)
u
hessian
cov


###################################################################
# Q4, draw posterior using a discrete grid approximation
###################################################################

num_draws = 1000
sample_posterior = function(beta) exp(log_posterior(beta))

# use MCMC to draw samples from log density with multiple parameters
# library(DPpackage)
# mcmc = list(nburn=2000, nsave=num_draws, ndisplay=500)
# grids = expand.grid(seq(0, 0.5, 0.01), seq(0, 0.5, 0.01))
# support = cbind(grids$Var1, grids$Var2)
# fit = PTsampler(sample_posterior, dim.theta=2, mcmc=mcmc, support=support)
# samples_mcmc_posterior = fit$save.state$thetasave

inner_trials = 100000
seq0 = sample(seq(u[1] - sqrt(cov[1,1]) * 3, u[1] + sqrt(cov[1,1]) * 3, 0.0001), replace=TRUE, size = inner_trials)
seq1 = sample(seq(u[2] - sqrt(cov[2,2]) * 3, u[2] + sqrt(cov[2,2]) * 3, 0.0001), replace=TRUE, size = inner_trials)
grid = rbind(seq0,seq1)
den_prop = apply(grid, 2, sample_posterior)
idx = sample(seq(1,inner_trials), size =1000, replace = TRUE, prob = den_prop)

samples_posterior = t(grid[,idx])

library(RColorBrewer)
rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
r = rf(32)
library(MASS)
k = kde2d(samples_posterior[,1], samples_posterior[,2], n=200)
image(k, col=r)

################################################################### 
# Q5, draw 2-d normal approximation
###################################################################

# multivariate normal distribution
# install.packages("mvtnorm")
C = chol(cov) 
mu = matrix(rep(u, num_draws), nrow=2)
yi = rbind(rnorm(num_draws), rnorm(num_draws))

samples_2d_normal = t(mu + C %*%  yi)
#library(mvtnorm)
#samples_2d_normal = rmvnorm(num_draws, u, covariance)

par(mfrow=c(1,1))
k = kde2d(samples_2d_normal[,1], samples_2d_normal[,2], n=200)
image(k, col=r)
###################################################################
# Q6, compare the summaries of the two draws
###################################################################

# to compute skewness, etc..
#install.packages("moments")
library(moments)

mySummary = function(dat) {
  values = round(c(quantile(dat, probs=c(0.025,0.25,0.5,0.75,0.975)), 
             mean(dat), sd(dat), skewness(dat), kurtosis(dat)), 4)
  return(setNames(values, c("2.5%", "25%", "50%", "75%", "97.5%", 
                            "mean", "sd", "skew", "kurt")))
}

mySummary(samples_posterior[,1])
mySummary(samples_2d_normal[,1])

###################################################################
# Q7, summarize the inference on exp(beta_1) 
#     describe merits for normal approximation
###################################################################

mySummary = function(dat) {
  values = round(c(quantile(dat, probs=c(0.025,0.25,0.5,0.75,0.975)), 
             mean(dat), sd(dat), skewness(dat), kurtosis(dat)), 4)
  return(setNames(values, c("2.5%", "25%", "50%", "75%", "97.5%", 
                            "mean", "sd", "skew", "kurt")))
}

# inference for exp(beta_1)
mySummary(exp(samples_posterior[,2]))
mySummary(exp(samples_2d_normal[,2]))

# drawing samples from multivariate normal is super fast

###################################################################
# Q8, evaluate the posterior
###################################################################

# use the samples drawn from posterior to simulate data from the function of binary response.

simulated_odd = matrix(NA, nrow = num_draws, ncol = length(y))
for (i in 1: num_draws)
  simulated_odd[i,] = p_ratio(samples_posterior[i,1], samples_posterior[i,2])
  
simulation = list(mu = c(), up = c(), low = c())

for (i in 1:length(y)) {
  simulation$mu = c(simulation$mu, rbinom(1, 1, mean(simulated_odd[,i])))
  simulation$up = c(simulation$up, rbinom(1, 1, quantile(simulated_odd[,i], probs=0.975)))
  simulation$low = c(simulation$low, rbinom(1, 1, quantile(simulated_odd[,i], probs=0.025)))
}


simulation$x = x
simulation$y = y
eval = list(p_real=c(), p_mu=c(), p_up=c(), p_low=c())

unique_v = c(simulation$x[1], 0, simulation$x[71])

for (v in unique_v) {
    eval$p_real = c(eval$p_real, sum(simulation$y[simulation$x==v]) / length(simulation$y[simulation$x==v]))
    eval$p_mu = c(eval$p_mu, sum(simulation$mu[simulation$x==v]) / length(simulation$mu[simulation$x==v]))
    eval$p_up = c(eval$p_up, sum(simulation$up[simulation$x==v]) / length(simulation$up[simulation$x==v]))
    eval$p_low = c(eval$p_low, sum(simulation$low[simulation$x==v]) / length(simulation$low[simulation$x==v]))
}

plot(eval$p_real, pch=21, ylim=c(0, 1), xaxt="n", main="Probabilities of success at different x")
axis(1, at=1:3, labels=unique_v)
points(eval$p_mu, type="l", col="red", lty=3)
points(eval$p_up, type="o", col="blue", lty=2)
points(eval$p_low, type="o", col="blue", lty=2)

# to plot histogram with binary response
#install.packages("popbio")
# library(popbio)

# par(mfrow=c(2,2))
# logi.hist.plot(x, y, boxp=FALSE, type="hist", col="gray", main="Data v.s. Simulations")
# logi.hist.plot(x, simulation$mu, boxp=FALSE, type="hist", col="gray", main="Simulation mean")
# logi.hist.plot(x, simulation$up, boxp=FALSE, type="hist", col="gray", main="Simulation upper")
# logi.hist.plot(x, simulation$low, boxp=FALSE, type="hist", col="gray", main="Simulation lower")

# par(mfrow=c(2,2))
# plot(x, y, main="Data v.s. Simulations")
# points(x, simulation$mu, col="blue", lty = 1)
# points(x, simulation$up, col="blue", lty = 1)
# points(x, simulation$low, col="blue", lty = 1)

###################################################################
# Q9, comment on the result and make suggestions
###################################################################

# The model performs well at point -1.216553, and fails at point 0 and 1.216553
# my suggestion of modifying this model is change a new type of prior

prior_beta_0 = function(beta_0) dnorm(beta_0, log(0.3), sigma_0)
prior_beta_1 = function(beta_1) dnorm(beta_1, log(2.0), sigma_1)

sample_posterior = function(beta) {
  return (exp(log(prior_beta_0(beta[1])) 
        + log(prior_beta_1(beta[2])) 
        + log_likelihood(beta[1], beta[2])))
}

inner_trials = 100000
seq0 = sample(seq(u[1] - sqrt(cov[1,1]) * 3, u[1] + sqrt(cov[1,1]) * 3, 0.0001), replace=TRUE, size = inner_trials)
seq1 = sample(seq(u[2] - sqrt(cov[2,2]) * 3, u[2] + sqrt(cov[2,2]) * 3, 0.0001), replace=TRUE, size = inner_trials)
grid = rbind(seq0,seq1)
den_prop = apply(grid, 2, sample_posterior)
idx = sample(seq(1,inner_trials), size =1000, replace = TRUE, prob = den_prop)

samples_posterior = t(grid[,idx])


library(RColorBrewer)
rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
r = rf(32)
library(MASS)
k = kde2d(samples_posterior[,1], samples_posterior[,2], n=200)
image(k, col=r)


simulation_samples = matrix(NA, nrow = num_draws, ncol = length(y))

for (i in 1: num_draws)
  simulation_samples[i,] = p_ratio(samples_posterior[i,1], samples_posterior[i,2])
  
simulation = list(mu = c(), up = c(), low = c())

for (i in 1:length(y)) {
  simulation$mu = c(simulation$mu, rbinom(1, 1, mean(simulation_samples[,i])))
  simulation$up = c(simulation$up, rbinom(1, 1, quantile(simulation_samples[,i], probs=0.975)))
  simulation$low = c(simulation$low, rbinom(1, 1, quantile(simulation_samples[,i], probs=0.025)))
}


simulation$x = x
simulation$y = y
eval = list(p_real=c(), p_mu=c(), p_up=c(), p_low=c())

unique_v = c(simulation$x[1], 0, simulation$x[71])

for (v in unique_v) {
    eval$p_real = c(eval$p_real, sum(simulation$y[simulation$x==v]) / length(simulation$y[simulation$x==v]))
    eval$p_mu = c(eval$p_mu, sum(simulation$mu[simulation$x==v]) / length(simulation$mu[simulation$x==v]))
    eval$p_up = c(eval$p_up, sum(simulation$up[simulation$x==v]) / length(simulation$up[simulation$x==v]))
    eval$p_low = c(eval$p_low, sum(simulation$low[simulation$x==v]) / length(simulation$low[simulation$x==v]))
}

plot(eval$p_real, pch=13, col="red", ylim=c(0, 1), xaxt="n", main="Probabilities at different x")
axis(1, at=1:3, labels=unique_v)
points(eval$p_mu, type="l", col="blue", lty=3)
points(eval$p_up, type="o", col="grey", lty=2)
points(eval$p_low, type="o", col="grey", lty=2)


######################################################################################################################################
### Computational problem
######################################################################################################################################

# "mother_op_sex" "night_wo_perm" "preg"  

###################################################################
# Q1
###################################################################

# Because Wave II happens one year later after Wave I, there would be more pregnant respondent. In addition, the predictor variable was collected 
# from Wave II, therefore, using response variable in wave II is more related to variable for whether the respondent spent a night away from home
# without permission in the last 12 months.

###################################################################
# Q2
###################################################################

# Sufficient, even though we don't have too much knowledge about priors on variables_mother_op_sex and night_wo_perm, we could still use non-informative prior, e.g. t/Cauchy/Lapace prior, to solve this problem.


###################################################################
# Q3
###################################################################

# Use the response and predictors are all binary, we will use binary logistic regression. Since we don't have too much information about the priors, weekly informative prior Cauchy prior is used.

# Pr (yi = 1 | beta; xi) = exp(0 + 1xi) / 1 + exp(0 + 1xi);

###################################################################
# Q4
###################################################################

# The estimand of interest is exp (beta_1), which captures the multiplicative change in the odds of success for a
# level change in the treatment.


###################################################################
# Q5
###################################################################

# This reference suggested to use Cauchy prior, since we don't have too much information about priors,
# https://tgmstat.wordpress.com/2014/03/19/weakly-informative-priors-for-logistic-regression/

# set.seed(8) # can warm up but the result is bad
set.seed(26)
#setwd("//myhome.itap.purdue.edu/puhome/pu.data/Desktop/STAT 695/HW1")
setwd("C:/Users/Wei/Desktop/STAT 695/HW1")
data = read.csv(file="pregnancy_data_subset.csv", header=TRUE)

names(data)
prior_intercept = function(beta_0) dcauchy(beta_0)
prior_mother_op_sex = function(beta_1) dcauchy(beta_1)
prior_night_wo_perm = function(beta_2) dcauchy(beta_2)

p_ratio = function(beta_0, beta_1, beta_2) 1 / (1 + exp(- beta_0 - beta_1 * data$mother_op_sex - beta_2 * data$night_wo_perm))
log_likelihood = function(beta_0, beta_1, beta_2) sum(log(dbinom(data$preg, 1, p_ratio(beta_0, beta_1, beta_2))))

log_posterior = function(beta) {
  return (log(prior_intercept(beta[1]))
        + log(prior_mother_op_sex(beta[2]))
        + log(prior_night_wo_perm(beta[3]))
        + log_likelihood(beta[1], beta[2], beta[3]))
}



optimal = optim(c(0, 1, 1), log_posterior, control=list(fnscale=-1), method="L-BFGS-B", hessian=TRUE)
optimal

u = optimal$par
hessian = optimal$hessian
cov = solve(-hessian)


###################################################################
# Q6
###################################################################
trials = 10000
seq1 = sample(seq(u[1] - sqrt(cov[1,1]) * 3, u[1] + sqrt(cov[1,1]) * 3, 0.0001), replace=TRUE, size = trials)
seq2 = sample(seq(u[2] - sqrt(cov[2,2]) * 3, u[2] + sqrt(cov[2,2]) * 3, 0.0001), replace=TRUE, size = trials)
seq3 = sample(seq(u[3] - sqrt(cov[3,3]) * 3, u[3] + sqrt(cov[3,3]) * 3, 0.0001), replace=TRUE, size = trials)
# grid = rbind(rep(u[1], trials), seq2, seq3)
grid = rbind(seq1, seq2, seq3)
dLogProp = apply(grid, 2, log_posterior)
cNormorlize = max(dLogProp)

sample_posterior = function(beta) exp(log_posterior(beta) - cNormorlize)
dProp = apply(grid, 2, sample_posterior)
idx = sample(seq(1, trials), size =1000, replace = TRUE, prob = dProp)
samples_posterior = t(grid[,idx])
vals = dProp[idx]

library(threejs)
ra = ceiling(256 * vals / max(vals))
col = rainbow(256, 2/3)

scatterplot3js(x=samples_posterior[,1], y=samples_posterior[,2], z=samples_posterior[,3], size=0.4, color = col[ra])


###################################################################
# Q7, evaluate the posterior
###################################################################

# use the samples drawn from posterior to simulate data from the function of binary response.
num_draws = 1000
simulated_p = matrix(NA, nrow = num_draws, ncol = dim(data)[1])

for (i in 1: num_draws)
  simulated_p[i,] = p_ratio(samples_posterior[i,1], samples_posterior[i,2], samples_posterior[i,3])
  
simulation = list(mu = c(), up = c(), low = c())

for (i in 1: dim(data)[1]) {
  simulation$mu = c(simulation$mu, mean(simulated_p[,i]))
  simulation$up = c(simulation$up, quantile(simulated_p[,i], probs=0.975))
  simulation$low = c(simulation$low, quantile(simulated_p[,i], probs=0.025))
}


xx = data$mother_op_sex
yy = data$night_wo_perm
simulation$real = data$preg

eval = list(p_real=c(), p_mu=c(), p_low=c(), p_up=c())
unique_v = c(0, 1)
for (v in unique_v) {
  for (w in unique_v) {
    eval$p_real = c(eval$p_real, mean(simulation$real[xx == v & yy == w]))
    eval$p_mu = c(eval$p_mu, mean(simulation$mu[xx == v & yy == w]))
    eval$p_low = c(eval$p_low, mean(simulation$low[xx == v & yy == w]))
    eval$p_up = c(eval$p_up, mean(simulation$up[xx == v & yy == w]))
  }
}



plot(eval$p_real, pch=13, col="red", ylim=c(0, 0.5), xaxt="n", ylab="", xlab="0: disapprove, 1: approve", main="Probabilities of pregnancy at different (mother_op_sex, night_wo_perm)")
axis(1, at=1:4, labels=c("(0,0)", "(0,1)", "(1,0)", "(1,1)"))
points(eval$p_mu, type="l", col="blue", lty=3)
points(eval$p_up, type="o", col="grey", lty=2)
points(eval$p_low, type="o", col="grey", lty=2)


###################################################################
# Q8, explain
###################################################################

# Parental supervision could reduce teenage pregnancy rate as much as 30% if her mother disapprove her having sex at this time in her life and disapprove teenage spent a night away from home.


