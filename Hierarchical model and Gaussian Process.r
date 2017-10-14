
######################################################################################################################################
### Computational problem
######################################################################################################################################

set.seed(26)
setwd("C:/Users/Wei/Desktop/STAT 695/HW2")

computation_data = read.csv(file="AirPassengers.csv", header=TRUE)

X = computation_data$X
Y = computation_data$AirPassengers
n = length(Y)

recursion = function(y, alpha) {
    if (y > 0) return((y + alpha - 1) / y * recursion(y - 1, alpha))
    else return(1)
}


# Q1 
log_likelihood = function(alpha, beta) {
  n * (alpha * log(beta)) - (sum(Y) + n * alpha) * log(beta + 1) + sum(log(mapply(recursion, Y, alpha)))
}


nums = 10
alpha_grid = seq(0.01, 1, length.out=nums)
beta_grid = seq(0.01, 1, length.out=nums)


grids = expand.grid(alpha_grid, beta_grid)
val_log_likelihood = matrix(mapply(log_likelihood, grids$Var1, grids$Var2), nrow=nums)


par(mfrow=c(1,2))
contour(alpha_grid, beta_grid, val_log_likelihood,
        xlab=bquote(paste(alpha)),
        ylab=bquote(paste(beta)),
        main=bquote(paste("Log-likelihood based on "*alpha*" and "*beta)))

# Q2
library(KScorrect)

LIMIT = 10
alpha_grid = seq(0.01, 10, 1)
log_beta_grid = seq(-LIMIT, LIMIT, 0.1)
grids = expand.grid(alpha_grid, beta_grid)
val_log_likelihood = matrix(mapply(log_likelihood, grids$Var1, grids$Var2), nrow=nums)

log_prior = matrix(NA, nrow=length(alpha_grid), ncol=length(log_beta_grid))

for (i in 1: length(alpha_grid)) {
  for (j in 1: length(log_beta_grid)) {
    log_prior[i,j] = log(1/1000) + log(dlunif(exp(log_beta_grid[j]), exp(-LIMIT), exp(LIMIT), base = exp(1)))
  }
}

par(mfrow=c(1,2))
contour(alpha_grid, exp(log_beta_grid), log_prior,
        xlab=bquote(paste(alpha)),
        ylab=bquote(paste(beta)),
        main=bquote(paste("Log-likelihood based on "*alpha*" and "*beta)))

# Q3
library(scatterplot3d)
log_posterior_alpha_log_beta = function(alpha, log_beta) {
  log_likelihood(alpha, exp(log_beta)) + log(dlunif(exp(log_beta), exp(-LIMIT), exp(LIMIT), base = exp(1)))
}

nums = 100
alpha_grid = seq(3.5, 8.5, length.out=nums)
log_beta_grid = seq(-4.3, -3.5, length.out=nums)
grids = expand.grid(alpha_grid, log_beta_grid) # create a grid for the 2-d tuples
val_log_posterior = mapply(log_posterior_alpha_log_beta, grids$Var1, grids$Var2) 
normalization_item = max(val_log_posterior)
val_log_posterior = val_log_posterior - normalization_item # normalize

par(mfrow=c(1,2))
scatterplot3d(grids$Var1, exp(grids$Var2), val_log_posterior, color="lightblue",pch=21, xlab=expression(alpha), ylab=expression(beta), zlab="", main="Log Posterior after normalization")
scatterplot3d(grids$Var1, exp(grids$Var2), exp(val_log_posterior), color="lightblue",pch=21, xlab=expression(alpha), ylab=expression(beta), zlab="", main="Posterior after normalization")


# Q4
log_posterior_alpha_beta = function(par) log_posterior_alpha_log_beta(par[1], log(par[2]))
optimal = optim(c(1, 1), log_posterior_alpha_beta, control=list(fnscale=-1), hessian=TRUE)
u = optimal$par
hessian = optimal$hessian
cov = -solve(hessian)


trials = 1000
sample_posterior = function(par) exp(log_posterior_alpha_beta(par) - normalization_item) # for fear that the highest likelihood is still too small
alpha_grid = sample(seq(u[1] - sqrt(cov[1,1]) * 3, u[1] + sqrt(cov[1,1]) * 3, length.out=trials), replace=TRUE, size = trials)
beta_grid = sample(seq(u[2] - sqrt(cov[2,2]) * 3, u[2] + sqrt(cov[2,2]) * 3, length.out=trials), replace=TRUE, size = trials)
grid = rbind(alpha_grid, beta_grid)
den_prop = apply(grid, 2, sample_posterior)
idx = sample(seq(1, trials), size =1000,replace = TRUE, prob = den_prop)
samples_posterior = t(grid[,idx])


C = chol(cov) 
num_draws = 1000
mu = matrix(rep(u, num_draws), nrow=2)
yi = rbind(rnorm(num_draws), rnorm(num_draws))
samples_2d_normal = t(mu + C %*%  yi)


library(RColorBrewer)
library(MASS)
rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
r = rf(32)
par(mfrow=c(1,2))
k = kde2d(samples_posterior[,1], samples_posterior[,2], n=200)
image(k, col=r, xlab=expression(alpha), ylab=expression(beta), main="Sampling from direct simulation")
k = kde2d(samples_2d_normal[,1], samples_2d_normal[,2], n=200)
image(k, col=r, xlab=expression(alpha), ylab=expression(beta), main="Sampling from multivariate normal")
# Q5



# Q6
library(moments)
mySummary = function(dat) {
  r = list()
  r$quantitle = rbind(quantile(dat[,1], probs=c(0.025,0.25,0.5,0.75,0.975)), quantile(dat[,2], probs=c(0.025,0.25,0.5,0.75,0.975)))
  r$mean = c(mean(dat[,1]), mean(dat[,2]))
  r$sd = c(sd(dat[,1]), sd(dat[,2]))
  r$skewness = skewness(dat)
  r$kurtosis = kurtosis(dat)
  r$correlation_coefficient = cor(dat[,1], dat[,2])
  rownames(r$quantitle) = names(r$mean) = names(r$sd) = names(r$skewness) = names(r$kurtosis) = c("alpha", "beta")
  return(mapply(round, r, 3))
}


mySummary(samples_posterior)
mySummary(samples_2d_normal)

# Q7
expectation = function(alpha, beta) (Y[1] + alpha) / (beta + 1)

gnum = 10
alpha_grid = seq(0.01, 10, length.out=gnum)
beta_grid = exp(seq(-10, 4, length.out=gnum))

conditional_mean_lambda_evaluations = matrix(NA, nrow=length(alpha_grid), ncol=length(beta_grid))
for (i in 1: gnum) {
  for (j in 1: gnum) {
    conditional_mean_lambda_evaluations[i,j] = expectation(alpha_grid[i], beta_grid[j])
  }
}

contour(alpha_grid, beta_grid, conditional_mean_lambda_evaluations,
        xlab=bquote(paste(alpha)),
        ylab=bquote(paste(beta)),
        main=bquote(paste("Posterior Expectation of "*lambda*" based on "*alpha*" and "*beta)))




######################################################################################################################################
### Applied problem
######################################################################################################################################
######################################################################################################################################
### Applied problem
######################################################################################################################################
set.seed(26)
require(ggplot2)
setwd("C:/Users/Wei/Desktop/STAT 695/HW2")
data = read.csv(file="rslt_V2.csv", header=FALSE, sep = " ")
x = data[,1:2]
Y = data[,3]

Y = Y[Y > 10]
x = x[Y > 10,]

sigma_n = 5
sigma_f = 20
l = 0.5

X = matrix(rep(0, 2 * nrow(x)), ncol=2)
for (i in 1: ncol(X)) X[, i] = (x[, i] - mean(x[, i])) / sd(x[, i])


# http://www.cs.toronto.edu/~duvenaud/cookbook/
# Calculate squared exponential kernel
# l: lengthscale, determines the length of the 'wiggles' in your function
# sigma_f: average distance of your function away from its mean
rbfKernel = function(X1, X2, l = 1, sigma_f = 1) {
  Kse = matrix(rep(0, nrow(X1) * nrow(X2)), nrow = nrow(X1))
  for (i in 1: nrow(Kse)) {
    for (j in 1: ncol(Kse)) {
      Kse[i,j] = exp(-0.5 * sum((X1[i, ] - X2[j, ])^2) / l^2) * sigma_f^2
    }
  }
  return(Kse)
}

X_pred = X

var_obs = rbfKernel(X, X, l, sigma_f)
Kconv = rbfKernel(X, X_pred, l, sigma_f)
var_pred = rbfKernel(X_pred, X_pred, l, sigma_f)
f_pred = as.vector(mean(Y) + t(Kconv) %*% solve(var_obs + sigma_n^2 * diag(1, ncol(var_obs))) %*% (Y - mean(Y)))
fvar_pred = var_pred - t(Kconv) %*% solve(var_obs + sigma_n^2 * diag(1, ncol(var_obs))) %*% Kconv

f_pred
fvar_pred

L = chol(var_obs + sigma_n^2 * diag(1, ncol(var_obs)))
invL = solve(L)
alpha = invL %*% t(invL) %*% (Y - mean(Y))
f_pred = as.vector(mean(Y) + t(Kconv) %*% alpha)
v_ = t(invL) %*% Kconv
fvar_pred = var_pred - t(v_) %*% v_

log_marginal_likelihood = function(X, Y) {
  -0.5 * t(Y - mean(Y)) %*% alpha - sum(log(diag(L))) - length(Y) / 2 * log(2 * pi)
}




nums = length(Y)
expected_Y = rep(NA, nums)
outlier_Y = rep(NA, nums)
error_bar = 2*sqrt(diag(fvar_pred))
iff = (Y >= f_pred - error_bar) & (Y <= f_pred + error_bar)
expected_Y[iff] = Y[iff]
outlier_Y[!iff] = Y[!iff]


ggplot() + 
geom_ribbon(data=NULL, aes(x=seq(1, nums), y=f_pred, ymin=f_pred - error_bar, ymax=f_pred + error_bar), fill="grey80") + # error bar
geom_line(data=NULL, aes(x=seq(1, nums), y=f_pred), size=1) + # mean
geom_point(data=NULL, aes(x=seq(1, nums), y=expected_Y), size=3, shape=23, fill="blue") +  # as expected
geom_point(data=NULL, aes(x=seq(1, nums), y=outlier_Y), size=3, shape=23, fill="red") +  # not fitting well
theme_bw() +
scale_y_continuous(lim=c(-10, 90), name="output, f(x)") +
xlab("input, x")


# par(mfrow=c(1,2))
# significance = (Y - f_pred) / sqrt(diag(fvar_pred))
# ggplot(NULL, aes(x = x[,1], y = x[,2], color = significance, size=5)) + 
# geom_point() + 
# labs(x = "modulator 1") + 
# labs(y = "modulator 2")

# ggplot(NULL, aes(x = x[,1], y = x[,2], color = Y, size=5)) + 
# geom_point() + 
# labs(x = "modulator 1") + 
# labs(y = "modulator 2")

# iff = Y > 50
# ggplot(NULL, aes(x = x[,1][iff], y = x[,2][iff], color = Y[iff], size=5)) + 
# geom_point() + 
# labs(x = "modulator 1") + 
# labs(y = "modulator 2")

X_pred = expand.grid(seq(10, 30, len=50), seq(5, 10, len=50)) 

for (i in 1: ncol(X)) X_pred[, i] = (X_pred[, i] - mean(x[, i])) / sd(x[, i])

var_obs = rbfKernel(X, X, l, sigma_f)
Kconv = rbfKernel(X, X_pred, l, sigma_f)
var_pred = rbfKernel(X_pred, X_pred, l, sigma_f)


L = chol(var_obs + sigma_n^2 * diag(1, ncol(var_obs)))
invL = solve(L)
alpha = invL %*% t(invL) %*% (Y - mean(Y))
v_ = t(invL) %*% Kconv
f_pred = as.vector(mean(Y) + t(Kconv) %*% alpha)
fvar_pred = var_pred - t(v_) %*% v_

plot(f_pred)

df_pred = data.frame(x1=X_pred[, 1] * sd(x[, 1]) + mean(x[, 1]),
                     x2=X_pred[, 2] * sd(x[, 2]) + mean(x[, 2]),
                     y=f_pred)

df_pred[which.max(df_pred[,3]), ]
scatterplot3d(df_pred$x1, df_pred$x2, f_pred, color="lightblue",pch=21, xlab="modulator 1", ylab="modulator 2", zlab="", main="Posterior predictive distribution")
