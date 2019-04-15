# A simple demo of Bayesian Optimization for Model Selection on Elastic Nets
# Author: Alex Fout
# Date: April 13 2019

############################
########## Setup ###########
############################
rm(list=ls())
library(mvtnorm)
library(glmnet)
library(nloptr)
library(ggplot2)

### Generate Data
get_data <- function(){
  set.seed(utf8ToInt("my other prior is a posterior"))
  n <- 1000
  p <- 1000
  r <- 100
  x <- matrix(runif(n*p)<0.5, n, p)
  b <- matrix(0, p)
  b[sample(p, r)] <- 4*(runif(r)-0.5)
  pi <- 1/(1+exp(-x %*% b))
  y <- runif(n)<pi
  return(list(x=x,y=y,b=b))
}
data <- get_data()
print(sum(data$y)/length(data$y))  # incidence rate


### Objective function (fit an elastic net and compute mse)
get_g <- function(data){
  g <- function(hp){
    # hp: matrix of hyperparameters for elastic net
    # each row in hp is a set of hp's (lambda, alpha)
    # where lambda >0 is the shrinkage magnitude and 
    # alpha in (0,1) is the balance between L1 and L2 penalties
    mse <- array(NA, c(dim(hp)[1], 1))
    for(i in 1:dim(hp)[1]){
      m <- glmnet(data$x, data$y, family = "binomial", lambda=c(hp[i,1]) , alpha=hp[i,2], standardize=F)
      mse[i,] <- sqrt(sum((data$y - predict(m, data$x))^2))
    }
    return(-mse)
  }
}
g <- get_g(data)

### Gaussian Process Components
get_kernel <- function(a){
  # a: list of hyperparameters for the kernel
  # a[1] is coefficient outside exponent
  # a[2:end] are weights for each parameter in the model
  a0 <- a[1]
  a <- a[2:length(a)]
  kernel <- function(x1,x2){
    # compute the covariance matrix between two sets of d-dimensional points
    # x1 is n1 by d array of n1 points in d dimensions
    # x2 is n2 by d array of n2 points in d dimensions
    n1 <- dim(x1)[1]
    # if only set of points, find cov with itself
    if(missing(x2)){
      K <- array(NA, c(n1,n1))
      for(i in 1:n1){
        K[i,i] <- a0
        for(j in 1:(i-1)){
          K[j,i] <- K[i,j] <- a0*exp(-sum(a^2*(x1[i,]-x1[j,])^2))
        }
      }
    } else {
      n2 <- dim(x2)[1]
      K <- array(NA, c(n1,n2))
      for(i in 1:n1){
        for(j in 1:n2){
          K[i,j] <- a0*exp(-sum(a^2*(x1[i,]-x2[j,])^2))
        }
      }
    }
    return(K)
  }
  return(kernel)
}
get_mean_fn <- function(a0){
  # get a constant mean function
  mu <- function(x){
    array(a0, c(dim(x)[1], 1))
  }
  return(mu)
}
minus_LL <- function(hp, knots, get_kernel, get_mean_fn){
  # return negative LL of gaussian process at knots
  kernel_hp <- hp[1:3]
  mean_fn_hp <- hp[4]
  kernel <- get_kernel(kernel_hp)
  mu <- get_mean_fn(mean_fn_hp)
  muk <- mu(knots)
  sigmak <- kernel(knots, knots)
  log_f_prior <- dmvnorm(gknots[,1], muk, sigmak, log=T)
  return(-log_f_prior)
}
f_post <- function(x, knots, gknots, sigma0kk){
  # compute posterior of GP given function value at knots
  # x: lattice (or whatever) of points to compute the posterior
  # knots: n by d matrix of n d-dimensional points
  # gknots: n by 1 matrix of function values at the knots
  # sigma0kk: optional covariance matrix for knots
  mux <- mu(x)
  if(missing(sigma0kk)){
    sigma0kk <- kernel(knots) + diag(10^-6, dim(knots)[1])
  }
  sigma0kk_inv <- chol2inv(chol(sigma0kk))
  mun <- array(NA, c(dim(x)[1], 1))
  sigman <- array(NA, c(dim(x)[1], 1))
  for(i in 1:dim(x)[1]){
    xi <- x[i,,drop=F]
    sigma0ik <- kernel(xi, knots)
    sigma0ii <- kernel(xi)
    sigma0ikkinv <- sigma0ik %*% sigma0kk_inv
    mun[i] <- sigma0ikkinv %*%(gknots-mu(knots)) + mu(xi)
    sigman[i] <- sqrt(max(sigma0ii - sigma0ikkinv %*% t(sigma0ik), 0))
  }
  return(list(mun=mun, sigman=sigman, sigma0kk=sigma0kk))
}
### Bayesian Optimization functions
# Expected Improvement (Acquisition Function)
EI <- function(x, knots, gknots, sigma0kk){
  # Evaluate aqcuistion function at points x
  fp <- f_post(x, knots, gknots, sigma0kk)
  fstar <- max(gknots)
  deln <- fp$mun-fstar
  EI <- pmax(deln, 0) + fp$sigman*dnorm(deln/fp$sigman)-abs(deln)*pnorm(-abs(deln)/fp$sigman)
  return(EI)
}
# Maximize the Aqcuisition function
max_EI <- function(knots, gknots, sigma0kk){
    nloptr(colMeans(knots), function(y) -EI(array(y, c(1,2)), knots, gknots, sigma0kk), lb=c(0,0), ub=c(1,1), opts=list("algorithm"="NLOPT_GN_CRS2_LM","xtol_rel"=1e-5,"maxeval"=2000))
}


############################
### Bayesian Optimization ##
############################

### Initial Knots
knots <- array(c(0.1, 0.7, 0.3, 0.2, 0.5, 0.9), c(3,2))
#knots <- array(runif(20), c(10,2))
gknots <- g(knots)  # evaluate g at knots
#plot(knots, xlim=c(0,1), ylim=c(0,1))

### MLE to fit Gaussian Process hyperparameters
hp_MLE <- optim(c(1,1,1,1), minus_LL, gr=NULL, knots, get_kernel, get_mean_fn, method="Nelder-Mead")$par
print(hp_MLE)
# Let's cheat a little by adjusting the hyperparameters from the MLE
hp_MLE[2:3] <- hp_MLE[2:3]/2
hp_MLE[1] <- hp_MLE[1]*2
print(hp_MLE)
# get kernel and mean fn for MLE hyperparameters
kernel <- get_kernel(hp_MLE[1:3])
mu <- get_mean_fn(hp_MLE[4])

### Do the Bayesian Optimization
sigma0kk <- kernel(knots) + diag(1e-6, dim(knots)[1])
n_iters <- 150
for(i in 1:n_iters){
  if(n_iters >=10 & (i %% round(n_iters/10)==0)){print(i)}
  EI_max <- max_EI(knots, gknots, sigma0kk)
  knot_new <- array(EI_max$solution, c(1,2))
  gnew <- g(knot_new)
  sigma0ik <- kernel(knot_new, knots)
  sigma0ii <- kernel(knot_new)
  sigma0kk <- rbind(cbind(sigma0kk, t(sigma0ik)),cbind(sigma0ik, sigma0ii))+diag(1e-6, dim(sigma0kk)[1]+1)
  knots <- rbind(knots, knot_new)
  gknots <- rbind(gknots, gnew)
}


############################
###### Model Assessment ####
############################
gbest <- cummax(gknots)
plot(gbest, type="l")
# get coeff estimates for the best model at each iter
nonzero_cfs <- array(NA, c(dim(knots)[1]))
cfs <- matrix(NA, dim(knots)[1], dim(data$x)[2]+1)
best_hp <- matrix(NA, dim(knots)[1], 2)
for(i in 3:dim(knots)[1]){
  best_hp[i,] <- knots[which.max(gknots[1:i]),]
  m <- glmnet(data$x, data$y, family = "binomial", lambda=c(best_hp[i,1]) , alpha=best_hp[i,2], standardize=F)
  cfs[i,] <- coef(m)[,1]
  nonzero_cfs[i] <- sum(cfs[i,]!=0)
}


# lattice of points to use for heatmaps
xres <- 100
x1 <- seq(0, 1, 1/(xres-1))
x2 <- seq(0, 1, 1/(xres-1))
x <- cbind(rep(x1, length(x2)), rep(x2, each=length(x1)))

# true value of function
greal <- g(x)


############################
###### Visualization #######
############################
# scaling functions so we can use same scale for all three heatmaps
range_mu <- function(x){
  x <- pmin(pmax(x, -20), -17.55)
  ((x+20)/(20-17.55))
}
range_sigma <- function(x){
  x <- x^2
  x/0.6^2
}
range_acq <- function(x){
  hist(acq, breaks=100)
  x <- pmin(20, x)
  x/20
}

# Plot gbest vs iter
png("graphics\\best_mse.png", width=6, height=6, units="in", res=300)
plot(-gbest, type="l", lwd=2, main="Best MSE So Far", xlab="Iteration", ylab="MSE")
dev.off()

# plot coefficients vs. iter
png("graphics\\gwas_coeficient_matplot.png", width=6, height=6, units="in", res=300)
matplot(cfs[1:50,2:1001], type="l", main="Regression Coefficient Estimates", ylab="Estimate", xlab="Iteration")
dev.off()
# plot num nonzero coeffs vs iter
png("graphics\\gwas_nonzero_coeficient.png", width=6, height=6, units="in", res=300)
plot(nonzero_cfs, type="l", lwd=2, main="Number of Nonzero Coefficients", xlab="Iteration", ylab="Count")
dev.off()
# plot percentiles vs iter
tiles <- c(0.1, 0.5, 0.9, 0.99)
qtiles <- apply(abs(cfs[1:50,]), 1, function(x) quantile(x, tiles, na.rm=T))
png("graphics\\gwas_coef_mag_quantiles.png", width=6, height=6, units="in", res=300)
matplot(t(qtiles), type="o", lty=1, pch=20, cex=0.8, main="Coefficient Magnitudes (Quantiles)", ylab="Coefficient Magnitude", xlab="Iteration")
legend(x="topright", legend = paste0(100*rev(tiles), "th %"), lty=1, pch=20, col=4:1, bty="n")
dev.off()
# true vs abs coefficients
for(i in 3:dim(cfs)[1]){
  png(sprintf("graphics\\gwas_animation2\\coef%03d.png", i-3), width=6, height=6, units="in", res=300)
  plot(data$b, cfs[i,2:1001], main=paste0("True vs. Estimated Coefficients, i=", i-3), xlab="True Value", ylab="Estimate", xlim=c(-2, 2), ylim=c(-0.07, 0.07), pch=20, cex=0.7)
  dev.off()
}
# model parameters vs. iter
cols <- c(rgb(1,0,0), rgb(0,0.7,0))
png("graphics\\gwas_en_hyperparams.png", width=6, height=6, units="in", res=300)
matplot(best_hp[1:50,], type="o", lty=1, pch=20, cex=0.8, col=cols, main="Elastic Net Hyperparameters", xlab="Iteration", ylab="Value")
legend(x="topleft", legend = c(expression(lambda), expression(alpha)), lty=1, pch=20, col=cols, bty="n")
dev.off()


# look at the knots
plot(knots, xlim=c(0,1), ylim=c(0,1))
p <- ggplot(data.frame(lambda=x[,1], alpha=x[,2], g=range_mu(greal))) + geom_tile(aes(lambda,alpha,alpha=g), fill="white", show.legend=F) + 
  geom_point(data=as.data.frame(knots), aes(V1,V2), color="green") + 
  labs(title="True Objective Function (With Knots)", x=expression(lambda), y=expression(alpha))+
  theme(panel.background=element_rect(fill="black"), panel.grid=element_blank()) + 
  theme(axis.title.y = element_text(angle = 0, vjust=0.5), strip.background =element_rect(fill="white")) + 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  lims(x=range(x[,1]), y=range(x[,2]))
ggsave(filename=paste0("graphics\\true_g.png"), 
       plot=p, width=6, height=6, units="in")
  


## make a heatmap for each iteration, with mean fn, sigma, and acq fn side by side
nx <- dim(x)[1]
for(i in 3:dim(knots)[1]){
  # compute posterior and acquisition fn at this iteration
  iknots <- as.matrix(knots[1:i,,drop=F])
  igknots <- as.matrix(gknots[1:i,drop=F])
  sigma0kk <- kernel(iknots) + diag(1e-6, dim(iknots)[1])
  fp <- f_post(x, iknots, igknots, sigma0kk)
  acq <- EI(x, iknots, igknots)
  df_plot <- data.frame(var=rep(c(1,2,3), each=nx), 
                        lambda=rep(x[,1],3), 
                        alpha=rep(x[,2],3), 
                        val=c(range_mu(fp$mun), range_sigma(fp$sigman), range01(acq)))
  df_plot$var <- factor(df_plot$var, levels=c(1,2,3), labels=c("mu", "sigma", "Acquisition~Function"))
  # Plot heatmaps
  p <- ggplot(df_plot) + geom_tile(aes(lambda,alpha,fill=var,alpha=val), show.legend=F) + 
    facet_grid(.~var, labeller=label_parsed) + 
    geom_point(data=as.data.frame(iknots), aes(V1,V2), color="green", size=1.6, show.legend=F) + 
    scale_fill_manual(values=c("red", "blue", "darkgreen")) + 
    scale_alpha_continuous(range=c(0,1)) + 
    theme(panel.background=element_rect(fill="black"), panel.grid=element_blank()) + 
    theme(axis.title.y = element_text(angle = 0, vjust=0.5), strip.background =element_rect(fill="white")) + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + 
    lims(x=range(x[,1]), y=range(x[,2])) + 
    labs(x=expression(lambda), y=expression(alpha), title=paste0("Elastic Net Model Selection: GWAS Example, iter=", i-3))
  #print(p)
  ggsave(filename=paste0("graphics\\gwas_animation\\iter", sprintf("%03d", i-3), ".png"), 
         plot=p, width=12, height=4, units="in")
  print(i-3)
}

