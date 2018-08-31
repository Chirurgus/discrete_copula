# Created by Oleksandr Sorochynskyi
# On 14 06 18

source("n-copula.R")
source("copula.R")


#######################################################################################
# Test speed, and maybe optimize
library(profvis)

nobs <- 10000;
dim <- 3;

x <- neurons[1:dim,1:nobs]

profvis({
  system.time(fit <- ncopula_fit(t(x),normalCopula(dim=dim,dispstr = "un")));
})

# results
# for 100 obs
# dim => time (sec)
# 2 => 0.4
# 3 => 27.9
# 4 => 106.9
# 5 => 371.06

# for 500 obs
# dim => time (sec)
# 2 => 0.8
# 3 => 39.7
# 4 => 
# 5 => 

# for 10000 obs
# dim => time (sec)
# 2 => 5
# 3 => 178  (655 for 44000 obs)
# 4 => 
# 5 => 

# So I don't see a way to optimize it since basically all the time is spend in prob
# and pCopula.
# One possibility would be to try to paralellize this, but I think it will only
# be useful if we can estimate multiple PAIRS at the time, since i suspect
# paralellizing loglikelihood will induce more overhead than speedup.
# Also maybe move mvnorm to GPU ??? :)

#######################################################################################
# Debug ncopula by comparit it to 2 dimentional case

t1 <- neurons[3,1:100]
t2 <- neurons[4,1:100]
cdf1 <- ecdf(t1);
cdf2 <- ecdf(t2);

# fit_copula works fine
# ncopula_fit diverges
(fit <- fit_copula(t1,t2,normalCopula(dim=2,dispstr = "un")));
(fit <- ncopula_fit(cbind(t1,t2),normalCopula(dim=2,dispstr = "un")));

# maybe give ncopula_fit a good inital value
# ... if i give it exacly the target value it does fine, otherwise it diverges to the
#     clostest bound, ie +/-1
# All the while fit_copula converges no matter teh inital value
(fit <- ncopula_fit(cbind(t1,t2),normalCopula(dim=2,dispstr = "un"), inital_val = -8));


# "set" contains row indices where likelihoods are negative or null
# "param" contains the parameter given to loglik that caused negative likelihoods

#usan all the data
cop_loglik <- copula_log_lik_maker(t1,t2,cdf1,cdf2,normalCopula(dim=2,dispstr = "un"));
ncop_loglik <- ncopula_make_loglik(cbind(t1,t2),c(cdf1,cdf2),normalCopula(dim=2,dispstr = "un"))

# using the troublsome subset
ts1 <- t1[set]
ts2 <- t2[set]
cop_loglik <- copula_log_lik_maker(ts1,ts2,cdf1,cdf2,normalCopula(dim=2,dispstr = "un"));
ncop_loglik <- ncopula_make_loglik(cbind(ts1,ts2),c(cdf1,cdf2),normalCopula(dim=2,dispstr = "un"))

cop_loglik(param)
ncop_loglik(param)
#######################################################################################
# Compare my combination version vs simple numerical integration
dim <- 2
nobs <- 10

x <- neurons[1:dim,1:nobs]
x <- matrix(rpois(dim*nobs,0.5),ncol=dim,nrow=nobs)
marginals <- apply(x,2,ecdf)
cop <- normalCopula(rep(0.5,max(dim*(dim-1)/2,0)),dim=dim,dispstr="un");

pseudo_density(x[1,1],x[1,2],marginals[[1]],marginals[[2]],cop)
generalized_pseudo_density(x[1,],marginals,cop)
generalized_pseudo_density_ingegration(x[1,],marginals,cop)$integral


#######################################################################################
# Debug negative probabilites

t1 <- neurons[3,]
t2 <- neurons[4,]
t3 <- neurons[5,]
marginals <- c(ecdf(t1),ecdf(t2),ecdf(t3))
# obtained while debuging
# set is the subset of neurons where this happpned
set <- 2
par <- c(1,-1,1)
x <- matrix(cbind(t1,t2,t3)[set,],ncol=3)

cop <- normalCopula(par, 3, dispstr = "un");
loglik <- ncopula_make_loglik(x,marginals,cop);
loglik(par);

oppd <- make_optim_generalized_pseudo_density(4,3,cop,marginals)

# While in loglik this gives a negative number, here it plays (somewhat nicely)
# it just gives 0, this is coherent with the result from loglik, but troubling for
# the said loglikelihood as log(0) = -Inf

oppd(as.vector(x))
generalized_pseudo_density(as.vector(x),marginals,cop)
generalized_pseudo_density_ingegration(as.vector(x),marginals,cop)

# So the problem must be deeper
u <- mapply(function(x,f) f(x), x=x, f=marginals)
u_minus <- mapply(function(x,f) f(x), x=(x-1), f=marginals)

pd <- function(u) {
  pCopula(u, cop);
}

get <- function(a,b,c) {
  ret <- u_minus
  if (a == 1) { ret[1] = u[1] }
  if (b == 1) { ret[2] = u[2] }
  if (c == 1) { ret[3] = u[3] }
  ret
}

pd(get(0,0,0))
pd(get(0,0,1))
pd(get(0,1,0))
pd(get(0,1,1))
pd(get(1,0,0))
pd(get(1,0,1))
pd(get(1,0,0))
pd(get(1,1,1))
  
##########################################################################################
# Test on neuron data

t1 <- neurons[3,1:100]
t2 <- neurons[4,1:100]
t3 <- neurons[5,1:100]
#t4 <- neurons[6,1:100]

fit <- ncopula_fit(cbind(t1,t2,t3), normalCopula(dim=3,dispstr = "un"))
fit <- ncopula_fit(cbind(t1,t2,t3), tCopula(dim=3,dispstr = "un",df.fixed=TRUE))
fit <- ncopula_fit(cbind(t1,t2,t3), frankCopula(dim=3),inital_val = 5)

nfit <- ncopula_fit(cbind(t1,t2), gumbelCopula(),hessian = TRUE );
fit <- fit_copula(t1,t2,gumbelCopula())
##########################################################################################
# Test 3 dim case

u <- matrix(runif(3*100),ncol=3)

fit <- ncopula_fit(u,normalCopula(dim=3))

cdf1 <- ecdf(u[,1]);
cdf2 <- ecdf(u[,2]);
cdf3 <- ecdf(u[,3]);

loglik <- ncopula_make_loglik(u, c(cdf1,cdf2,cdf3),normalCopula(3))
loglik(P2p(cov(u)))
##########################################################################################
# Test generalized_pseudo_desnity
sig <- matrix(c(1,0.1,0.1,1),nrow=2);
cop <- normalCopula(P2p(sig),2);
u <- cbind(c(0,0,1,2,3),c(0,1,2,0,3));
cdf1 <- ecdf(u[,1]);
cdf2 <- ecdf(u[,2]);
pseudo_density(u[1,1],u[1,2],cdf1,cdf2,cop);
generalized_pseudo_density(u[1,],c(cdf1,cdf2),cop)

pd <- make_pseudo_density(max(u),cop,cdf1,cdf2);
gpd <- make_optim_generalized_pseudo_density(2,2,cop,c(cdf1,cdf2))
pd(u[1,1],u[1,2])
gpd(u[1,])


loglik <- ncopula_make_loglik(u, c(cdf1,cdf2),cop)
loglik(P2p(sig))

fit <- ncopula_fit(u, normalCopula(),c(cdf1,cdf2),lower=-0.99999,upper=0.9)

# Plot loglike as it give non-finite difference
x <- seq(-1,1,by=0.01,)
plot(x,lapply(x,function(x) -loglik(x)))


f <- function(x) {
  if(x >= 0) {
    return(x^2)
  }
  else {
    return(-Inf)
  }
}
library(optimx)
fit <- optimx(10, f, lower=0, upper=20,hessian=TRUE)
