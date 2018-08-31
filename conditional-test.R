# Created by Oleksandr Sorochynskyi
# On 19 07 18

source('~/IdV/R/copula.R')


################################################################################
# Different mean for data generated form unifrom distributions

data <- as.vector(bin_checkerboard(1))

cdf <- ecdf(data);
qu <- make_quantile(data);

all(data == qu(cdf(data)))

mean(data);
mean(estim_data);

plot(function(x) qu(x),xlim=c(0,1))

################################################################################
# Estimate noise_cov variance
i <- 23
j <- 12
i <- pair_near[1];
j <- pair_near[2];

m <- bin_checkerboard(c(j,i))
load("C:/Documents/IdV/R/params/conditional_1_25_off_check_frank_nonzero_mean.RData")

cond_noise_cov(1,2,m)
n.rep <- 1000
noise_cov <- replicate(200,{
  rand <- cond_model2.2_rand(1, 2, m, model2.2$param[i,j], n.rep);
  
  # noise_cov_real[i,j] <- cond_noise_cov(i,j+n0, m);
  # stimul_cov_real[i,j] <- cond_stimul_cov(i, n0+j, m);
  
  cond_noise_cov(1,2,rand);
  #stimul_cov_simu[i,j] <- cond_stimul_cov(1,2,rand);
});

################################################################################

# Take another look at independence tests
m <- bin_barmovie(pair_median);
x1 <- as.vector(m[1,,]);
x2 <- as.vector(m[2,,]);
copula <-frankCopula(model3_bar[pair_median[1],pair_median[2]])

copula_kramer_indep_test(x1,x2,nrep=10);
copula_pillow_test(x1, x2, copula);
copula_loglik_test(x1,x2, indepCopula(), copula, one.sided = TRUE);

################################################################################
# Noise cov of generated data vs real
 
n1 <- 1;
n2 <- 2;

bar_m <- bin_barmovie(c(n1,n2));
full_m <- bin_fullfield(c(n1,n2));
check_m <- bin_checkerboard(c(n1,n2));

# Fit on checkerboard data
selected.bins <- cond_select_bins(n1, n2, bar_m, cond_mean_prod_criteria);
copula <- frankCopula();
fit <- cond_model2.2(n1,n2,bar_m,copula,selected.bins);

rand <- cond_model2.2_rand(n1, n2, bar_m, fit$param, 1000);

cond_total_cov(n1,n2,bar_m);
cond_noise_cov(n1,n2,bar_m)
cond_stimul_cov(n1, n2, bar_m)

cond_total_cov(n1,n2,rand);
cond_noise_cov(n1,n2,rand);
cond_stimul_cov(n1, n2,rand);


################################################################################
# Try nested model
data <- bin_neurons(c(1,2), 4, 0);
t1 <- data[1,];
t2 <- data[2,];

stimul <- runif(length(t1));
stimul_pobs <- ecdf(stimul)(stimul);

inner_cop <- fit_copula(t1,t2, claytonCopula())$cop;
inner_param <- inner_cop@parameters;


inner_pobs <- pCopula(cbind(ecdf(t1)(t1),ecdf(t2)(t2)),inner_cop);
outter_cop <- fitCopula(claytonCopula(),cbind(stimul_pobs,inner_pobs))@copula;
outter_param <- outter_cop@parameters;

cop <- onacopula("Clayton", C(outter_param, 1, C(inner_param, 2:3)))

################################################################################
# Normalize data by removing the mean of the bin from the data and fit.
# The result should give the dependance structure due to "noise",
# and not due to stimulus.
m <- get_rep_neurons();
n1 <- 8;
n2 <- 16


std_x <- sapply(1:dim(m)[2], function(t) {
  m[n1,t,] - mean(m[n1,t,]);
});
std_y <- sapply(1:dim(m)[2], function(t) {
  m[n2,t,] - mean(m[n2,t,]);
});

x <- as.vector(std_x);
y <- as.vector(std_y);

copula <- normalCopula();
fit <- fit_copula(x, y, copula);

copula_plot(x, y);
copula_plot_cop_density(fit$cop, x, y);
copula_plot_diff(x, y, fit$cop);

# Now let's see if this model can generate intresting statistics
# Correlation
rand <- copula_rand_emp(length(x), fit$cop, x1=x, x2=y);
r_x <- rand[,1];
r_y <- rand[,2];
copula_plot(r_x,r_y);

r_m <- array(NA, dim=c(2, dim(m)[-1]));
r_m[1,,] <- r_x;
r_m[2,,] <- r_y;
cond_noise_cov(1, 2, r_m);
cond_stimul_cov(1, 2, r_m);
cond_total_cov(1, 2, r_m);
# VS
n_m <- array(NA, dim=c(2, dim(m)[-1]));
n_m[1,,] <- x;
n_m[2,,] <- y;
cond_noise_cov(n1, n2, n_m);
cond_stimul_cov(n1, n2, n_m);
cond_total_cov(n1, n2, n_m);

# This basically failed, although both stimul cov are basically 0,
# the noise covariance is quite different

# Condional on one neuron
plot(x,y)
copula_plot(x,y);
n <- 10000

u1 <- 0.05
rand <- copula_rand_cond(rep(u1,n),runif(n),cc,x,y);
hist(rand[,2], main=as.character(u1),xlim=c(-2,2));

u2 <- 0.5
rand <- copula_rand_cond(rep(u2,n),runif(n),cc,x,y);
hist(rand[,2], main=as.character(u2),xlim=c(-2,2));

u3 <- 0.95
rand <- copula_rand_cond(rep(u3,n),runif(n),cc,x,y);
hist(rand[,2], main=as.character(u3),xlim=c(-2,2));

# Unsurpising results
################################################################################
# Compare tau infered from bilinear vs discrete coopulas
data <- bin_neurons(c(1,2), 4,0);
dfit <- fit_copula(data[1,],data[2,],normalCopula())
cfit <- fitCopula(frankCopula(), cbind(copula_distributional_transfrom(data[1,]),copula_distributional_transfrom(data[2,])))


################################################################################
# Does copula theta vary over timebins
#   => Theata does seem to depend on time, howerver the many parameters are close to 0
#      Anoter considerable set of them are negative

m <- get_rep_neurons();
# these are 129 whatever appart
n1 <- 6;
n2 <- 11;
# these are 1080 whatever appart
n1 <- 8;
n2 <- 16

copula <- normalCopula();
pars <- sapply(1:dim(m)[2], function(t) {
                              # We don't want to fit all time bins, only those that fire.
                              if (mean("*"(m[n1,t,],m[n2,t,])) == 0
                                  #|| abs(cor(m[n1,t,],m[n2,t,],method="s")) < 0.2
                                  #|| !copula_kramer_indep_test(m[n1,t,],m[n2,t,])$signif
                              ) {
                                NA;
                              }
                              else {
                               fit_copula(m[n1,t,], m[n2,t,], copula)$par
                              }
})
# like half of these are NA, WITHOUT filtering out independant ones
plot(pars,ylim=c(-1,1))

# Let's see how "FERRARI COEF" influences the parameter
ferrari_coefs <- sapply(1:dim(m)[2], function(t) {
  sqrt("*"(mean(m[n1,t,]),mean(m[n2,t,])));
});

# How close are the two means?
mean_nspikes1 <- sapply(1:dim(m)[2], function(t) {
  mean(m[n1,t,]);
});

mean_nspikes2 <- sapply(1:dim(m)[2], function(t) {
  mean(m[n2,t,]);
});

cor(mean_nspikes1, mean_nspikes2)

mean(sapply(1:dim(m)[2], function(t) {
  cor(m[n1,t,] - mean_nspikes1,m[n2,t,] - mean_nspikes2)
}));

# cov between 2 neurns for a given timebin
covar <- sapply(1:dim(m)[2], function(t) {
  cov(m[n1,t,],m[n2,t,]);
});

plot(mean_nspikes,covar);

################################################################################
# "Scatterplot" of num spikes of one neuron vs another
n1 <- 1;
n2 <- 2;
bin <- 755;
count <- xyTable(m[n1,bin,],m[n2,bin,]);
plot(count$x, count$y, cex=count$number)
################################################################################
# noise cov vs distance => got a nice 1/x behavior

m <- get_rep_neurons();

noise_cov <- function(x,y) {
  mean(sapply(1:dim(m)[2], function(i) cov(m[x,i,],m[y,i,])));
}

m.noise.cov <- outer(1:25, 1:25, function(x,y) { mapply(noise_cov, x, y); })
plot(neurons_distance[lower.tri(neurons_distance)], m.noise.cov[lower.tri(m.noise.cov)],
     ylab="noise cov",xlab="distance")
