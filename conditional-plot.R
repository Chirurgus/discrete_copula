# Created by Oleksandr Sorochynskyi
# On 02 08 18

library(ggplot2)

################################################################################
# Noise cov /time bin, in data VS model

cov_over_time <- function(m) {
  sapply(1:dim(m)[2], function(t) { cov(m[1,t,],m[2,t,]); });
}
n1 <- 8;
n2 <- 16;
data <- bin_checkerboard(c(n1,n2));
bins <- cond_select_bins(1,2,data,cond_mean_prod_criteria);
fit <- cond_model2.2(1,2,data,frankCopula(),bins);
simu <- cond_model2.2_rand(1,2,data,fit$param,200)

cov_data <- cov_over_time(data);
cov_simu <- cov_over_time(simu);

plot(cov_data, cov_simu, main="Dist = 1050");
abline(h=mean(cov_simu))
abline(v=mean(cov_data))
abline(0,1)

################################################################################
# Model 3 vs model 2.2 param (these should be very simular for uncorrelated stimuli),
# but with a constant bias for those with correlated stimuli (fullfield)

filename <- "conditional_1_25_off_fullfield_frank_nonzero_mean.RData";
load(filename)
model2.2_fullfield <- model2.2$param;
model3_fullfield <- model3$param;

filename <- "conditional_1_25_off_check_frank_nonzero_mean.RData";
load(filename)
model2.2_check <- model2.2$param;
model3_check <- model3$param;

# This should be almost identical
plot(model2.2_check, model3_check);
abline(0,1);

# These should be biased
plot(model2.2_fullfield, model3_fullfield);
abline(0,1);#

################################################################################
# noise/stimul cov in data vs simulation for the 3rd model

n <- 10;

bar_m <- bin_barmovie(1:n);
full_m <- bin_fullfield(1:n);
check_m <- bin_checkerboard(1:n);

copula <- frankCopula();

noise_cov_real <- matrix(NA, ncol=n, nrow=n);
noise_cov_simu <- matrix(NA, ncol=n, nrow=n);
stimul_cov_real <- matrix(NA, ncol=n, nrow=n);
stimul_cov_simu <- matrix(NA, ncol=n, nrow=n);

for (i in 2:n) {
  for (j in 1:(i-1)) {
    # Fit on checkerboard data
    selected.bins <- cond_select_bins(i, j, bar_m, cond_mean_prod_criteria);
    fit <- cond_model3(i,j,check_m,copula,selected.bins);
    
    rand_raw <- copula_rand_emp(200*dim(bar_m)[2],
                                frankCopula(fit$param),
                                x1=as.vector(bar_m[i,,]),
                                x2= as.vector(bar_m[j,,]));
    rand <- array(NA, dim=c(2, dim(bar_m)[2], 200));
    rand[1,,] <- rand_raw[,1];
    rand[2,,] <- rand_raw[,2];
    
    noise_cov_real[i,j] <- cond_noise_cov(i,j,bar_m);
    stimul_cov_real[i,j] <- cond_stimul_cov(i, j, bar_m);
    
    noise_cov_simu[i,j] <- cond_noise_cov(1,2,rand);
    stimul_cov_simu[i,j] <- cond_stimul_cov(1,2,rand);
  }
}

plot(noise_cov_real, noise_cov_simu);
abline(0,1);

plot(stimul_cov_real, stimul_cov_simu);
abline(0,1);

################################################################################
# noise/stimul cov in data vs simulation

n <- 10;

bar_m <- bin_barmovie(1:n);
full_m <- bin_fullfield(1:n);
check_m <- bin_checkerboard(1:n);

copula <- frankCopula();
  
noise_cov_real <- matrix(NA, ncol=n, nrow=n);
noise_cov_simu <- matrix(NA, ncol=n, nrow=n);
stimul_cov_real <- matrix(NA, ncol=n, nrow=n);
stimul_cov_simu <- matrix(NA, ncol=n, nrow=n);

for (i in 2:n) {
  for (j in 1:(i-1)) {
    # Fit on checkerboard data
    selected.bins <- cond_select_bins(i, j, bar_m, cond_mean_prod_criteria);
    fit <- cond_model2.2(i,j,check_m,copula,selected.bins);
    
    rand <- cond_model2.2_rand(i, j, bar_m, fit$param, 200);
    
    noise_cov_real[i,j] <- cond_noise_cov(i,j,bar_m);
    stimul_cov_real[i,j] <- cond_stimul_cov(i, j, bar_m);
    
    noise_cov_simu[i,j] <- cond_noise_cov(1,2,rand);
    stimul_cov_simu[i,j] <- cond_stimul_cov(1,2,rand);
  }
}

plot(noise_cov_real, noise_cov_simu);
abline(0,1);

plot(stimul_cov_real, stimul_cov_simu);
abline(0,1);
################################################################################
# compare coefficients from 2 different stimuli (Bar vs ...)

filename <- "conditional_1_10_off_bar_frank_nonzero_mean.RData";
load(filename)
model2.2_bar <- model2.2$param;

filename <- "conditional_1_17_off_fullfield_frank_nonzero_mean.RData";
load(filename)
model2.2_fullfield <- model2.2$param[1:10, 1:10];

filename <- "conditional_1_17_off_check_frank_nonzero_mean.RData";
load(filename)
model2.2_check <- model2.2$param[1:10, 1:10];

plot(model2.2_fullfield, model2.2_bar)
abline(0, 1);

plot(model2.2_bar, model2.2_check)
abline(0,1);

plot(model2.2_fullfield, model2.2_check);
abline(0,1);

################################################################################
# Conditional model theta vs distance

cond_plot_distance <- function(theta, dist, ...) {
  plot(dist, theta, xlab= "Distance (micro m)", ylab= "Theta", ylim= c(-10,10), ...);
  abline(h=0);
}

cond_plot_meanprod <- function(theta, meanprod, other_theta, size, ...) {
  df <- data.frame(theta=theta,meanprod=meanprod,size=size);
  ggplot(data= df, mapping = aes(x=meanprod, y= theta)) +
    geom_point(aes(colour = factor(size))) +
    geom_abline(slope=0, intercept= other_theta) +
    annotate(geom="text", label="theta_2.2", x=0.9, y=other_theta, vjust=-1) +
    geom_abline(slope=0, intercept= mean(theta,na.rm=TRUE)) +
    annotate(geom="text", label="mean(theta,na.rm=TRUE)", x=0.9, y=mean(theta,na.rm=TRUE), vjust=-1) +
    ylab("Theta(t)") +
    xlab("sqrt(mean(x)*mean(y))");
}

filename <- "conditional_1_8_off_check_frank_nonzero_mean.RData";
filename <- "conditional_1_10_off_bar_frank_nonzero_mean.RData";
filename <- "conditional_1_17_off_check_frank_nonzero_mean.RData";
filename <- "conditional_1_17_off_frank_nonzero_mean.RData";
filename <- "conditional_1_17_off_frank_nonzero_jointspike.RData";
load(filename);
n <- 10;

# Coefs VS distance
coefs <- apply(model1$param, c(1,2), function(x) mean(x, na.rm= TRUE));
cond_plot_distance(coefs[lower.tri(coefs)],
                   neurons_distance[1:n,1:n][lower.tri(neurons_distance[1:n,1:n])],
                   main= "Model 1");

coefs <- model2.1$param;
cond_plot_distance(coefs[lower.tri(coefs)],
                   neurons_distance[1:n,1:n][lower.tri(neurons_distance[1:n,1:n])],
                   main= "Model 2.1");


coefs <- model2.2$param;
cond_plot_distance(coefs[lower.tri(coefs)],
                   neurons_distance[1:n,1:n][lower.tri(neurons_distance[1:n,1:n])],
                   main= "Model 2.2");

coefs <- model3$param;
cond_plot_distance(coefs[lower.tri(coefs)],
                   neurons_distance[1:n,1:n][lower.tri(neurons_distance[1:n,1:n])],
                   main= "Model 3");


# coefs/bin vs total for pair 11 & 6 (they are close together)

m <- get_rep_neurons();
n1 <- 2;
n2 <- 6;

mean_prod <- sapply(1:dim(m)[2], function(t) {
  cond_mean_prod_coef(m[n1,t,], m[n2,t,]);  
})

coef <- model1$param[n1,n2,];
num_joint_spikes <- sapply(1:dim(m)[2], function(t) {
  length(which(m[n1,t,] != 0 & m[n2,t,] != 0));  
});

cond_plot_meanprod(coef, mean_prod, model2.2$param[n1,n2], size= num_joint_spikes);
