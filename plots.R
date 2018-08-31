# Created by Oleksandr Sorochynskyi
# On 10 8 18

# This file will contain all the plots for my Journal Report

source('~/IdV/R/conditional.R')

library(ggplot2)

my.theme <- theme(text= element_text(size=25));
write2pdf <- function(filename) dev.copy2pdf(file=paste("plots/",filename,sep=""), height=5, width=5);
tri <- function(m) as.vector(m[lower.tri(m)]);

pair_near <- c(23,12)
pair_far <- c(23,3)
pair_median <- c(23,2)

bin.size = 8/60; # 16*8=128

# params
load("summary/params.RData")
# covar 
load("summary/covar.RData")
# distances
load("summary/dist.RData")

#######################################################################################################
# Scatterplots for copulas
# Gaussin clayton fmg frank

# Fix the dependance level at tau=0.5
tau=0.6

fit_gauss <- normalCopula(iTau(normalCopula(),tau));
fit_clayton <- claytonCopula(iTau(claytonCopula(),tau));
fit_fmg <- fgmCopula(iTau(fgmCopula(),tau));
fit_frank <- frankCopula(iTau(frankCopula(),tau));

rand<- rCopula(1000, indepCopula());
ggplot(data=data.frame(rand=rand), aes(x=rand[,1], y=rand[,2])) +
  geom_point() +
  my.theme+
  theme( axis.title= element_blank() )
write2pdf("copula_scatter/indep.pdf");

rand<- rCopula(1000, fit_gauss);
ggplot(data=data.frame(rand=rand), aes(x=rand[,1], y=rand[,2])) +
  geom_point() +
  my.theme+
  theme( axis.title= element_blank() )
write2pdf("copula_scatter/gauss.pdf");

rand<- rCopula(1000, fit_clayton);
ggplot(data=data.frame(rand=rand), aes(x=rand[,1], y=rand[,2])) +
  geom_point() +
  my.theme+
  theme( axis.title= element_blank() )
write2pdf("copula_scatter/clayton.pdf");

rand<- rCopula(1000, fit_fmg);
ggplot(data=data.frame(rand=rand), aes(x=rand[,1], y=rand[,2])) +
  geom_point() +
  my.theme+
  theme( axis.title= element_blank() )
write2pdf("copula_scatter/fmg.pdf");

rand<- rCopula(1000, fit_frank);
ggplot(data=data.frame(rand=rand), aes(x=rand[,1], y=rand[,2])) +
  geom_point() +
  my.theme +
  theme( axis.title= element_blank() )
write2pdf("copula_scatter/frank.pdf");

################################################################################
# Barplot of the neuron distribution


m <- bin_fullfield(pair_near[1]);
ggplot(data.frame(spike.count= as.vector(m)), aes(x=spike.count)) +
  geom_bar() +
  labs(x="Number of spikes") +
  my.theme +
  theme(axis.title.y = element_blank());
write2pdf("neuron/barplot.pdf")

################################################################################
# Model3 for OFF neurons 

# Distance

ggplot(mapping=aes(x=tri(distance[1:25,1:25]), y=tri(params$model3$par$fullfield[1:25,1:25]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model3/off_dist_vs_fullfield.pdf")

ggplot(mapping= aes(x=tri(distance[1:25,1:25]), y=tri(params$model3$par$bar[1:25,1:25]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model3/off_dist_vs_bar.pdf")

ggplot(mapping= aes(x=tri(distance[1:25,1:25]), y=tri(params$model3$par$check[1:25,1:25]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model3/off_dist_vs_check.pdf")


# Covariance prediction

# Stimul correlation 

ggplot(mapping= aes(x=tri(covar$fullfield_data$check_param$model3$real$stimul[1:25,1:25]), y=tri(covar$fullfield_data$check_param$model3$simu$stimul[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  ylim(c(-.015,0.01)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/off_cov_stimul_fullfield_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$check_data$check_param$model3$real$stimul[1:25,1:25]), y=tri(covar$check_data$check_param$model3$simu$stimul[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  ylim(c(-.015,0.01)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/off_cov_stimul_check_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$bar_data$check_param$model3$real$stimul[1:25,1:25]), y=tri(covar$bar_data$check_param$model3$simu$stimul[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  ylim(c(-0.01,0.05)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/off_cov_stimul_bar_data_check_param.pdf")

# Noise correlation 

ggplot(mapping= aes(x=tri(covar$fullfield_data$check_param$model3$real$noise[1:25,1:25]), y=tri(covar$fullfield_data$check_param$model3$simu$noise[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/off_cov_noise_fullfield_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$check_data$check_param$model3$real$noise[1:25,1:25]), y=tri(covar$check_data$check_param$model3$simu$noise[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/off_cov_noise_check_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$bar_data$check_param$model3$real$noise[1:25,1:25]), y=tri(covar$bar_data$check_param$model3$simu$noise[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-0.01,0.05)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/off_cov_noise_bar_data_check_param.pdf")


# Parameter vs parameter

ggplot(mapping= aes(x=tri(params$model3$par$bar[1:25,1:25]), y=tri(params$model3$par$fullfield[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Fullfield stimulus", y="Bar movie stimulus") +
  my.theme;
write2pdf("model3/off_param_fullfield_vs_bar.pdf")

ggplot(mapping= aes(x=tri(params$model3$par$fullfield[1:25,1:25]), y=tri(params$model3$par$check[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Fullfield stimulus", y="Checkerboard stimulus") +
  my.theme;
write2pdf("model3/off_param_fullfield_vs_check.pdf")

ggplot(mapping= aes(x=tri(params$model3$par$bar[1:25,1:25]), y=tri(params$model3$par$check[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Bar movie stimulus", y="Checkerboard stimulus") +
  my.theme;
write2pdf("model3/off_param_bar_vs_check.pdf")

################################################################################
# Model3 for ON neurons 

# Distance

ggplot(mapping=aes(x=tri(distance[-(1:25),-(1:25)]), y=tri(params$model3$par$fullfield[-(1:25),-(1:25)]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model3/on_dist_vs_fullfield.pdf")

ggplot(mapping= aes(x=tri(distance[-(1:25),-(1:25)]), y=tri(params$model3$par$bar[-(1:25),-(1:25)]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model3/on_dist_vs_bar.pdf")

ggplot(mapping= aes(x=tri(distance[-(1:25),-(1:25)]), y=tri(params$model3$par$check[-(1:25),-(1:25)]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model3/on_dist_vs_check.pdf")


# Covariance prediction

# Stimul correlation 

ggplot(mapping= aes(x=tri(covar$fullfield_data$check_param$model3$real$stimul[-(1:25),-(1:25)]), y=tri(covar$fullfield_data$check_param$model3$simu$stimul[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  ylim(c(-.015,0.01)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/on_cov_stimul_fullfield_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$check_data$check_param$model3$real$stimul[-(1:25),-(1:25)]), y=tri(covar$check_data$check_param$model3$simu$stimul[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  ylim(c(-.015,0.01)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/on_cov_stimul_check_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$bar_data$check_param$model3$real$stimul[-(1:25),-(1:25)]), y=tri(covar$bar_data$check_param$model3$simu$stimul[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  ylim(c(-0.01,0.05)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/on_cov_stimul_bar_data_check_param.pdf")

# Noise correlation 

ggplot(mapping= aes(x=tri(covar$fullfield_data$check_param$model3$real$noise[-(1:25),-(1:25)]), y=tri(covar$fullfield_data$check_param$model3$simu$noise[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/on_cov_noise_fullfield_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$check_data$check_param$model3$real$noise[-(1:25),-(1:25)]), y=tri(covar$check_data$check_param$model3$simu$noise[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/on_cov_noise_check_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$bar_data$check_param$model3$real$noise[-(1:25),-(1:25)]), y=tri(covar$bar_data$check_param$model3$simu$noise[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-0.01,0.05)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/on_cov_noise_bar_data_check_param.pdf")


# Parameter vs parameter

ggplot(mapping= aes(x=tri(params$model3$par$bar[-(1:25),-(1:25)]), y=tri(params$model3$par$fullfield[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Fullfield stimulus", y="Bar movie stimulus") +
  my.theme;
write2pdf("model3/on_param_fullfield_vs_bar.pdf")

ggplot(mapping= aes(x=tri(params$model3$par$fullfield[-(1:25),-(1:25)]), y=tri(params$model3$par$check[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Fullfield stimulus", y="Checkerboard stimulus") +
  my.theme;
write2pdf("model3/on_param_fullfield_vs_check.pdf")

ggplot(mapping= aes(x=tri(params$model3$par$bar[-(1:25),-(1:25)]), y=tri(params$model3$par$check[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Bar movie stimulus", y="Checkerboard stimulus") +
  my.theme;
write2pdf("model3/on_param_bar_vs_check.pdf")

################################################################################
# Model3 for ON/OFF neuron pairs BAR and CHeck

# Distance

ggplot(mapping=aes(x=tri(distance[-(1:25),1:25]), y=tri(params$model3$par$fullfield[-(1:25),1:25]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model3/on_off_dist_vs_fullfield.pdf")

ggplot(mapping= aes(x=tri(distance[-(1:25),1:25]), y=tri(params$model3$par$bar[-(1:25),1:25]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model3/on_off_dist_vs_bar.pdf")

ggplot(mapping= aes(x=tri(distance[-(1:25),1:25]), y=tri(params$model3$par$check[-(1:25),1:25]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model3/on_off_dist_vs_check.pdf")


# Covariance prediction

# Stimul correlation 

ggplot(mapping= aes(x=tri(covar$fullfield_data$check_param$model3$real$stimul[-(1:25),1:25]), y=tri(covar$fullfield_data$check_param$model3$simu$stimul[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  ylim(c(-.015,0.01)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/on_off_cov_stimul_fullfield_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$check_data$check_param$model3$real$stimul[-(1:25),1:25]), y=tri(covar$check_data$check_param$model3$simu$stimul[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  ylim(c(-.015,0.01)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/on_off_cov_stimul_check_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$bar_data$check_param$model3$real$stimul[-(1:25),1:25]), y=tri(covar$bar_data$check_param$model3$simu$stimul[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  ylim(c(-0.01,0.05)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/on_off_cov_stimul_bar_data_check_param.pdf")

# Noise correlation 

ggplot(mapping= aes(x=tri(covar$fullfield_data$check_param$model3$real$noise[-(1:25),1:25]), y=tri(covar$fullfield_data$check_param$model3$simu$noise[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/on_off_cov_noise_fullfield_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$check_data$check_param$model3$real$noise[-(1:25),1:25]), y=tri(covar$check_data$check_param$model3$simu$noise[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/on_off_cov_noise_check_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$bar_data$check_param$model3$real$noise[-(1:25),1:25]), y=tri(covar$bar_data$check_param$model3$simu$noise[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-0.01,0.05)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model3/on_off_cov_noise_bar_data_check_param.pdf")


# Parameter vs parameter

ggplot(mapping= aes(x=tri(params$model3$par$bar[-(1:25),1:25]), y=tri(params$model3$par$fullfield[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Fullfield stimulus", y="Bar movie stimulus") +
  my.theme;
write2pdf("model3/on_off_param_fullfield_vs_bar.pdf")

ggplot(mapping= aes(x=tri(params$model3$par$fullfield[-(1:25),1:25]), y=tri(params$model3$par$check[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Fullfield stimulus", y="Checkerboard stimulus") +
  my.theme;
write2pdf("model3/on_off_param_fullfield_vs_check.pdf")

ggplot(mapping= aes(x=tri(params$model3$par$bar[-(1:25),1:25]), y=tri(params$model3$par$check[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Bar movie stimulus", y="Checkerboard stimulus") +
  my.theme;
write2pdf("model3/on_off_param_bar_vs_check.pdf")

################################################################################
# Model2.2 for OFF neurons 

# Distance

ggplot(mapping=aes(x=tri(distance[1:25,1:25]), y=tri(params$model2.2$par$fullfield[1:25,1:25]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model2.2/off_dist_vs_fullfield.pdf")

ggplot(mapping= aes(x=tri(distance[1:25,1:25]), y=tri(params$model2.2$par$bar[1:25,1:25]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model2.2/off_dist_vs_bar.pdf")

ggplot(mapping= aes(x=tri(distance[1:25,1:25]), y=tri(params$model2.2$par$check[1:25,1:25]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model2.2/off_dist_vs_check.pdf")


# Covariance prediction

# Stimul correlation 

ggplot(mapping= aes(x=tri(covar$fullfield_data$check_param$model2.2$real$stimul[1:25,1:25]), y=tri(covar$fullfield_data$check_param$model2.2$simu$stimul[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/off_cov_stimul_fullfield_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$check_data$check_param$model2.2$real$stimul[1:25,1:25]), y=tri(covar$check_data$check_param$model2.2$simu$stimul[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/off_cov_stimul_check_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$bar_data$check_param$model2.2$real$stimul[1:25,1:25]), y=tri(covar$bar_data$check_param$model2.2$simu$stimul[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-0.01,0.05)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/off_cov_stimul_bar_data_check_param.pdf")

# Noise correlation 

ggplot(mapping= aes(x=tri(covar$fullfield_data$check_param$model2.2$real$noise[1:25,1:25]), y=tri(covar$fullfield_data$check_param$model2.2$simu$noise[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/off_cov_noise_fullfield_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$check_data$check_param$model2.2$real$noise[1:25,1:25]), y=tri(covar$check_data$check_param$model2.2$simu$noise[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/off_cov_noise_check_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$bar_data$check_param$model2.2$real$noise[1:25,1:25]), y=tri(covar$bar_data$check_param$model2.2$simu$noise[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-0.01,0.05)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/off_cov_noise_bar_data_check_param.pdf")


# Parameter vs parameter

ggplot(mapping= aes(x=tri(params$model2.2$par$bar[1:25,1:25]), y=tri(params$model2.2$par$fullfield[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Fullfield stimulus", y="Bar movie stimulus") +
  my.theme;
write2pdf("model2.2/off_param_fullfield_vs_bar.pdf")

ggplot(mapping= aes(x=tri(params$model2.2$par$fullfield[1:25,1:25]), y=tri(params$model2.2$par$check[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Fullfield stimulus", y="Checkerboard stimulus") +
  my.theme;
write2pdf("model2.2/off_param_fullfield_vs_check.pdf")

ggplot(mapping= aes(x=tri(params$model2.2$par$bar[1:25,1:25]), y=tri(params$model2.2$par$check[1:25,1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Bar movie stimulus", y="Checkerboard stimulus") +
  my.theme;
write2pdf("model2.2/off_param_bar_vs_check.pdf")

################################################################################
# Model2.2 for ON neurons 

# Distance

ggplot(mapping=aes(x=tri(distance[-(1:25),-(1:25)]), y=tri(params$model2.2$par$fullfield[-(1:25),-(1:25)]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model2.2/on_dist_vs_fullfield.pdf")

ggplot(mapping= aes(x=tri(distance[-(1:25),-(1:25)]), y=tri(params$model2.2$par$bar[-(1:25),-(1:25)]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model2.2/on_dist_vs_bar.pdf")

ggplot(mapping= aes(x=tri(distance[-(1:25),-(1:25)]), y=tri(params$model2.2$par$check[-(1:25),-(1:25)]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model2.2/on_dist_vs_check.pdf")


# Covariance prediction

# Stimul correlation 

ggplot(mapping= aes(x=tri(covar$fullfield_data$check_param$model2.2$real$stimul[-(1:25),-(1:25)]), y=tri(covar$fullfield_data$check_param$model2.2$simu$stimul[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/on_cov_stimul_fullfield_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$check_data$check_param$model2.2$real$stimul[-(1:25),-(1:25)]), y=tri(covar$check_data$check_param$model2.2$simu$stimul[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/on_cov_stimul_check_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$bar_data$check_param$model2.2$real$stimul[-(1:25),-(1:25)]), y=tri(covar$bar_data$check_param$model2.2$simu$stimul[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-0.01,0.05)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/on_cov_stimul_bar_data_check_param.pdf")

# Noise correlation 

ggplot(mapping= aes(x=tri(covar$fullfield_data$check_param$model2.2$real$noise[-(1:25),-(1:25)]), y=tri(covar$fullfield_data$check_param$model2.2$simu$noise[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/on_cov_noise_fullfield_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$check_data$check_param$model2.2$real$noise[-(1:25),-(1:25)]), y=tri(covar$check_data$check_param$model2.2$simu$noise[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/on_cov_noise_check_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$bar_data$check_param$model2.2$real$noise[-(1:25),-(1:25)]), y=tri(covar$bar_data$check_param$model2.2$simu$noise[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-0.01,0.05)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/on_cov_noise_bar_data_check_param.pdf")


# Parameter vs parameter

ggplot(mapping= aes(x=tri(params$model2.2$par$bar[-(1:25),-(1:25)]), y=tri(params$model2.2$par$fullfield[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Fullfield stimulus", y="Bar movie stimulus") +
  my.theme;
write2pdf("model2.2/on_param_fullfield_vs_bar.pdf")

ggplot(mapping= aes(x=tri(params$model2.2$par$fullfield[-(1:25),-(1:25)]), y=tri(params$model2.2$par$check[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Fullfield stimulus", y="Checkerboard stimulus") +
  my.theme;
write2pdf("model2.2/on_param_fullfield_vs_check.pdf")

ggplot(mapping= aes(x=tri(params$model2.2$par$bar[-(1:25),-(1:25)]), y=tri(params$model2.2$par$check[-(1:25),-(1:25)]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Bar movie stimulus", y="Checkerboard stimulus") +
  my.theme;
write2pdf("model2.2/on_param_bar_vs_check.pdf")

################################################################################
# Model2.2 for ON/OFF neuron pairs BAR and CHeck

# Distance

ggplot(mapping=aes(x=tri(distance[-(1:25),1:25]), y=tri(params$model2.2$par$fullfield[-(1:25),1:25]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model2.2/on_off_dist_vs_fullfield.pdf")

ggplot(mapping= aes(x=tri(distance[-(1:25),1:25]), y=tri(params$model2.2$par$bar[-(1:25),1:25]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model2.2/on_off_dist_vs_bar.pdf")

ggplot(mapping= aes(x=tri(distance[-(1:25),1:25]), y=tri(params$model2.2$par$check[-(1:25),1:25]))) +
  geom_point() +
  labs(x="Distance (in micrometers)", y="Model parameter") +
  my.theme;
write2pdf("model2.2/on_off_dist_vs_check.pdf")


# Covariance prediction

# Stimul correlation 

ggplot(mapping= aes(x=tri(covar$fullfield_data$check_param$model2.2$real$stimul[-(1:25),1:25]), y=tri(covar$fullfield_data$check_param$model2.2$simu$stimul[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/on_off_cov_stimul_fullfield_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$check_data$check_param$model2.2$real$stimul[-(1:25),1:25]), y=tri(covar$check_data$check_param$model2.2$simu$stimul[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/on_off_cov_stimul_check_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$bar_data$check_param$model2.2$real$stimul[-(1:25),1:25]), y=tri(covar$bar_data$check_param$model2.2$simu$stimul[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-0.01,0.05)) +
  labs(x="Observed stimulation cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/on_off_cov_stimul_bar_data_check_param.pdf")

# Noise correlation 

ggplot(mapping= aes(x=tri(covar$fullfield_data$check_param$model2.2$real$noise[-(1:25),1:25]), y=tri(covar$fullfield_data$check_param$model2.2$simu$noise[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/on_off_cov_noise_fullfield_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$check_data$check_param$model2.2$real$noise[-(1:25),1:25]), y=tri(covar$check_data$check_param$model2.2$simu$noise[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-.015,0.01)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/on_off_cov_noise_check_data_check_param.pdf")

ggplot(mapping= aes(x=tri(covar$bar_data$check_param$model2.2$real$noise[-(1:25),1:25]), y=tri(covar$bar_data$check_param$model2.2$simu$noise[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  #ylim(c(-0.01,0.05)) +
  labs(x="Observed noise cov", y="Cov predicted by the model") +
  my.theme;
write2pdf("model2.2/on_off_cov_noise_bar_data_check_param.pdf")


# Parameter vs parameter

ggplot(mapping= aes(x=tri(params$model2.2$par$bar[-(1:25),1:25]), y=tri(params$model2.2$par$fullfield[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Fullfield stimulus", y="Bar movie stimulus") +
  my.theme;
write2pdf("model2.2/on_off_param_fullfield_vs_bar.pdf")

ggplot(mapping= aes(x=tri(params$model2.2$par$fullfield[-(1:25),1:25]), y=tri(params$model2.2$par$check[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Fullfield stimulus", y="Checkerboard stimulus") +
  my.theme;
write2pdf("model2.2/on_off_param_fullfield_vs_check.pdf")

ggplot(mapping= aes(x=tri(params$model2.2$par$bar[-(1:25),1:25]), y=tri(params$model2.2$par$check[-(1:25),1:25]))) +
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  labs(x="Bar movie stimulus", y="Checkerboard stimulus") +
  my.theme;
write2pdf("model2.2/on_off_param_bar_vs_check.pdf")

################################################################################
# Model 3 density plots

library(gridExtra);
plot_copula <- function(x1, x2, copula, params) {
  cdf1 <- ecdf(x1);
  cdf2 <- ecdf(x2);
  options = params;
  options[["x1"]] = x1;
  options[["x2"]] = x2;
  options[["copula"]] = copula;
  fit <- do.call(fit_copula,options);
  
  
  p1 <- copula_plot_cop_density(fit$cop, x1, x2)# + ggtitle(class(copula)); 
  p2 <- copula_plot_diff(x1,x2, fit$cop)# + ggtitle("log(data+1/rand+1)")
 
  list(fit= p1, dif= p2, cop_fit = fit); 
}

# fullfield stimulusm

m <- bin_fullfield(pair_near,bin.size);

plots <- mapply(function(cop,params) plot_copula(as.vector(m[1,,]),as.vector(m[2,,]),cop,params),
       list(
         indepCopula(),
         claytonCopula(),
         gumbelCopula(),
         frankCopula(),
         normalCopula(),
         rotCopula(gumbelCopula()),
         rotCopula(claytonCopula())
       ),
       list(
         list(),
         list(),
         list(upper=6),
         list(),
         list(),
         list(upper=6),
         list(upper=6)
       )
);

data_p <- copula_plot(as.vector(m[1,,]),as.vector(m[2,,]));

settings <- my.theme +
               theme(axis.text.x = element_text(angle=90, hjust=1),
               legend.position = "bottom",
               legend.key.width = unit(0.08, "npc"),
               legend.justification = "left"
               );

plot(data_p + settings)
write2pdf("copula_compare/data_logdensity_fullfield_near.pdf")

plot(plots[,1]$fit + settings)
write2pdf("copula_compare/indep_logdensity_fullfield_near.pdf")
plot(plots[,1]$dif + settings)
write2pdf("copula_compare/indep_logratio_fullfield_near.pdf")

plot(plots[,2]$fit + settings)
write2pdf("copula_compare/clayton_logdensity_fullfield_near.pdf")
plot(plots[,2]$dif + settings)
write2pdf("copula_compare/clayton_logratio_fullfield_near.pdf")

plot(plots[,3]$fit + settings)
write2pdf("copula_compare/gumbel_logdensity_fullfield_near.pdf")
plot(plots[,3]$dif + settings)
write2pdf("copula_compare/gumbel_logratio_fullfield_near.pdf")

plot(plots[,4]$fit + settings)
write2pdf("copula_compare/frank_logdensity_fullfield_near.pdf")
plot(plots[,4]$dif + settings)
write2pdf("copula_compare/frank_logratio_fullfield_near.pdf")

plot(plots[,5]$fit + settings)
write2pdf("copula_compare/normal_logdensity_fullfield_near.pdf")
plot(plots[,5]$dif + settings)
write2pdf("copula_compare/normal_logratio_fullfield_near.pdf")

plot(plots[,6]$fit + settings)
write2pdf("copula_compare/rot_gumbel_logdensity_fullfield_near.pdf")
plot(plots[,6]$dif + settings)
write2pdf("copula_compare/rot_gumbel_logratio_fullfield_near.pdf")

plot(plots[,7]$fit + settings)
write2pdf("copula_compare/rot_clayton_logdensity_fullfield_near.pdf")
plot(plots[,7]$dif + settings)
write2pdf("copula_compare/rot_clayon_logratio_fullfield_near.pdf")

