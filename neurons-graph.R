# Created by Oleksandr Sorochynskyi
# On 14 06 18

library(plotrix)

source("open_neurons.R")
source("neurons.R")
source("marginals.R")
source("copula.R")

#########################################################################################
# Barplot of a neuron
x <- bin_neurons(2, 4, 0);
df <- data.frame(nspikes= factor(x))
ggplot(df, aes(x= nspikes)) +
  geom_bar(stat="count",
           width=0.7) +
  my.theme +
  xlab("Number of spikes") +
  ylab("Frequency");

barplot(table(x))

#########################################################################################
# Pretty copula plot 

copula_plot(data[1,], data[2,]);
copula_plot_diff(data[1,], data[2,], fit_copula(data[1,],data[2,],claytonCopula())$cop)

#######################################################################################

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

data <- bin_neurons(c(4,3),4,0);
t1 <- data[1,];
t2 <- data[2,];

plots <- mapply(function(cop,params) plot_copula(t1,t2,cop,params),
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

data_p <- copula_plot(t1,t2);


grid.arrange(data_p, data_p,
             plots[,1]$fit, plots[,1]$dif,
             plots[,2]$fit, plots[,2]$dif, 
             plots[,3]$fit, plots[,3]$dif,
             plots[,4]$fit, plots[,4]$dif,
             plots[,5]$fit, plots[,5]$dif,
             plots[,6]$fit, plots[,6]$dif,
             plots[,7]$fit, plots[,7]$dif,
             ncol=2, nrow=8);

my.theme <- theme(text= element_text(size=25),
               axis.text.x = element_text(angle=90, hjust=1),
               legend.position = "bottom",
               legend.key.width = unit(0.08, "npc"),
               legend.justification = "left"
               );

plot(data_p + my.theme)
dev.copy2pdf(file="plots/data_logdensity.pdf")

plot(plots[,1]$fit + my.theme)
dev.copy2pdf(file="plots/indep_logdensity.pdf")
plot(plots[,1]$dif + my.theme)
dev.copy2pdf(file="plots/indep_logratio.pdf")

plot(plots[,2]$fit + my.theme)
dev.copy2pdf(file="plots/clayton_logdensity.pdf")
plot(plots[,2]$dif + my.theme)
dev.copy2pdf(file="plots/clayton_logratio.pdf")

plot(plots[,3]$fit + my.theme)
dev.copy2pdf(file="plots/gumbel_logdensity.pdf")
plot(plots[,3]$dif + my.theme)
dev.copy2pdf(file="plots/gumbel_logratio.pdf")

plot(plots[,4]$fit + my.theme)
dev.copy2pdf(file="plots/frank_logdensity.pdf")
plot(plots[,4]$dif + my.theme)
dev.copy2pdf(file="plots/frank_logratio.pdf")

plot(plots[,5]$fit + my.theme)
dev.copy2pdf(file="plots/normal_logdensity.pdf")
plot(plots[,5]$dif + my.theme)
dev.copy2pdf(file="plots/normal_logratio.pdf")

plot(plots[,6]$fit + my.theme)
dev.copy2pdf(file="plots/rot_gumbel_logdensity.pdf")
plot(plots[,6]$dif + my.theme)
dev.copy2pdf(file="plots/rot_gumbel_logratio.pdf")

plot(plots[,7]$fit + my.theme)
dev.copy2pdf(file="plots/rot_clayton_logdensity.pdf")
plot(plots[,7]$dif + my.theme)
dev.copy2pdf(file="plots/rot_clayon_logratio.pdf")

#######################################################################################
# Plot densities using xyTable

data <- bin_neurons(c(1,2),20,0);
t1 <- data[1,];
t2 <- data[2,];
cdf1 <- ecdf(t1);
cdf2 <- ecdf(t2);

#plot data
tbl <- xyTable(cdf1(t1),cdf2(t2));
weight <- tbl$number*0.001
plot(tbl$x, tbl$y, cex=weight, pch=16,xlim=c(0,1),ylim=c(0,1));

#######################################################################################
# Graph using physical distance between neurons
library(igraph)

g <- graph_from_adjacency_matrix(neurons_distance,mode="undirected",weighted="w")
tkplot(g)
# Futurhman - Reigold layout should tak into account weignts, but it doesn't seem
# to give anything beter than random

# it's still just a blob
#########################################################################################
# Compare copulas for binned models

data <- bin_neurons(c(1,2),4);
x1 <- data[1,];
x2 <- data[2,];

logliks <- function(copula, cuts=1) {
  data <- bin_neurons(c(1,2),4);
  x1 <- data[1,];
  x2 <- data[2,];
  index <- cut(1:length(x1),cuts,labels=FALSE);
  fit <- fit_copula(x1[index==1],x2[index==1],copula,upper=10);
  logliks <- sapply(1:cuts, function(i) copula_loglik(x1[index == i],x2[index == i],setTheta(copula,fit$par)));
  logliks;
}
cops <- list(
  indepCopula(),
  frankCopula(),
  claytonCopula(),
  gumbelCopula()
)
lls <- lapply(cops, logliks)
lls <- lapply(cops, function(c) fit_copula(x1,x2,c,upper=10));
stripchart(lls);

#########################################################################################
# Compare copulas for binned models

# Let's start with simple correlation
cor_per_binsz <- function(binsz) {
  data <- bin_neurons(c(5,8), binsz);
  x1 <- data[1,]
  x2 <- data[2,]
  cor(x1,x2,method="sp");
}
cors <- sapply(1:10,cor_per_binsz);
plot(1:10,cors,ylim=c(0,1));

cop_per_binsz <- function(binsz) {
  data <- bin_neurons(c(1,2), binsz);  
  x1 <- data[1,]
  x2 <- data[2,]
  fit_copula(x1,x2,claytonCopula())$par;
}
cops <- sapply(1:10,cop_per_binsz)
plot(1:10, cops);
# Wow for claytonCopula this decreases
# For frank it does it first rises then drops off

loglik_per_binsz <- function(binsz) {
  data <- bin_neurons(c(8,4), binsz);  
  x1 <- data[1,]
  x2 <- data[2,]
  llc <- fit_copula(x1,x2,claytonCopula(),upper=6)$loglik/length(x1);
  llg <- fit_copula(x1,x2,gumbelCopula(),upper=6)$loglik/length(x1);
  llc-llg
}
y <- sapply(1:10,loglik_per_binsz);
plot(1:10,y);

#########################################################################################
# K-fold tests

x1 <- bin_neurons(1, 4);
x2 <- bin_neurons(2, 4);

res <- lapply(
  c(gumbelCopula(),
    claytonCopula(),
    frankCopula(),
    indepCopula()
    ),
  function(copula) { 
    res <- neuron.k_fold(x1,x2,copula,20,upper=10);
  });

xlim.min <- min(sapply(res, function(x) min(x$loglik)));
xlim.max <- max(sapply(res, function(x) max(x$loglik)));

lapply(res, function(x) hist(x$loglik, xlim=c(xlim.min,xlim.max)))
#########################################################################################
# Out of sample likelihood comparison between models

pairwise_loglik <- function(x,cop,params) {
  m <- matrix(NA, nrow=nrow(x),ncol=nrow(x))
  
  for (i in 2:nrow(x)) {
    for (j in 1:(i-1)) {
      copula <- setTheta(cop, params[i,j]);
      m[i,j] <- copula_loglik(x[i,],x[j,],copula);
    }
  }
  m
}

data <- bin_neurons(1:10, 4);

train_data <- data[,1:5000]
test_data <- data[,-(1:5000)]

clayton_fit <- pairwise_copula_fit(train_data,claytonCopula())
#fgm_fit <-pairwise_copula_fit(train_data,fgmCopula(),do.hessian=FALSE)
frank_fit <- pairwise_copula_fit(train_data,frankCopula())
gumbel_fit <- pairwise_copula_fit(train_data,gumbelCopula(),upper=6)
normal_fit <- pairwise_copula_fit(train_data,normalCopula())
#rot_clayton_fit <- pairwise_copula_fit(train_data,rotCopula(claytonCopula()),upper=10)
#t_fit <- pairwise_copula_fit(train_data,tCopula(dispstr = "un",df.fixed = TRUE),inital_val=0)

indep_loglik <- pairwise_loglik(test_data,indepCopula(),NULL)
clayton_loglik <- pairwise_loglik(test_data,claytonCopula(),clayton_fit$par)
#fgm_loglik <- pairwise_loglik(test_data,fgmCopula(), fgm_fit$par)
frank_loglik <- pairwise_loglik(test_data,frankCopula(),frank_fit$par)
gumbel_loglik <- pairwise_loglik(test_data,gumbelCopula(),gumbel_fit$par)
normal_loglik <- pairwise_loglik(test_data,normalCopula(),normal_fit$par)
#rot_clayton_loglik <- pairwise_loglik(test_data,rotCopula(claytonCopula()), normal_fit$par)
#t_loglik <- pairwise_loglik(test_data,tCopula(df.fixed = TRUE),t_fit$par)


par(mfrow=c(4,1))
hist(indep_loglik,main="indepe",xlim=c(-3000,-0),col=sample(colors(),1),ylim=c(0,15))
hist(clayton_loglik,main="clayton",xlim=c(-3000,-0),col=sample(colors(),1),ylim=c(0,15))
#hist(fgm_loglik,main="fgm",xlim=c(-3000,-0),col=sample(colors(),1),ylim=c(0,15))
hist(frank_loglik,main="frank",xlim=c(-3000,-0),col=sample(colors(),1),ylim=c(0,15))
hist(gumbel_loglik,main="gumbel",xlim=c(-3000,-0),col=sample(colors(),1),ylim=c(0,15))
hist(normal_loglik,main="normal",xlim=c(-3000,-0),col=sample(colors(),1),ylim=c(0,15))
#hist(rot_clayton_loglik,main="rot_clayton",xlim=c(-3000,-0),col=sample(colors(),1),ylim=c(0,15))
#hist(t_loglik,main="t",xlim=c(-30000,-10000),col=sample(colors(),1),ylim=c(0,15))


cor(as.vector(clayton_loglik),as.vector(fgm_loglik),use="pair");
cor(as.vector(clayton_loglik),as.vector(normal_loglik),use="pair");

sapply(list(clayton_loglik,frank_loglik,gumbel_loglik,normal_loglik),function(x) x[3,1])

copula_loglik_test(neurons[3,],
                   neurons[1,],
                   frankCopula(frank_fit$par[3,1]),
                   normalCopula(normal_fit$par[3,1]))
#########################################################################################
# compare all copula fits
read_fit <- function(copula_name) {
  par <- unname(as.matrix(read.table(paste(copula_name,"10fit-par.txt",sep=""))));
  hess <- unname(as.matrix(read.table(paste(copula_name,"10fit-hessian.txt",sep=""))));
  loglik_max <- unname(as.matrix(read.table(paste(copula_name,"10fit-loglik_max.txt",sep=""))));
  
  list(par=par,hessian=hess,loglik_max=loglik_max) 
}

clayton_fit <- read_fit("clayton")
fgm_fit <-read_fit("fgm");
frank_fit <- read_fit("frank");
gumbel_fit <- read_fit("gumbel");
indep_loglik <- unname(as.matrix(read.table("independant10fit-loglik_max.txt")));
normal_fit <- read_fit("normal");
rot_clayton_fit <- read_fit("rot_clayton");
t_fit <- read_fit("t");

par(mfrow=c(4,1))
hist(clayton_fit$loglik_max,main="clayton",xlim=c(-70000,-10000),col=sample(colors(),1),ylim=c(0,15))
hist(indep_loglik,main="indepe",xlim=c(-70000,-10000),col=sample(colors(),1),ylim=c(0,15))
hist(fgm_fit$loglik_max,main="fgm",xlim=c(-70000,-10000),col=sample(colors(),1),ylim=c(0,15))
hist(frank_fit$loglik_max,main="frank",xlim=c(-70000,-10000),col=sample(colors(),1),ylim=c(0,15))
hist(gumbel_fit$loglik_max,main="gumbel",xlim=c(-70000,-10000),col=sample(colors(),1),ylim=c(0,15))
hist(normal_fit$loglik_max,main="normal",xlim=c(-70000,-10000),col=sample(colors(),1),ylim=c(0,15))
hist(rot_clayton_fit$loglik_max,main="rot_clayton",xlim=c(-70000,-10000),col=sample(colors(),1),ylim=c(0,15))
hist(t_fit$loglik_max,main="t",xlim=c(-70000,-10000),col=sample(colors(),1),ylim=c(0,15))

# fgm is clearly misspecified, as for all estimation it just hits the boundairy, fgm
# copulas can only capture a mild dependance so this sin't surprising.

# let's see if we get simular loglikelioods for different copulas
cor(as.vector(clayton_fit$loglik_max),as.vector(fgm_fit$loglik_max),use="pair")
cor(as.vector(clayton_fit$loglik_max),as.vector(gumbel_fit$loglik_max),use="pair")
cor(as.vector(normal_fit$loglik_max),as.vector(gumbel_fit$loglik_max),use="pair")


# looking at this it seems that that the loglikelihoods for a given copula are very!
# correlated
# Maybe the choice of copula doesn't even matter that much. 

#######################################################################################
# Pairwise independence tests

data <- neurons[1:10,]

m <- matrix(NA,ncol=nrow(data),nrow=nrow(data))
for (i in 2:nrow(data)) {
  for (j in 1:(i-1)) {
    m[i,j] <- copula_kramer_indep_test(data[i,],data[j,])$p.value 
  }
}
copula_kramer_indep_test(data[1,],data[2,])
#######################################################################################
# show that copulas actually ned differantiating

data <- bin_neurons(c(9,10), 4);

t1 <- data[1,];
t2 <- data[2,];

sample_index <- sample(1:length(t1), floor(length(t1)/4))
t1s <- t1[sample_index];
t2s <- t2[sample_index];

indep_fit <- fit_copula(t1s,t2s, indepCopula())
rot_clayton_fit <- fit_copula(t1s,t2s, rotCopula(claytonCopula(),rep(TRUE,2)),upper=6);
clayton_fit <- fit_copula(t1s,t2s, claytonCopula())
gumbel_fit <- fit_copula(t1s,t2s, gumbelCopula(), upper=6)
frank_fit <- fit_copula(t1s,t2s, frankCopula())
normal_fit <- fit_copula(t1s,t2s, normalCopula())
#t_fit <- fit_copula(t1s,t2s, tCopula())
#galambos_fit <- fit_copula(t1s,t2s, galambosCopula(),upper=2)
#huslerReiss_fit <- fit_copula(t1s,t2s, huslerReissCopula(),upper=1)
tawn_fit <- fit_copula(t1s,t2s, tawnCopula())
tev_fit <- fit_copula(t1s,t2s, tevCopula())
#amh_fit <- fit_copula(t1s,t2s, amhCopula(),lower=0,upper=1)
#joe_fit <- fit_copula(t1s,t2s, joeCopula(),upper=6)
fgm_fit <- fit_copula(t1s,t2s, fgmCopula())

t1t <- t1[-sample_index]
t2t <- t2[-sample_index]

# test sample
par(mfrow=c(5,1))
mapply(function(copula) {
         data <- replicate(200, {
            sample_index <- sample(1:length(t1t),floor(length(t1t)/4))
            t1t <- t1[-sample_index]
            t2t <- t2[-sample_index]
            copula_loglik(t1t, t2t, copula)
         });
         hist(
            data,
            breaks=10,
            xlab="Loglik",
            main=class(copula),
            col=sample(colors(),1)
            ,xlim=c(-23500, -22000)
          );
       },
       list(#indepCopula(),
            claytonCopula(clayton_fit$par),
            gumbelCopula(gumbel_fit$par),
            frankCopula(frank_fit$par),
            rotCopula(claytonCopula(rot_clayton_fit$par)),
            normalCopula(normal_fit$par)#,
            #tawnCopula(tawn_fit$par),
            #tevCopula(tev_fit$par),
            #fgmCopula(fgm_fit$par)
            )
);

# See how theese graphs compare to our tests
#######################################################################################
# How different copulas compare to one another

# Spearman's rho cor matrix
cor_m <- spearman_cor_matrix(neurons[1:10,]);
# Copula parameters
gumbel_m <- as.matrix(read.table("gumbel10fit-par.txt"));
rot_clayton_m <- as.matrix(read.table("rot_clayton10fit-par.txt"));
frank_m <- as.matrix(read.table("frank10fit-par.txt"))
normal_m <- as.matrix(read.table("normal10fit-par.txt"))
t_m <- as.matrix(read.table("t10fit-par.txt"))

color2D.matplot(cor_m);
color2D.matplot(gumbel_m);
color2D.matplot(rot_clayton_m);
color2D.matplot(frank_m);
color2D.matplot(normal_m);
color2D.matplot(t_m);

#######################################################################################
# How Spearmans rho compares to gumbel theta

# Spearman's rho cor matrix
cor_m <- spearman_cor_matrix(neurons[1:10,]);
# Estimated parametsrs using clayton copula
cop_m <- as.matrix(read.table("gumbel10fit-par.txt"));

# plot cop vs cor coefs, 
plot(cop_m, cor_m);

# Let's compare that to what rho should be equal to given gumbel theta
x <- seq(1,6,by=0.1);
y <- lapply(x, function(t) rho(gumbelCopula(t,2)))

# what mle gives us
plot(cop_m, cor_m, xlim=c(1,6));
# what we're supposed to get if we listen to Spearmans rho
points(x,y,col="red")

#######################################################################################

# pretty color matrix cov

# Spearman's rho cor matrix
cor_m <- spearman_cor_matrix(neurons[1:10,]);
# Estimated parametsrs using clayton copula
cop_m <- as.matrix(read.table("gumbel10fit-par.txt"));

color2D.matplot(cor_m);
color2D.matplot(cop_m);

#######################################################################################
# joint hist
t1 <- neurons [1,]
t2 <- neurons [2,]
joint_freq <- matrix(0, nrow= max(c(t1,t2)), ncol= max(c(t1,t2)))
for (i in seq_along(t1)) {
  joint_freq[t1[i],t2[i]] = joint_freq[t1[i],t2[i]] + 1
}
color2D.matplot(joint_freq,cs1=c(0,256))

#######################################################################################
# Compare how marginals fit

t <- neurons[1,]

# Poisson fit
fit_poiss <- margin_poiss(t);

barplot(table(t), col=rgb(0.1,0.1,0.1,0.5))
barplot(table(fit_poiss$ran(length(t))), add= TRUE, col=rgb(0.8,0.8,0.8,0.5))

# Negative binomeal
fit_nbinom <- margin_nbinom(t);

barplot(table(t), col=rgb(0.1,0.1,0.1,0.5))
barplot(table(rnbinom(length(t), fit_nbinom$param[1],fit_nbinom$param[2])), add= TRUE, col=rgb(0.8,0.8,0.8,0.5))

#######################################################################################
# Try to apply community search algorithms to our data
# Result: it's on big blob, or one big interconnected comunity

# look for community ;)
library(igraph)

g <- graph_from_adjacency_matrix(cop_m,mode="undirected",weighted="w")
community <- cluster_walktrap(g)
plot(community, g)
# Nice, there's one big community, this is useless

#######################################################################################
# Again community search but this time with thresholds 

library(igraph)

t <- neurons[1:10,]
cop_m <- as.matrix(read.table("gumbel10fit-par.txt"));
cop_m <- spearman_cor_matrix(t)
cop_m[upper.tri(cop_m)] <- t(cop_m)[upper.tri(cop_m)]
diag(cop_m)<-NA

threshold <- 1.7;
cop_m[cop_m < threshold] <- 0;
# OR
sigmoid <- function(x) (1)*exp(-(1)*exp(-(20)*(x-.5)))
plot(sigmoid,xlim=c(0,1))
cop_m <- sigmoid(cop_m)


g <- graph_from_adjacency_matrix(cop_m,mode="undirected",weighted="w",diag=FALSE)
community <- cluster_edge_betweenness(g)
plot(community, g)

