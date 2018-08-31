# Created by Oleksandr Sorochynskyi
# On 14 06 18

source("open_neurons.R")
source("neurons.R")


#########################################################################################
# Fit model (with marginals) to one pair, and then try it on another

pair1 <- bin_neurons(c(1,2),4,0);
pair2 <- bin_neurons(c(3,4),4,0);
x1 <- pair1[1,]
x2 <- pair1[2,]
y1 <- pair2[1,]
y2 <- pair2[2,]

# rerun same thing but on y1, y2
# x1 <- y1
# x2 <- y2

mfit1 <- margin_nbinom(x1);
mfit2 <- margin_nbinom(x2);
fit <- fit_copula(x1,x2,frankCopula(),mfit1$cdf,mfit2$cdf,upper=100)
fit <- fit_copula(x1,x2,gumbelCopula(),mfit1$cdf,mfit2$cdf,upper=6)
fit <- fit_copula(x1,x2,claytonCopula(),mfit1$cdf,mfit2$cdf)
fit <- fit_copula(x1,x2,normalCopula(),mfit1$cdf,mfit2$cdf)

copula_loglik(x1=y1,x2=y2,cop=fit$cop,cdf1=mfit1$cdf,cdf2=mfit2$cdf)


#########################################################################################
# distance

c <- cor(t(bin_neurons(1:25,1)),t(bin_neurons(1:25,1)),method="p")
c_vec <- c[upper.tri(c)];
nd_vec <- neurons_distance[upper.tri(neurons_distance)];
cor(c_vec,nd_vec);
plot(c_vec, nd_vec);

index <- which(nd_vec <= 400);
cor(c_vec[index],nd_vec[index]);
# Correlation seems to be most present when distance goes to zero

inv.distance.lm <- lm(c_vec ~ 0 + I(1/nd_vec));
idist.fitted.val <- fitted.values(inv.distance.lm);
cor.corected <- c_vec - idist.fitted.val;
#########################################################################################

data <- bin_neurons(c(1,2),4);

cfit <- fit_copula(data[1,],data[2,], claytonCopula());
gfit <- fit_copula(data[1,],data[2,], gumbelCopula(),upper=10);
rcfit <- fit_copula(data[1,],data[2,], rotCopula(claytonCopula()),upper=10);

test_res <- copula_loglik_test(data[1,],data[2,],gumbelCopula(gfit$par),claytonCopula(cfit$par));
#########################################################################################
# Independance copula loglikes
ind_loglik <- function (x1,x2,cdf1,cdf2) {
  c.density <- pseudo_density
  likelihoods <- mapply(c.density,
                        x1,
                        x2,
                        MoreArgs = list(Fb=cdf1,
                                        Fa=cdf2,
                                        cop=indepCopula(2)
                                        )
                        );
  sum(log(likelihoods));
}

# This takes a long time
m <- neurons[1:10,]
nrow <- nrow(m)
ind_loglik_max <- matrix(NA, nrow=nrow, ncol=nrow);
for (i in 2:nrow) {
  for (j in 1:(i-1)) {
    ind_loglik_max[i,j] <- ind_loglik(m[i,],m[j,],ecdf(m[i,]),ecdf(m[j,]));
  }
}

write.table(ind_loglik_max, "independant10fit-loglik_max.txt")

#########################################################################################
# Let's estimate all possible copulas
pairwise_amh_fit <- pairwise_copula_fit(neurons[1:10,], amhCopula(),inital_val=0)

write.table(pairwise_clayton_fit$par, "clayton10fit-par.txt");
write.table(pairwise_clayton_fit$hessian, "clayton10fit-hessian.txt");
write.table(pairwise_clayton_fit$loglik_max, "clayton10fit-loglik_max.txt");

pairwise_rotClayton_fit <- pairwise_copula_fit(neurons[1:10,],
                                     copula=rotCopula(claytonCopula()),
                                     upper=6)

write.table(pairwise_rotClayton_fit$par, "rot_clayton10fit-par.txt");
write.table(pairwise_rotClayton_fit$hessian, "rot_clayton10fit-hessian.txt");
write.table(pairwise_rotClayton_fit$loglik_max, "rot_clayton10fit-loglik_max.txt");

pairwise_frank_fit <- pairwise_copula_fit(neurons[1:10,],
                                     copula=frankCopula(),
                                     upper=30)

write.table(pairwise_frank_fit$par, "frank10fit-par.txt");
write.table(pairwise_frank_fit$hessian, "frank10fit-hessian.txt");
write.table(pairwise_frank_fit$loglik_max, "frank10fit-loglik_max.txt");

# Can't fit it, for some reason optim goes out of param range
pairwise_normal_fit <- pairwise_copula_fit(neurons[1:10,],
                                     copula=normalCopula)

write.table(pairwise_normal_fit$par, "normal10fit-par.txt");
write.table(pairwise_normal_fit$hessian, "normal10fit-hessian.txt");
write.table(pairwise_normal_fit$loglik_max, "normal10fit-loglik_max.txt");

#Can't fit this eightr optim goes out of range
# t-Copula fit with 4 df
pairwise_t_fit <- pairwise_copula_fit(neurons[1:10,],
                                     copula=tCopula,
                                     inital_val = 0)

write.table(pairwise_t_fit$par, "t10fit-par.txt");
write.table(pairwise_t_fit$hessian, "t10fit-hessian.txt");
write.table(pairwise_t_fit$loglik_max, "t10fit-loglik_max.txt");

# fgm mode
# Doen'st work for some reason :0 
pairwise_fgm_fit <- pairwise_copula_fit(neurons[1:10,],
                                     copula=fgmCopula())

write.table(pairwise_fgm_fit$par, "fgm10fit-par.txt");
write.table(pairwise_fgm_fit$hessian, "fgm10fit-hessian.txt");
write.table(pairwise_fgm_fit$loglik_max, "fgm10fit-loglik_max.txt");



#########################################################################################
# Gumbel seems to do verry well, 6 as an upper bound is working great too. 

pairewise_fit <- pairwise_copula_fit(neurons[1:10,], copula=gumbelCopula(), upper=6)

write.table(pairewise_fit$par, "gumbel10fit-par.txt");
write.table(pairewise_fit$hessian, "gumbel10fit-hessian.txt");
write.table(pairewise_fit$loglik_max, "gumbel10fit-loglik_max.txt");
