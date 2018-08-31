# Created by Oleksandr Sorochynskyi
# On 14 06 18

source("open_neurons.R")
source("copula.R");


library(ggplot2)


#######################################################################################
# Scatterplots for copulas


my.theme <- theme(text= element_text(size=25));

rand<- rCopula(1000, indepCopula());
ggplot(data=data.frame(rand=rand), aes(x=rand[,1], y=rand[,2])) +
  geom_point() +
  my.theme;

dev.copy2pdf(file="plots/test.pdf", height=5, width=5);


#######################################################################################
# Loglikelihoods, bootstraped for one pair, multiple models

boot_loglik <- function(nrep,x1,x2,copula) {
  replicate(nrep, {
    samp <- sample(1:length(x1),length(x1),replace = TRUE);
    copula_loglik(x1[samp],x2[samp],copula);
    }
  );
}

hist_loglik_across_models <- function(x1,x2,cops,cop.options) {
  train_i <- sample(1:length(x1), floor(length(x1)/2));
  ll <- list();
  for (i in 1:length(cops)) {
    fit <- tryCatch({
      options <- cop.options[[i]];
      options[["x1"]] <- x1[train_i];
      options[["x2"]] <- x2[train_i];
      options[["copula"]] <- cops[[i]];
      do.call(fit_copula,options);
    },error=function(err) {
      warning(paste("Coult not fit copula: ", class(cops[[i]])));
      NULL;
      }
    );
    
    if (!is.null(fit)) {
      ll[[i]] <- 
        list(loglik=boot_loglik(100,x1,x2,setTheta(cops[[i]], fit$par)),
             name=class(cops[[i]])
        );
    }
  }
  
  xlim.min = 1.1*min(sapply(ll,function(l) min(l$loglik,na.rm=TRUE)),na.rm=TRUE);
  xlim.max = 0.9*max(sapply(ll,function(l) max(l$loglik,na.rm=TRUE)),na.rm=TRUE);
  xlim = c(xlim.min,xlim.max); 
  
  par(mfrow=c(length(ll),1));
 
  lapply(ll, function(l) {
      tryCatch(hist(l$loglik,xlim=xlim,main=l$name,col=sample(colors(),1)),
               error=function(x){}
               );
    }
  ); 
  return();
}

data <- bin_neurons(c(3,2),4);
hist_loglik_across_models(data[1,],data[2,],
                          list(
                            indepCopula(),
                            gumbelCopula(),
                            claytonCopula(),
                            frankCopula(),
                            normalCopula(),
                            rotCopula(claytonCopula()),
                            rotCopula(gumbelCopula())
                            ),
                          list(list(),
                               list(upper=6),
                               list(),
                               list(),
                               list(),
                               list(upper=6),
                               list()
                              )
                          );

#######
# Analyse the distributions of booted loglikes
gfit <- fit_copula(data[1,],data[2,],gumbelCopula(),upper=6);
cfit <- fit_copula(data[1,],data[2,],claytonCopula(),upper=6);
rotcfit <- fit_copula(data[1,],data[2,],rotCopula(claytonCopula()),upper=6);

gcop <- setTheta(gumbelCopula(),gfit$par);
ccop <- setTheta(claytonCopula(),cfit$par);
crcop <- setTheta(rotCopula(claytonCopula()),rotcfit$par);

gll <- boot_loglik(100,data[1,],data[2,],gcop)
rcll <- boot_loglik(100,data[1,],data[2,],crcop)
cll <- boot_loglik(100,data[1,],data[2,],ccop)

#######################################################################################
# Plot loglikelihood distributins from different copula families

independance_m <- as.matrix(read.table("independant10fit-loglik_max.txt"))
gumbel_m <- as.matrix(read.table("gumbel10fit-loglik_max.txt"));
rot_clayton_m <- as.matrix(read.table("rot_clayton10fit-loglik_max.txt"));
frank_m <- as.matrix(read.table("frank10fit-loglik_max.txt"))
normal_m <- as.matrix(read.table("normal10fit-loglik_max.txt"))
t_m <- as.matrix(read.table("t10fit-loglik_max.txt"))

hist(na.omit(as.vector(gumbel_m)),xlim=c(-70000,-10000))
hist(na.omit(as.vector(rot_clayton_m)),add=TRUE,col="red")
hist(na.omit(as.vector(frank_m)),add=TRUE,col="black")
hist(na.omit(as.vector(normal_m)),add=TRUE,col="green")
hist(na.omit(as.vector(t_m)),add=TRUE,col="blue")
hist(na.omit(as.vector(independance_m)),add=TRUE,col="yellow")

# Conclusion: they are basically the same, except independace, that is bad
#             half the time it's just -Inf sooooooooo, yeah

#######################################################################################
# plot loglik
x1 <- neurons[1,]
x2 <- neurons[2,]
F1 <- margin_empir(x1)$cdf;
F2 <- margin_empir(x2)$cdf;
cop <- gumbelCopula

loglik <- copula_log_lik_maker(x1,x2,F1,F2,cop);
x <- seq(1,6,by=0.01);
y <- lapply(x,loglik)
plot(x,y)

######################################################################################
# empirical copula
library(plotrix)

data <- bin_neurons(c(1,2), 20);
t1 <- data[1,]
t2 <- data[2,]
x <- seq(0,1,by=0.01);
y <- x;
z <- matrix(NA, ncol = length(x), nrow= length(y));
for (i in seq_along(x)) {
  for (j in seq_along(y)) {
    z[i,j] <- C.n(cbind(x[i],y[j]), cbind(pobs(t1),pobs(t2)));
  }
}
persp(x,y,z,theta=0,main="Empirical copula");
color2D.matplot(z,main="Empirical copula");

tb <- xyTable(t1,t2)
#tb <- lapply(tb, function(x) x/length(t1));
plot(tb$x,tb$y, cex=(tb$number-tb$x*tb$y)*100)

cdf1 <- ecdf(t1);
cdf2 <- ecdf(t2);
tbu <- xyTable(cdf1(t1),cdf2(t2));
plot(tbu$x,tbu$y,cex=tbu$number*0.05,xlim=c(0,1),ylim=c(0,1));

#######################################################################################
# 3d graps of fitted copulas
library(plotrix)

t1 <- neurons[9,]
t2 <- neurons[10,]

ffit <- fit_copula(t1,t2, frankCopula());
cfit <- fit_copula(t1,t2, claytonCopula());
rcfit <- fit_copula(t1,t2, rotCopula(claytonCopula()),upper=6)

lapply(c(frankCopula(ffit$par),
         claytonCopula(cfit$par),
         rotCopula(claytonCopula(rcfit$par))),
       function(cop) {
          x <- seq(0,1,by=0.01);
          y <- x;
          z <- matrix(NA, ncol = length(x), nrow= length(y));
          for (i in seq_along(x)) {
            for (j in seq_along(y)) {
              z[i,j] <- dCopula(c(x[i], y[j]), cop);
            }
          }
          persp(x,y,z,theta=0,main=paste(class(cop)," density"));
          color2D.matplot(z,main=paste(class(cop)," density"));
       });

#######################################################################################
# Pretty copula graphs

plot(x=indepCopula(),n=1000,main="Independance copula");
plot(x=claytonCopula(3),n=1000,main="Clayton copula");
plot(x=gumbelCopula(2),n=1000,main="Gumbel copula");
plot(x=rotCopula(claytonCopula(4)),n=1000,main="Rotated clayton copula");
plot(x=frankCopula(0.1),n=10000,main="Frank copula");
plot(x=normalCopula(0.7),n=1000,main="Normal copula");
plot(x=tCopula(0.7),n=1000, main= "t copula")
plot(x=fgmCopula(1),n=2000, main= "FGM copula");
plot(x=)