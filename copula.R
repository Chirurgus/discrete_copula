# created by Oleksandr Sorochynskyi
# on 7 06 18

library(copula)
library(ggplot2)

source('~/IdV/R/open_neurons.R')

########################################################################################
# Utility

copula_distributional_transfrom <- function(x) {
  cdf <- ecdf(x);
  cdf(x-1) + runif(length(x)) * (cdf(x) - cdf(x-1));
}

copula_plot_cop_density <- function(copula, x1, x2) {
  stopifnot(length(x1) == length(x2));
  
  x1 <- ecdf(x1)(x1);
  x2 <- ecdf(x2)(x2);
  
  #plot data
  x <- sort(unique(c(0,x1,1)));
  y <- sort(unique(c(0,x2,1)));
  count <- matrix(0,ncol=length(x),nrow=length(y));
  for (i in 2:length(x)) {
    for (j in 2:length(y)) {
      count[j,i] <- prob(copula, c(x[i-1],y[j-1]), c(x[i],y[j]))
    }
  }
  pos_x <- matrix(x, nrow=nrow(count), ncol=ncol(count), byrow=TRUE)
  pos_y <- matrix(y, nrow=nrow(count), ncol=ncol(count), byrow=FALSE)
  
  xmin <- pos_x[-nrow(pos_x),-ncol(pos_x)];
  ymin <- pos_y[-nrow(pos_y),-ncol(pos_y)];
  xmax <- pos_x[-1,-1];
  ymax <- pos_y[-1,-1];
  
  ret <- ggplot(data.frame(), aes()) +
   geom_rect(
      aes(
        fill= log(as.vector(count[-1,-1])),
        xmax= as.vector(xmax),
        xmin= as.vector(xmin),
        ymax= as.vector(ymax),
        ymin= as.vector(ymin)
      )
    ) +
    labs(fill = "Log density  ") +
    scale_fill_gradientn(colors=c("blue", "red"),
                        values=scales::rescale(c(-25, -5, -3, 0)),
                        limits=c(-25, 0)
                        );
  return(ret);

}

copula_plot_circle <- function(x1, x2, coef=0.01) {
  x1 <- ecdf(x1)(x1)
  x2 <- ecdf(x2)(x2)
  tbl <- xyTable(x1,x2);
  size <- tbl$number*coef; 
  plot(tbl$x, tbl$y, cex=size, xlim= c(0,1), ylim= c(0,1));
}

copula_plot_diff <- function(x1, x2, copula) {
  stopifnot(length(x1) == length(x2));
  
  # samp <- copula_rand_emp(length(x1), copula, x1 = x1, x2=x2);
  # s1 <- ecdf(samp[,1])(samp[,1]);
  # s2 <- ecdf(samp[,2])(samp[,2]);
  
  x1 <- ecdf(x1)(x1);
  x2 <- ecdf(x2)(x2);
 
  # Random data generated with copula 
  # tbl_samp <- xyTable(s1,s2,digits=7);
  # x_samp <- sort(unique(c(0,s1,1)));
  # y_samp <- sort(unique(c(0,s2,1)));
  # count_samp <- matrix(0,ncol=length(x_samp),nrow=length(y_samp));
  # for (i in seq_along(x_samp)) {
  #   for (j in seq_along(y_samp)) {
  #     index <- abs(tbl_samp$x - x_samp[i]) < 1e-7 & abs(tbl_samp$y - y_samp[j]) < 1e-7;
  #     count_samp[j,i] <- max(0,tbl_samp$number[index]);
  #   }
  # }
  # Actual data
  tbl <- xyTable(x1,x2,digits=7);
  x <- sort(unique(c(0,x1,1)));
  y <- sort(unique(c(0,x2,1)));
  freq <- matrix(0,ncol=length(x),nrow=length(y));
  dens <- matrix(0,ncol=length(x),nrow=length(y));
  for (i in 2:length(x)) {
    for (j in 2:length(y)) {
      index <- abs(tbl$x - x[i]) < 1e-7 & abs(tbl$y - y[j]) < 1e-7;
      freq[j,i] <- max(0,tbl$number[index])/length(x1);
      dens[j,i] <- prob(copula, c(x[i-1],y[j-1]), c(x[i],y[j]))
    }
  }
  pos_x <- matrix(x, nrow=nrow(freq), ncol=ncol(freq), byrow=TRUE)
  pos_y <- matrix(y, nrow=nrow(freq), ncol=ncol(freq), byrow=FALSE)
  
  # Sanity check
  stopifnot(sum(freq) == 1, sum(dens) == 1);

  xmin <- pos_x[-nrow(pos_x),-ncol(pos_x)];
  ymin <- pos_y[-nrow(pos_y),-ncol(pos_y)];
  xmax <- pos_x[-1,-1];
  ymax <- pos_y[-1,-1];
  
  ret <- ggplot(data.frame(), aes()) +
   geom_rect(
      aes(
        #fill= log((as.vector(count[-1,-1])+1) / (as.vector(count_samp[-1,-1])+1)),
        fill= log(as.vector(freq[-1,-1]) / as.vector(dens[-1,-1])),
        xmax= as.vector(xmax),
        xmin= as.vector(xmin),
        ymax= as.vector(ymax),
        ymin= as.vector(ymin)
      )
    ) +
    scale_fill_gradientn(colors=c('blue','white','red'),
                         values=scales::rescale(c(-5,-.5,0,.5,4)),
                         limits=c(-5,4)
                         )+
    labs(fill = "Log density ratio  ");
  return(ret);
}

# x1 and x2 will be transformed to pseudo obs
copula_plot <- function(x1, x2) {
  stopifnot(length(x1) == length(x2));
  
  x1 <- ecdf(x1)(x1);
  x2 <- ecdf(x2)(x2);
  
  #plot data
  tbl <- xyTable(x1,x2,digits=7);
  x <- sort(unique(c(0,x1,1)));
  y <- sort(unique(c(0,x2,1)));
  count <- matrix(0,ncol=length(x),nrow=length(y));
  colnames(count) <- x;
  rownames(count) <- y;
  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      index <- abs(tbl$x - x[i]) < 1e-7 & abs(tbl$y - y[j]) < 1e-7;
      count[j,i] <- max(0,tbl$number[index])/length(x1);
    }
  }
  pos_x <- matrix(x, nrow=nrow(count), ncol=ncol(count), byrow=TRUE)
  pos_y <- matrix(y, nrow=nrow(count), ncol=ncol(count), byrow=FALSE)
  
  # Sanity check
  stopifnot(sum(count) == 1);
  
  xmin <- pos_x[-nrow(pos_x),-ncol(pos_x)];
  ymin <- pos_y[-nrow(pos_y),-ncol(pos_y)];
  xmax <- pos_x[-1,-1];
  ymax <- pos_y[-1,-1];
  
  ret <- ggplot(data.frame(), aes()) +
   geom_rect(
      aes(
        fill= log(as.vector(count[-1,-1])),
        xmax= as.vector(xmax),
        xmin= as.vector(xmin),
        ymax= as.vector(ymax),
        ymin= as.vector(ymin)
      )
    ) +
    labs(fill = "Log frequency  ") +
    scale_fill_gradientn(colors=c("blue", "red"),
                        values=scales::rescale(c(-25, -5, -3, 0)),
                        limits=c(-25, 0)
                        );
  return(ret);
}

# test soubld be a function taking x1 x2, cop0 cop2
# type decides if it tests for type 1 or type 2 error
estim_test_error <- function(rep, nobs, x1, x2, cop0, cop1, test, type =1) {
  nobs <- length(x1); 
 
  # Display progresss, since this can be quite slow 
  pb <- txtProgressBar(min=0,max=rep,style=3);
 
  q1 <- make_quantile(x1);
  q2 <- make_quantile(x2);
  signif <- replicate(rep, {
    if (type == 1) {
      data <- copula_rand_emp(nobs,cop0,q1,q2);
    }
    else {
      data <- copula_rand_emp(nobs,cop1,q1,q2);
    }
    ret <- test(data[,1],data[,2],cop0,cop1)$signif;
    
    setTxtProgressBar(pb,getTxtProgressBar(pb)+1);
    
    ret;
  });
  
  close(pb)
  
  p <- length(which(signif==TRUE))/nrep
  ifelse(type==1,p,1-p)
}

copula_loglik_vector <- function(x1, x2, cop, cdf1, cdf2) {
  if (missing(cdf1)) cdf1 = ecdf(x1);
  if (missing(cdf2)) cdf2 = ecdf(x2);
  
  
  u1 <- cdf1(x1);
  u2 <- cdf2(x2);
  break.points1 <- sort(unique(c(0,u1,1)));
  break.points2 <- sort(unique(c(0,u2,1)));

  c.density <- make_pseudo_density(cop, break.points1, break.points2);
  likelihoods <- mapply(c.density,u1,u2);
  log(likelihoods);
}

copula_loglik <- function(...) {
  sum(copula_loglik_vector(...));  
}

make_quantile <- function(x) {
  # Break points
  u <- unique(c(0,cumsum(table(x)/length(x))))
  x <- sort(unique(x));
  
  function(p) {
    index <- findInterval(p, u,
                          left.open= TRUE,
                          rightmost.closed= TRUE,
                          all.inside= TRUE)
    x[index]
  }
}

copula_rand_cond <- function(u1, u2, cop, x2, x1, q1, q2) {
  if (missing(q1)) {
    if (missing(x1)){
      stop("copula_rand_emp needs either the quantile function or x");
    }
    q1 <- make_quantile(x1);
  }
  if (missing(q2)) {
    if (missing(x2)){
      stop("copula_rand_emp needs either the quantile function or x");
    }
    q2 <- make_quantile(x2);
  }
  
  rand <- cCopula(cbind(u1, u2), copula = cop, inverse = TRUE)
  cbind(q1(rand[,1]),q2(rand[,2]));
}

# Random numer generator using emprical cdf, and quantile
copula_rand_emp <- function(n, cop, q1, q2, x1, x2) {
  if (missing(q1)) {
    if (missing(x1)){
      stop("copula_rand_emp needs either the quantile function or x");
    }
    q1 <- make_quantile(x1);
  }
  if (missing(q2)) {
    if (missing(x2)){
      stop("copula_rand_emp needs either the quantile function or x");
    }
    q2 <- make_quantile(x2);
  }
  
  rand <- rCopula(n,cop);
  cbind(q1(rand[,1]),q2(rand[,2]));
}

#######################################################################################
# Kramer von Misses type test (by Copula independance test for bivaiate descrete data)

Sn_statistic <- function(x1,x2){
  ########OPTION WITHOUT ZERO ROWS AND COLUMNS#######
  #If data comes in a matrix form (instead of two columns) comment
  #out the first line and replace it with Data <-data
  Data<-table(x1,x2)
  rows<-nrow(Data)
  cols<-ncol(Data)
  
  n<-sum(Data)
  T<-apply(t(apply(Data,1,cumsum)),2,cumsum)
  C<-cbind(rep(0,rows+1),rbind(rep(0,cols),T/n))
  h<-cbind(rep(0,rows+1),rbind(rep(0,cols),Data/n))
  #Marginals:
  F<-c(0,T[,cols]/n)
  G<-c(0,T[rows,]/n)
  
  #Using dudv
  matrix<- C-F%o%G
  ab.bar<-(F-c(0,F[1:(length(F)-1)]))%o%(G-c(0,G[1:(length(G)-1)]))
  K1<-matrix[2:(rows+1),2:(cols+1)]
  K2<-matrix[1:rows,2:(cols+1)]
  K3<-matrix[2:(rows+1),1:cols]
  K4<-matrix[1:rows,1:cols]
  Sn.2<-sum((1/18)*(ab.bar[2:(rows+1),2:(cols+1)]) *(2*K1^2+2*K2^2+2*K3^2+2*K4^2+2*K1*K2+2*K1*K3+K1*K4+K2*K3+2*K2*K4+2*K3*K4))
  
  Sn.2;
}
  
copula_kramer_indep_test <- function(x1,x2,alpha=0.05,nrep=200) {
  q1 <- make_quantile(x1)
  q2 <- make_quantile(x2)
  
  Sn_sample <- replicate(nrep, {
    t <- copula_rand_emp(length(x1),indepCopula(),q1,q2)
    Sn_statistic(t[,1],t[,2])
  });
  
  crit <- quantile(Sn_sample, 1-alpha);
  Sn <- Sn_statistic(x1,x2);
  
  list(crit.value=crit,
       Sn.stat = Sn,
       p.value= 1- ecdf(Sn_sample)(Sn),
       signif=Sn>=crit
       )
}

#######################################################################################
# Euclidian distance test

copula_distance <- function(x1, x2, copula, cdf1, cdf2, empir) {
  if (missing(cdf1)) cdf1 <- ecdf(x1);
  if (missing(cdf2)) cdf2 <- ecdf(x2);
  if (missing(empir)) empir <- C.n(cbind(cdf1(x1),cdf2(x2)),cbind(cdf1(x1),cdf2(x2)))
  u <- cbind(cdf1(x1),cdf2(x2))
  estim <- pCopula(u,copula)
  sum((empir-estim)^2);
  
}
# compares the given copula to the emprirical copula
copula_distance_test <- function(x1, x2, copula0, copula1, cdf1, cdf2) {
  if (missing(cdf1)) cdf1 <- ecdf(x1);
  if (missing(cdf2)) cdf2 <- ecdf(x2);
  empir <- C.n(cbind(cdf1(x1),cdf2(x2)),cbind(cdf1(x1),cdf2(x2)))
  diff0 <- copula_distance(x1, x2, copula0, cdf1, cdf2, empir);
  diff1 <- copula_distance(x1, x2, copula1, cdf1, cdf2, empir);
  res <- diff0 <= diff1;
  list(choice=ifelse(res, "Copula0", "Copula1"),
       signif=!res
  );
}

#######################################################################################
# Pillow independence tests (from Characterizing neural dependencies)

copula_pillow_lld_dist <- function(x1, x2, copula) {
  q1 <- make_quantile(x1);
  q2 <- make_quantile(x2);
  
  samp <- replicate(20, {
    rand <- copula_rand_emp(length(x1), indepCopula(), x1= x1, x2= x2);
    indep_loglik <- mean(copula_loglik_vector(rand[,1], rand[,2], indepCopula()));
    model_loglik <- fit_copula(rand[,1], rand[,2], copula)$loglik_max/length(x1);
    model_loglik - indep_loglik;
  });
  
  samp;
}

copula_pillow_ll_threshold <- function(x1, x2, copula, lld_dist) {
  ll <- copula_log_lik_maker(x1, x2, ecdf(x1), ecdf(x2), copula);
  
  indep_loglik <- mean(copula_loglik_vector(x1, x2, indepCopula()));
  
  threshold <- quantile(lld_dist, 0.95);
  
  uniroot(function(theta) indep_loglik - ll(theta)/length(x1) - threshold,interval=c(-20,20)) ;
}

copula_pillow_test <- function(x1, x2, copula, lld_dist)
{
  if (missing(lld_dist)) {
    lld_dist <- copula_pillow_lld_dist(x1,x2,copula);
  }
  
  indep_loglik <- mean(copula_loglik_vector(x1, x2, indepCopula()));
  model_loglik <- fit_copula(x2, x2, copula)$loglik_max/length(x1);
  lld <- model_loglik - indep_loglik;
  
  list(signif=quantile(lld_dist,0.95) < lld,
       p.value= 1-ecdf(lld_dist)(lld));
}

#######################################################################################
# Likelihood ratio test (by your one and only me(with a little help of Pearson and Neyman)
  # And aslo Vuong (something likethat)

copula_loglik_test <- function(x1,x2,copula0,copula1,cdf1,cdf2,alpha,one.sided= FALSE) {
  if (missing(cdf1)) cdf1 <- ecdf(x1);
  if (missing(cdf2)) cdf2 <- ecdf(x2);
  if (missing(alpha)) alpha = 0.05;

  ll0 <- copula_loglik_vector(x1,x2,copula0,cdf1,cdf2);
  ll1 <- copula_loglik_vector(x1,x2,copula1,cdf1,cdf2);
  
  dll <- ll1 - ll0;
  mu_hat <- mean(dll);
  o_hat <- sum(dll^2)/length(dll) - mu_hat^2
  lr <- mu_hat/o_hat;
 
  if (one.sided) { 
    threshold <- qnorm(1-alpha);
  }
  else {
    threshold <- qnorm(1-alpha/2);
  }
  if (one.sided) {
    signif <- lr > threshold;
  }
  else {
    signif <- (abs(lr) > threshold);
  }
  p.value <- NULL;
  
  if (one.sided) {
    p.value <- pnorm(lr);
  }
  else {
    if (signif && lr > threshold) {
      p.value <- 2*pnorm(lr);
    }
    if (signif && lr < -threshold) {
      p.value <- 2*pnorm(lr,lower.tail = TRUE);
    }
  }
  
  list(signif= signif,lr=lr, t=threshold, p.value= p.value);
}
   

# copula_loglik_test <- function(x1,x2,copula0,copula1,cdf1,cdf2,alpha,nrep) {
#   if (missing(cdf1)) cdf1 <- ecdf(x1);
#   if (missing(cdf2)) cdf2 <- ecdf(x2);
#   if (missing(alpha)) alpha = 0.05;
#   if (missing(nrep)) nrep = 200;
#   
#   q1 <- make_quantile(x1);
#   q2 <- make_quantile(x2);
#   loglik_s <- replicate(nrep, {
#     s <- copula_rand_emp(length(x1), copula0, q1, q2);
#     ll <- sapply(c(copula1,copula0), function(cop) copula_loglik(s[,1],s[,2],cop,cdf1,cdf2))
#     ll[1] - ll[2]; 
#   });
#  
#   loglik1 <- copula_loglik(x1,x2,copula1,cdf1,cdf2) 
#   loglik0 <- copula_loglik(x1,x2,copula0,cdf1,cdf2)
#   loglik_dif <- loglik1 - loglik0;
#   
#   c <- quantile(loglik_s, 1-alpha);
#   
#   pvalue <- 1- ecdf(loglik_s)(loglik_dif); 
#   res <- c <= loglik_dif;
#   
#   list(llik_H0=loglik0, llik_H1=loglik1, c=c, llik_dif=loglik_dif,
#        result=ifelse(res,"H1","H0"),sign_value=alpha,"p-value"=pvalue,
#        signif=res);
# }

#######################################################################################
# Copula estimation

# Optimized (memoized version)
make_pseudo_density <- function(cop, u1, u2) {
  m <- list();
  u1 <- unique(u1);
  u2 <- unique(u2);
  if (is.unsorted(u1)) u1 <- sort(u1);
  if (is.unsorted(u2)) u2 <- sort(u2);
  function(a,b) {
    index <- paste(as.character(a),as.character(b));
    if (is.null(m[[index]])) { 
      index1 <- which(u1 == a)
      index2 <- which(u2 == b)
      lower <- c(u1[index1 - 1], u2[index2 - 1]);
      upper <- c(u1[index1], u2[index2]);
      m[[index]] <<- prob(cop,lower,upper);
    }
    m[[index]];
  }
}

copula_log_lik_maker <- function(x1, x2, cdf1, cdf2, copula) {
  nobs <- length(x1)
  
  u1 <- cdf1(x1);
  u2 <- cdf2(x2);
  break.points1 <- sort(unique(c(0,u1,1)));
  break.points2 <- sort(unique(c(0,u2,1)));
  
  function (theta) {
    copula <- setTheta(copula, theta);
    c.density <- make_pseudo_density(copula, break.points1, break.points2);
    likelihoods <- mapply(c.density,u1,u2);
    sum(log(likelihoods));
  }
}

fit_copula <- function(x1, x2, copula, cdf1, cdf2, inital_val, lower, upper, hessian= FALSE) {
  if (missing(cdf1)) cdf1 <- ecdf(x1);
  if (missing(cdf2)) cdf2 <- ecdf(x2);
  if (is(copula,"indepCopula")) {
    return(list(par= NA, hessian = NA, loglik_max= copula_loglik(x1,x2,copula,cdf1,cdf2),cop=copula));
  }
  if (missing(inital_val)) inital_val <- iRho(copula, cor(x1,x2,method="s"));
  if (missing(lower)) {
    if (is(copula,"rotCopula")) lower = max(-20, copula@copula@param.lowbnd)
    else lower = max(-20, copula@param.lowbnd)
  }
  if (missing(upper)) {
    if (is(copula, "rotCopula")) upper = min(20, copula@copula@param.upbnd)
    else upper = min(20, copula@param.upbnd)
  }
  
  copula_loglik_opt <- copula_log_lik_maker(x1,x2,cdf1, cdf2, copula)
  fit <- optim(inital_val,function(x) -copula_loglik_opt(x), method="Brent", lower=lower, upper=upper, hessian= hessian)
  if (fit$convergence) {
    warning("Optimization may not have converged: non zero convergencce code from optim()",
              paste("\noptim(): ", fit$message));
  }
  list(par= fit$par,
       hessian = fit$hessian,
       loglik_max= -fit$value,
       cop= setTheta(copula, fit$par)); 
}
