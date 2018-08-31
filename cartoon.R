# Created by Oleksandr Sorochynskyi
# On 30 08 18

source('~/IdV/R/conditional.R')

plot_icdf <- function(x) {
  icdf <- make_quantile(x)
  u <- sort(ecdf(x)(x))
  xu <- sort(unique(u))
  
  plot <- ggplot();
  
  for(i in 2:length(xu)) {
    plot <- plot + geom_segment(aes_string(x=xu[i-1],xend=xu[i],y=icdf(xu[i]),yend=icdf(xu[i])),color="black",size=2)
  }
  
  plot <- plot +
    geom_point(mapping= aes(x=xu[-length(xu)],y=icdf(xu)[-1]),size=4) +
    geom_segment(aes(x=0,xend=xu[1],y=0,yend=0),color="black",size=2) +
    xlim(0,1) +
    ylim(min(x),max(x)) +
    #geom_vline(mapping= aes(xintercept= 0)) +
    #geom_vline(mapping= aes(xintercept= 1)) +
    #geom_hline(mapping= aes(yintercept= max(x))) +
    #geom_hline(mapping= aes(yintercept= min(x))) +
    theme(axis.title =element_blank(),
        axis.line= element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none");
  plot;
}


####
# How copula models work illustration (carton)

# Maginal distributions (x2) => marginal cdfs (x2) => pseudo obs     => Empirical copula => Estimation
                    #                                                     copula family  =>

copula <- frankCopula(5);

u <- rCopula(10000, copula)

x1 <- qpois(u[,1], 20)
x2 <- qbinom(u[,2], 10, 0.3)

plot(factor(x1))
plot(factor(x2))

# cdf 
cdf1 <- ecdf(x1);
cdf2 <- ecdf(x2);

# icdf
q1 <- make_quantile(x1);
q2 <- make_quantile(x2);

# pseudo obs 
u1 <- cdf1(x1)
u2 <- cdf2(x2)

# Inverse CDF
plot_icdf(x1)
plot_icdf(x2)


# density + scatter
x <- expand.grid(seq(0,1,0.01),seq(0,1,0.01))
z <- outer(seq(0,1,0.01), seq(0,1,0.01), function(x,y) dCopula(cbind(x,y), copula))

ggplot() +
  xlim(0,1) +
  ylim(0,1) +
  # geom_vline(mapping= aes(xintercept= 0)) +
  # geom_vline(mapping= aes(xintercept= 1)) +
  # geom_hline(mapping= aes(yintercept= 1)) +
  # geom_hline(mapping= aes(yintercept= 0)) +
  stat_contour(mapping= aes(x= x[,1], y= x[,2], z= as.vector(z), fill= ..level..),bins= 30, geom="polygon", color= "white") +
  scale_fill_distiller(palette="Spectral") +
  geom_point(mapping= aes(u[1:100,1], u[1:100,2])) +
  theme(axis.title =element_blank(),
        axis.line= element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none");
  

# Empirical Copula
copula_plot(x1,x2) +
  # geom_vline(mapping= aes(xintercept= 0)) +
  # geom_vline(mapping= aes(xintercept= 1)) +
  # geom_hline(mapping= aes(yintercept= 1)) +
  # geom_hline(mapping= aes(yintercept= 0)) +
  theme(axis.title =element_blank(),
        axis.line= element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none");

ggplot() +
  geom_point(aes(x,y,size=number),data=data.frame(xyTable(x1,x2))) +
  # geom_vline(mapping= aes(xintercept= 0)) +
  # geom_vline(mapping= aes(xintercept= max(x1))) +
  # geom_hline(mapping= aes(yintercept= max(x2))) +
  # geom_hline(mapping= aes(yintercept= 0)) +
  theme(axis.title =element_blank(),
        axis.line= element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none");

  





