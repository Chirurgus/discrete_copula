# Created by Oleksandr Sorochynskyi
# On 19 07 18

################################################################################
#

m <- bin_neurons(c(1,2), 1, 0);
x <- ifelse(m[1,] > 0, 1, 0);
y <- ifelse(m[2,] > 0, 1, 0);

x.full <- function(x,lag) {
  stopifnot(lag > 0);
  index <- length(x) - 0:(lag-1);
  x[-index]
}
x.lag <- function(x,lag) {
  stopifnot(lag > 0);
  index <- 1:lag
  x[-index]
}

# estimate marginals with logit regression
x.glm <- glm(x.full(x,1) ~ x.lag(x,1), family=binomial)
summary(x.glm);
