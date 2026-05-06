library(evd)
library(mixqev)

# Gumbel
n <- 200
w1 <- matrix(0, 100, 3)
eta1 <- rep(0, 100)
for (i in 1:100) {
  print(paste("1.", i))
  set.seed(i)
  y <- rgumbel(n, loc = 0)
  out <- fit_mixq_ls(y, tol = 1e-05, it_max = 500)
  w1[i, ] <- round(out$w, 2)
  fit <- fgev(y)[[1]]
  eta1[i] <- fit[3]
}

# Frechet
n <- 200
w2 <- matrix(0, 100, 3)
eta2 <- rep(NA, 100)
for (i in 1:100) {
  print(paste("2.",i))
  set.seed(i)
  y <- rfrechet(n, loc = 0, scale = 1, shape = 2)
  out <- fit_mixq_ls(y, tol = 1e-05, it_max = 500)
  w2[i, ] <- round(out$w, 2)
  fit <- try(fgev(y)[[1]])
  if (!is.character(fit)) eta2[i]<-fit[3]
}

# Reversed Weibull
n <- 200
w3 <- matrix(0, 100, 3)
eta3 <- rep(NA, 100)
for (i in 1:100) {
  print(paste("3.",i))
  set.seed(i)
  y <- rrweibull(n, loc = 0, scale = 1, shape = 2)
  out <- fit_mixq_ls(y, tol = 1e-05, it_max = 500)
  w3[i, ] <- round(out$w, 2)
  fit <- try(fgev(y)[[1]])
  if (!is.character(fit)) eta3[i]<-fit[3]
}

# 1/3 mixed quantile
n <- 200
w4 <- matrix(0, 100, 3)
eta4 <- rep(NA, 100)
for (i in 1:100) {
  print(paste("4.",i))
  set.seed(i)
  y <- rgumbel(n, loc=0)+rfrechet(n, loc=0, scale=1, shape=2)+rrweibull(n, loc=5, scale=1, shape=2)
  out <- fit_mixq_ls(y, tol = 1e-05, it_max = 500)
  w4[i, ] <- round(out$w, 2)
  fit<-try(fgev(y)[[1]])
  if (!is.character(fit)) eta4[i]<-fit[3]
}

tab <- t(sapply(list(w1, w2, w3, w4), colMeans)) |>
  round(2)
tab
rowSums(tab)

# [,1]   [,2]   [,3]
# [1,] 0.6664 0.0395 0.2939
# [2,] 0.1892 0.7640 0.0470
# [3,] 0.2615 0.0000 0.7385
# [4,] 0.2871 0.3756 0.3380


