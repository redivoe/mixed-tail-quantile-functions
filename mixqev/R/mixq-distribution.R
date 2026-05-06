
# quantile functions
qtypeI <- function(u) {
  -log(1 - u) - (1 - u) / 2
}

qtypeII <- function(u, xi, alpha) {
  alpha + (1 - u)^(-1 / xi) * (1 - (1 - u) / (2 * xi))
}

qtypeIII <- function(u, xi, alpha) {
  alpha - (1 - u)^(1 / xi) * (1 + (1 - u) / (2 * xi))
}


# coefficients of expected order statistics
get_b_typeI <- function(n) {
  digamma(n + 1) - digamma(n:1) - (n:1) / (2 * (n + 1))
}

get_b_typeII <- function(n, xi, alpha) {
  alpha + (-1 + xi * (2 + (1:n / (xi * n + xi - 1)))) / (2 * xi) * exp(lgamma(n + 1) - lgamma(n:1) + lgamma(n:1 - 1 / xi) - lgamma(n + 1 - 1 / xi))
}

get_b_typeIII <- function(n, xi, alpha) {
  alpha - (1 + xi * (3 - 1:n + n + 2 * xi * (n + 1))) / (2 * xi^2) * exp(lgamma(n + 1) - lgamma(n:1) + lgamma(n:1 + 1 / xi) - lgamma(n + 2 + 1 / xi))
}

get_B_mix <- function(n, xi1, xi2, alpha1, alpha2, which_type) {
  out <- list()
  out[[1]] <- rep(1, n)
  if(1 %in% which_type){
    out <- c(out, list(get_b_typeI(n)))
  }
  if(2 %in% which_type){
    out <- c(out, list(get_b_typeII(n, xi1, alpha1)))
  }
  if(3 %in% which_type){
    out <- c(out, list(get_b_typeIII(n, xi2, alpha2)))
  }
  do.call(what = cbind, args = out)
}

# obj function xi
mix_origin2lin <- function(params, which_type){
  theta <- numeric(4)
  params$w <- params$w / sum(params$w)
  theta[1] <- params$mu
  theta[2] <- params$sigma * params$w[1]
  theta[3] <- params$sigma * params$w[2]
  theta[4] <- params$sigma * params$w[3]
  return(theta[c(1, which_type + 1)])
}

obj_fun_xi <- function(xi, sort_e, theta, alpha1, alpha2, which_type) {
  n <- length(sort_e)
  B <- get_B_mix(n = n, xi1 = xi[1], xi2 = xi[2], alpha1 = alpha1, alpha2 = alpha2, which_type)
  sum((sort_e - B %*% theta)^2)
}
obj_fun_xi1 <- function(xi1, sort_e, theta, alpha1, which_type) {
  n <- length(sort_e)
  B <- get_B_mix(n = n, xi1 = xi1, xi2 = NULL, alpha1 = alpha1, alpha2 = NULL, which_type)
  sum((sort_e - B %*% theta)^2)
}
obj_fun_xi2 <- function(xi2, sort_e, theta, alpha2, which_type) {
  n <- length(sort_e)
  B <- get_B_mix(n = n, xi1 = NULL, xi2 = xi2, alpha1 = NULL, alpha2 = alpha2, which_type)
  sum((sort_e - B %*% theta)^2)
}


# obj function for mu and sigma
fit_ls_locscale <- function(sort_e, params, which_type) {
  n <- length(sort_e)
  B <- array(dim = c(n, 2))
  B[, 1] <- 1
  B[, 2] <- 0
  if(1 %in% which_type){
    B[, 2] <- B[, 2] + params$w[1] * get_b_typeI(n)
  }
  if(2 %in% which_type){
    B[, 2] <- B[, 2] + params$w[2] * get_b_typeII(n, params$xi1, params$alpha1)
  }
  if(3 %in% which_type){
    B[, 2] <- B[, 2] + params$w[3] * get_b_typeIII(n, params$xi2, params$alpha2)
  }
  qr_out <- qr(B)
  theta <- qr.coef(qr_out, sort_e)
  return(theta)
}

# obj function for w
obj_fun_w <- function(m, sort_e, mu, sigma, xi1, xi2, alpha1, alpha2, which_type) {
  n <- length(sort_e)
  B <- get_B_mix(n, xi1, xi2, alpha1, alpha2, which_type)
  w <- exp(m) / sum(exp(m))
  theta <- c(mu, sigma * w)
  sum((sort_e - B %*% theta)^2)
}


mixq_pred <- function(u, X, params){
  a <- c(params$mu + X%*%params$beta)
  b <- 0
  if(params$w[1] != 0){
    b <- b + params$w[1] * qtypeI(u)
  }
  if(params$w[2] != 0){
    b <- b + params$w[2] * qtypeII(u, xi = params$xi1, alpha = params$alpha1)
  }
  if(params$w[3] != 0){
    b <- b + params$w[3] * qtypeIII(u, xi = params$xi2, alpha = params$alpha2)
  }
  b <- b * params$sigma
  outer(a, b, FUN = "+")
}

mixq_sim <- function(n, X, params){
  u <- runif(n)

  a <- c(params$mu + X%*%params$beta)
  b <- 0
  if(params$w[1] != 0){
    b <- b + params$w[1] * qtypeI(u)
  }
  if(params$w[2] != 0){
    b <- b + params$w[2] * qtypeII(u, xi = params$xi1, alpha = params$alpha1)
  }
  if(params$w[3] != 0){
    b <- b + params$w[3] * qtypeIII(u, xi = params$xi2, alpha = params$alpha2)
  }
  b <- b * params$sigma
  a + b
}

