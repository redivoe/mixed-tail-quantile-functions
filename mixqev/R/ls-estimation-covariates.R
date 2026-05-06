
fit_mixq_ls_X <- function(y,
                          X,
                          which_type = c(1, 2, 3),
                          tol = 1e-05,
                          it_max = 100,
                          verbose = FALSE) {

  n <- length(y)
  if (!is.matrix(X)) {
    X <- matrix(X)
    if(nrow(X) != n){
      stop("Number of rows of X is not the same length as y")
    }
  }

  # initialization
  out <- list()

  out$w <- c(0.8, 0.1, 0.1)
  out$w[-which_type] <- 0
  out$w <- out$w / sum(out$w)

  out$beta <- qr.solve(cbind(1, X), y)
  out$mu <- out$beta[1] * 0.7
  out$beta <- out$beta[-1] * 0.7
  out$sigma <- sd(y - cbind(1, X) %*% c(out$mu, out$beta))
  if(2 %in% which_type){
    out$alpha1 <- quantile(y, 0.05)
  }
  if(3 %in% which_type){
    out$alpha2 <- quantile(y, 0.95)
  }

  obj <- numeric(length = it_max)
  ratio <- 1
  h <- 0


  while ((ratio > tol) & (h < it_max)) {
    h <- h + 1

    e <- y - X %*% out$beta
    sort_e <- sort(e)
    index <- order(e)

    #1) estimate xi
    if(all(c(2, 3) %in% which_type)){
      out_xi <- optim(par = c(1.5, 0.5),
                      fn = obj_fun_xi,
                      lower = c(1 + 1e-4, 1e-4),
                      upper = c(3, 2),
                      method = "L-BFGS-B",
                      sort_e = sort_e,
                      theta = mix_origin2lin(out, which_type),
                      alpha1 = out$alpha1,
                      alpha2 = out$alpha2,
                      which_type = which_type)
      out$xi1 <- out_xi$par[1]
      out$xi2 <- out_xi$par[2]
    }else if(2 %in% which_type){
      out_xi1 <- optimise(f = obj_fun_xi1,
                          interval = c(1 + 1e-4, 3),
                          sort_e = sort_e,
                          theta = mix_origin2lin(out, which_type),
                          alpha1 = out$alpha1,
                          which_type = which_type)
      out$xi1 <- out_xi1$minimum
    }else if(3 %in% which_type){
      out_xi2 <- optimise(f = obj_fun_xi2,
                          interval = c(1e-4, 2),
                          sort_e = sort_e,
                          theta = mix_origin2lin(out, which_type),
                          alpha2 = out$alpha2,
                          which_type = which_type)
      out$xi2 <- out_xi2$minimum
    }

    # 2) estimate mu and sigma
    out_locscale <- fit_ls_locscale(sort_e = sort_e, params = out, which_type = which_type)
    out$mu <- out_locscale[1]
    out$sigma <- out_locscale[2]


    # 3) estimate weights
    if(length(which_type) > 1){
      out_w <- optim(par = rep(1, length(which_type)),
                     fn = obj_fun_w,
                     sort_e = sort_e,
                     mu = out$mu,
                     sigma = out$sigma,
                     xi1 = out$xi1,
                     xi2 = out$xi2,
                     alpha1 = out$alpha1,
                     alpha2 = out$alpha2,
                     which_type = which_type)
      out$w[which_type] <- exp(out_w$par) / sum(exp(out_w$par))
    }

    #4) estimate beta
    B <- get_B_mix(n, out$xi1, out$xi2, out$alpha1, out$alpha2, which_type)
    theta <- mix_origin2lin(out, which_type)
    y_res <- y[index] - B %*% theta

    qr_out <- qr(cbind(1, X[index, ]))
    out$beta <- matrix(qr.coef(qr_out, y_res)[-1], ncol = 1)

    #5)
    obj[h] <- sum((y[index] - X[index, ] %*% out$beta - B %*% theta) ^ 2)
    if (h > 10){
      ratio <- (obj[h - 1] - obj[h]) / obj[h - 1]
    }

    if (isTRUE(verbose)) {
      cat(sprintf("Iter %3d - Loss: %.6f\n", h, obj[h]))
    }

  }

  out$obj <- obj[1:h]
  out$which_type <- which_type

  return(out)
}
