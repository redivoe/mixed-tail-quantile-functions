

fit_mixq_ls <- function(y, tol = 1e-05, it_max = 500, upper_bound_xi1 = 3, upper_bound_xi2 = 2) {
  n <- length(y)
  sort_y <- sort(y)

  ## initialization
  mu <- mean(y)
  sigma <- sd(y)
  w <- c(1 / 3, 1 / 3, 1 / 3)

  alpha1 <- quantile(y, 0.05)
  alpha2 <- quantile(y, 0.95)
  which_type <- c(1, 2, 3)

  obj <- numeric(length = it_max)
  ratio <- 1
  iter <- 0

  params <- list(mu = mu, sigma = sigma, w = w, alpha1 = alpha1, alpha2 = alpha2)

  while ((ratio > tol) & (iter < it_max)) {
    iter <- iter + 1

    # 1) estimate xi
    out_xi <- optim(par = c(1.5, 0.5),
                    fn = obj_fun_xi,
                    lower = c(1 + 1e-4, 1e-4),
                    upper = c(upper_bound_xi1, upper_bound_xi2),
                    method = "L-BFGS-B",
                    sort_e = sort_y,
                    theta = mix_origin2lin(params, which_type),
                    alpha1 = params$alpha1,
                    alpha2 = params$alpha2,
                    which_type = which_type)
    params$xi1 <- out_xi$par[1]
    params$xi2 <- out_xi$par[2]

    # 2) estimate location and scale
    out_locscale <- fit_ls_locscale(sort_e = sort_y, params = params, which_type = which_type)
    params$mu <- out_locscale[1]
    params$sigma <- out_locscale[2]

    # 3) estimate weights
    out_w <- optim(par = rep(1, length(which_type)),
                   fn = obj_fun_w,
                   sort_e = sort_y,
                   mu = params$mu,
                   sigma = params$sigma,
                   xi1 = params$xi1,
                   xi2 = params$xi2,
                   alpha1 = params$alpha1,
                   alpha2 = params$alpha2,
                   which_type = which_type)
    params$w <- exp(out_w$par) / sum(exp(out_w$par))

    #4) evaluate objective function
    theta <- mix_origin2lin(params, which_type)
    B <- get_B_mix(n, params$xi1, params$xi2, params$alpha1, params$alpha2, which_type)
    obj[iter] <- sum((sort_y - B %*% theta) ^ 2)

    if (iter > 10){
      ratio <- (obj[iter - 1] - obj[iter]) / obj[iter - 1]
    }
  }

  params$obj <- obj[1:iter]

  return(params)
}


