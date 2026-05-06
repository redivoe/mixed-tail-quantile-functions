library(Qtools)
library(quantreg)
library(qrcm)
library(evd)
library(EnvStats)
library(mev)
library(evgam)
library(furrr)
library(purrr)

# install.packages("mixqev_0.1.0.tar.gz", path = NULL)
library(mixqev)

# functions for nonparametric estimation of conditional quantile functions
# (see Note 1 of the manuscript)
get_ccdf <- function(fx, n_points = 1e4){
  x <- seq(0, 1, length.out = n_points)
  y <- fx(x)
  return(list(
    out_cmid_ecdf = Qtools::cmidecdf(y ~ x, ecdf_est = 'logit'),
    x = x
  ))
}

get_q_true <- function(x, prob, true_ccdf){
  n <- length(x)
  Fhat <- matrix(NA, n, length(true_ccdf$x))
  for (i in 1:n) {
    idx_i <- which.min(abs(x[i] - true_ccdf$x))[1]
    Fhat[i, ] = true_ccdf$out_cmid_ecdf$Fhat[idx_i, ]
  }
  yo <- true_ccdf$out_cmid_ecdf$yo
  q_true <- yo[sapply(apply(Fhat >= prob, 1, which), min)]
  return(q_true)
}


set.seed(1)
fx <- list("pareto" = \(x) 3 * x + EnvStats::rpareto(length(x), location = 1, shape = 2) - 1 ,
           "logn" = \(x) 3 * x + rlnorm(length(x)),
           "egp" = \(x) 3 * x + mev::regp(n = length(x), scale = 1, shape = 1/10, kappa = 1.5))

qf <- list("pareto" = \(x, u) 3 * x + EnvStats::qpareto(u, location = 1, shape = 2) - 1 ,
           "logn" = \(x, u) 3 * x + qlnorm(p = u),
           "egp" = \(x, u) 3 * x + mev::qegp(p = u, scale = 1, shape = 1/10, kappa = 1.5))

dgp <- names(fx)

reps <- 100
n <- c(50, 100, 500)
cases <- purrr::cross(.l = list("n" = n, "dgp" = dgp, "rep" = 1:reps))
probs <- c(0.95, 0.99, 0.999)

set.seed(1)
data <- purrr::map(.x = cases, \(case) {
  x <- runif(case$n)
  y <- fx[[case$dgp]](x)
  q_true <- sapply(probs, \(prob) qf[[case$dgp]](x, prob))
  # q_true <- sapply(probs, \(prob) get_q_true(x, prob, true_ccdf = true_ccdf[[case$dgp]]))
  return(list(
    "x" = x,
    "y" = y,
    "q_true" = q_true
  ))
})


plan(multisession(workers = 16))

out <- future_map(.x = data,
                  .f = function(d){
                    x <- d[["x"]]
                    y <- d[["y"]]
                    q_true <- d[["q_true"]]

                    # quantile regression
                    out_qr <- quantreg::rq(formula = y ~ x, data = data.frame(x, y),
                                           tau = probs)
                    q_pred_qr <- predict(out_qr)

                    # quantile regression coefficient modeling
                    s <- matrix(data = 1, nrow = 2, ncol = 4)
                    s[2, 2:4] <- 0
                    out_qrcm <- qrcm::iqr(formula = y ~ x, data = data.frame(x, y), s = s)
                    q_pred_qrcm <- predict(out_qrcm, type = "QF", se = FALSE, p = probs)

                    # GEV regression
                    out_lm <- lm(formula = y ~ x, data = data.frame(x, y))
                    e <- residuals(out_lm)
                    out_fgev <- evd::fgev(e, std.err = FALSE)$estimate
                    q_pred_gev <- cbind(1, x) %*% out_lm[[1]] %*% rbind(rep(1, length(probs))) +
                      rep(1, length(x)) %*% rbind(qgev(probs, out_fgev[1], out_fgev[2], out_fgev[3]))

                    # EGP regression
                    q_pred_egp <- tryCatch(expr = {
                      out_egp <- mev::fit.egp(xdat = e, model = "pt-beta")
                      cbind(1, x) %*% out_lm[[1]] %*% rbind(rep(1, length(probs))) +
                        rep(1, length(x)) %*% rbind(qegp(p = probs, scale = out_egp$param["scale"], shape = out_egp$param["shape"], kappa = out_egp$param["kappa"]))
                      },
                      error = \(e) NA)

                    # {evgam}
                    out_evgam_gev <- evgam(list(y ~ x, ~ 1, ~ 1), data.frame(x, y), family = "gev")
                    q_pred_evgam_gev <- predict(out_evgam_gev, type = "quantile", prob = probs)

                    out_evgam_gpd <- evgam(list(y ~ x, ~ 1), data.frame(x, y), family = "gpd")
                    q_pred_evgam_gpd <- predict(out_evgam_gpd, type = "quantile", prob = probs)

                    # mixQ
                    # out_mixq <- fit_mixq_ls(y = e)
                    # q_pred_mixq <-  cbind(1, x) %*% out_lm[[1]] %*% rbind(rep(1, length(probs))) +
                    #   rep(1, length(x)) %*% rbind(qmix(probs, loc = out_mixq$loc, scale = out_mixq$scale, w = out_mixq$w, xi1 = out_mixq$xi1, xi2 = out_mixq$xi2, mu_lower = out_mixq$mu_lower, mu_upper = out_mixq$mu_upper))

                    out_mixq <- fit_mixq_ls_X(y = y, X = x, tol = 1e-05)
                    q_pred_mixq <- sapply(probs, \(p) mixq_pred(u = p, X = x, params = out_mixq))

                    rmse <- as.data.frame(rbind(
                      sqrt(colMeans((q_pred_qr - q_true) ^ 2)),
                      sqrt(colMeans((q_pred_qrcm - q_true) ^ 2)),
                      sqrt(colMeans((q_pred_gev - q_true) ^ 2)),
                      sqrt(colMeans((q_pred_egp - q_true) ^ 2)),
                      sqrt(colMeans((q_pred_evgam_gev - q_true) ^ 2)),
                      sqrt(colMeans((q_pred_evgam_gpd - q_true) ^ 2)),
                      sqrt(colMeans((q_pred_mixq - q_true) ^ 2))
                      ))
                    colnames(rmse) <- paste0("p_", probs)
                    rmse$method <- c("qr", "qrcm", "gev", "egp", "evgam_gev", "evgam_gpd", "mixq")

                    return(rmse)
                  },
                  .options = furrr_options(seed = 1))

save(list = c("cases", "data", "out"),
     file = "out-sim-mixq.RData")

