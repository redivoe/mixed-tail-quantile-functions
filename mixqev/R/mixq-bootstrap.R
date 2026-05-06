
mixq_bootstrap <- function(out_fit, y, X, which_beta, B = 200){
  n <- length(y)
  
  out_fit_H0 <- fit_mixq_ls_X(y, X[, -which_beta, drop = FALSE], which_type = out_fit$which_type, verbose = FALSE)
  t_obs <- tail(out_fit_H0$obj, 1) - tail(out_fit$obj, 1)
  
  t_sim <- numeric(length = B)

  for(b in 1:B){
    y_sim <- mixq_sim(n = n, X = X[, -which_beta, drop = FALSE], params = out_fit_H0)
    out_fit_H0_sim <- fit_mixq_ls_X(y_sim, X[, -which_beta, drop = FALSE], which_type = out_fit$which_type, verbose = FALSE)
    out_fit_H1_sim <- fit_mixq_ls_X(y_sim, X, which_type = out_fit$which_type, verbose = FALSE)
    t_sim[b] <- tail(out_fit_H0_sim$obj, 1) - tail(out_fit_H1_sim$obj, 1)
  }
  
  return(list("t_obs" = t_obs,
              "t_sim" = t_sim,
              "pvalue" = mean((t_sim > abs(t_obs)) | (t_sim < -abs(t_obs)))))
}
