#' Fit the non-linear errors in variables model.
#' 
#' @param x A vector of length `N_x` of repeated observations on `n` subjects with the reference method. 
#'          Note that `N_x = r_x1 + r_x2 + ... + r_xn` and the replicate measurements are ordered by subject 
#'          in tuples of sizes `r_x1`, `r_x2`, ..., `r_xn`.
#' @param r_x A vector of length `n` with elements `[r_x1, r_x2,..., r_xn]`, the numbers of replicate
#'            measurements made by the reference method on each subject.
#' @param y A vector of length `N_y` of repeated observations on `n` subjects with the comparator method. 
#'          Note that `N_y = r_y1 + r_y2 + ... + r_yn` and the replicate measurements are ordered by subject 
#'          in tuples of sizes `r_y1`, `r_y2`, ..., `r_yn`.
#' @param r_y A vector of length `n` with elements `[r_y1, r_y2,..., r_yn]`, the numbers of replicate
#'            measurements made by the comparator method on each subject.
#' @param orders A vector of length 3 with elements `[p, d_x, d_y]` representing the orders of the bias and 
#'               precision polynomials to be fit. If order_select = TRUE, these entries are taken to be maximum 
#'               orders and all polynomials up to these orders will be fit.
#' @param order_select A boolean indicating whether a single specific model should be fit (FALSE) or whether all 
#'                     polynomials up to the orders specified in `orders` are to be fit (TRUE). 
#'                     
#' @return A list with elements: 
#'            `Parameters` the names of the model's parameters.
#'            `Estimate` the model's parameters' estimates.
#'            `BLA` the best linear approximation for each subject's true underlying trait value.
#'            `Bias_BIC` the BIC value for each bias polynomial considered, when `order_select = TRUE`.
#'            `Prec_x_BIC` the BIC value for each reference method precision polynomial considered, when `order_select = TRUE`.
#'            `Prec_y_BIC` the BIC value for each comparator method precision polynomial considered, when `order_select = TRUE`.
#'            `Optimal_orders` the BIC-selected polynomial orders, when `order_select = TRUE`. Ordered as `[p, d_x, d_y]`.
#'            
fit_model <- function(x, r_x, y, r_y, orders = c(1,0,0), order_select = FALSE){
  # Calculate mean of X and Y
  mean_x <- mean(x)
  mean_y <- mean(y)
  
  # Calculate SSW for X and Y
  low_x <- (cumsum(r_x)-r_x)+1
  upp_x <- cumsum(r_x)
  low_y <- (cumsum(r_y)-r_y)+1
  upp_y <- cumsum(r_y)
  n <- length(r_x)
  subject_means_x <- rep(0, n)
  subject_means_y <- rep(0, n)
  subject_vars_x <- rep(0, n)
  subject_vars_y <- rep(0, n)
  for(i in 1:n){
    subject_means_x[i] <- mean(x[low_x[i]:upp_x[i]])
    subject_means_y[i] <- mean(y[low_y[i]:upp_y[i]])
    subject_vars_x[i] <- var(x[low_x[i]:upp_x[i]])
    subject_vars_y[i] <- var(y[low_y[i]:upp_y[i]])
  }
  subject_means_x_rep <- rep(x = subject_means_x, times = r_x) 
  subject_means_y_rep <- rep(x = subject_means_y, times = r_y) 
  SSW_x <- sum((x - subject_means_x_rep)^2)
  SSW_y <- sum((y - subject_means_y_rep)^2)
  
  # Calculate SSB for X and Y
  SSB_x <- sum((subject_means_x_rep - mean_x)^2)
  SSB_y <- sum((subject_means_y_rep - mean_y)^2)
  
  # Estimate mu, sigma2_s
  mu_hat <- mean_x
  N_x <- sum(r_x)
  sigma2_hat <- (SSB_x - sum((1-r_x/N_x) * subject_vars_x)) / (N_x - sum(r_x^2)/N_x)
  
  # Calculate BLA of S
  s_hat <- mu_hat + sigma2_hat * (subject_means_x - mu_hat) / (sigma2_hat + subject_vars_x/r_x)
  
  # Estimate precision parameters for reference method X
  subject_sds_x <- sqrt(subject_vars_x)
  if(order_select){
    bic_x <- rep(0, orders[2]+1)
    m_x <- lm(subject_sds_x ~ 1)
    bic_x[1] <- n*log(sum(resid(m_x)^2)/n) + log(n)
    for(d in 1:orders[2]){
      m_x <- lm(subject_sds_x ~ poly(s_hat, degree = d, raw = TRUE))
      bic_x[d+1] <- n*log(sum(resid(m_x)^2)/n) + (d+1)*log(n)
    }
    bic_x_star <- matrix(data = 0, nrow = 1000, ncol = orders[2]+1)
    indx_star <- sapply(1:1000, FUN = function(b){sample(n, n, replace = TRUE)})
    for(b in 1:1000){
      m_x <- lm(subject_sds_x[indx_star[, b]] ~ 1)
      bic_x_star[b, 1] <- n*log(sum(resid(m_x)^2)/n) + log(n)
      for(d in 1:orders[2]){
        m_x <- lm(subject_sds_x[indx_star[, b]] ~ poly(s_hat[indx_star[,b]], degree = d, raw = TRUE))
        bic_x_star[b, d+1] <- n*log(sum(resid(m_x)^2)/n) + (d+1)*log(n)
      }
    }
    bic_x_hat <- 2*bic_x - colMeans(bic_x_star)
    d_x <- which.min(bic_x_hat) - 1
  }else{
    d_x <- orders[2]
  }
  if(d_x == 0){
    # Estimate constant precision and re-estimate the BLA
    sigma2_x_hat <- SSW_x / (N_x - n)
    sigma2_hat <- (SSB_x - sum((1-r_x/N_x) * sigma2_x_hat)) / (N_x - sum(r_x^2)/N_x)
    s_hat <- mu_hat + sigma2_hat * (subject_means_x - mu_hat) / (sigma2_hat + sigma2_x_hat/r_x)
    omega_x_hat <- sqrt(sigma2_x_hat)
  }else{
    # Fit precision model
    prec_x_model <- lm(subject_sds_x ~ poly(s_hat, degree = d_x, raw = TRUE))
    omega_x_hat <- as.numeric(coef(prec_x_model))
  }
  
  # Estimate precision parameters for comparator method Y
  subject_sds_y <- sqrt(subject_vars_y)
  if(order_select){
    bic_y <- rep(0, orders[3]+1)
    m_y <- lm(subject_sds_y ~ 1)
    bic_y[1] <- n*log(sum(resid(m_y)^2)/n) + log(n)
    for(d in 1:orders[3]){
      m_y <- lm(subject_sds_y ~ poly(s_hat, degree = d, raw = TRUE))
      bic_y[d+1] <- n*log(sum(resid(m_y)^2)/n) + (d+1)*log(n)
    }
    bic_y_star <- matrix(data = 0, nrow = 1000, ncol = orders[3]+1)
    for(b in 1:1000){
      m_y <- lm(subject_sds_y[indx_star[, b]] ~ 1)
      bic_y_star[b, 1] <- n*log(sum(resid(m_y)^2)/n) + log(n)
      for(d in 1:orders[3]){
        m_y <- lm(subject_sds_y[indx_star[, b]] ~ poly(s_hat[indx_star[,b]], degree = d, raw = TRUE))
        bic_y_star[b, d+1] <- n*log(sum(resid(m_y)^2)/n) + (d+1)*log(n)
      }
    }
    bic_y_hat <- 2*bic_y - colMeans(bic_y_star)
    d_y <- which.min(bic_y_hat) - 1
  }else{
    d_y <- orders[3]
  }
  if(d_y == 0){
    # Estimate constant precision 
    sigma2_y_hat <- SSW_y / (N_y - n)
    omega_y_hat <- sqrt(sigma2_y_hat)
  }else{
    # Fit precision model
    prec_y_model <- lm(subject_sds_y ~ poly(s_hat, degree = d_y, raw = TRUE))
    omega_y_hat <- as.numeric(coef(prec_y_model))
  }
  
  # Estimate bias parameters
  s_hat_y <- rep(x = s_hat, times = r_y)
  if(d_y == 0){
    var_y_hat <- rep(sigma2_y_hat, N_y)
  }else{
    var_y_hat <- (t(sapply(X = s_hat_y, FUN = function(x){x^(0:d_y)})) %*% omega_y_hat)^2
  }
  if(order_select){
    bic <- rep(0, orders[1])
    for(p in 1:orders[1]){
      bias_model <- lm(y ~ poly(s_hat_y, degree = p, raw = TRUE), weights = 1/var_y_hat)
      bic[p] <- n*log(sum(resid(bias_model)^2)/n) + (p+1)*log(n)
    }
    bic_star <- matrix(data = 0, nrow = 1000, ncol = orders[1])
    id_y <- rep(1:n, times = r_y)
    for(b in 1:1000){
      id_boot <- sample(n, size = n, replace = TRUE)
      y_star <- as.vector(unlist(sapply(id_boot, FUN = function(u){y[which(id_y == u)]})))
      s_hat_y_star <- as.vector(unlist(sapply(id_boot, FUN = function(u){s_hat_y[which(id_y == u)]})))
      for(p in 1:orders[1]){
        if(d_y == 0){
          Ny_star <- length(y_star)
          var_y_hat_star <- rep(sigma2_y_hat, Ny_star)
        }else{
          var_y_hat_star <- (t(sapply(X = s_hat_y_star, FUN = function(x){x^(0:d_y)})) %*% omega_y_hat)^2
        }
        bias_model <- lm(y_star ~ poly(s_hat_y_star, degree = p, raw = TRUE), weights = 1/var_y_hat_star)
        bic_star[b, p] <- n*log(sum(resid(bias_model)^2)/n) + (p+1)*log(n)
      }
    }
    bic_hat <- 2*bic - colMeans(bic_star)
    p <- which.min(bic_hat)
  }else{
    p <- orders[1]
  }
  bias_model <- lm(y ~ poly(s_hat_y, degree = p, raw = TRUE), weights = 1/var_y_hat)
  beta_hat <- as.numeric(coef(bias_model))
  
  # Organize all parameter estimates
  param <- c("mu", "sigma", paste0("beta",0:p), paste0("omegax", 0:d_x), paste0("omegay", 0:d_y))
  theta_hat <- c(mu_hat, sqrt(sigma2_hat), beta_hat, omega_x_hat, omega_y_hat)
  
  # Return relevant quantities
  if(!order_select){
    return(list("Parameters" = param,
                "Estimates" = theta_hat, 
                "BLA" = s_hat))  
  }else{
    return(list("Parameters" = param,
                "Estimates" = theta_hat, 
                "BLA" = s_hat, 
                "Bias_BIC" = bic_hat,
                "Prec_x_BIC" = bic_x_hat, 
                "Prec_y_BIC" = bic_y_hat,
                "Optimal_orders" = c(p, d_x, d_y)))
  }
}


#' Take bootstrap samples of the data and fit the model to each.
#' 
#' @param B A scalar indicating how many bootstrap samples to draw
#' @param x A vector of length `N_x` of repeated observations on `n` subjects with the reference method. 
#'          Note that `N_x = r_x1 + r_x2 + ... + r_xn` and the replicate measurements are ordered by subject 
#'          in tuples of sizes `r_x1`, `r_x2`, ..., `r_xn`.
#' @param r_x A vector of length `n` with elements `[r_x1, r_x2,..., r_xn]`, the numbers of replicate
#'            measurements made by the reference method on each subject.
#' @param y A vector of length `N_y` of repeated observations on `n` subjects with the comparator method. 
#'          Note that `N_y = r_y1 + r_y2 + ... + r_yn` and the replicate measurements are ordered by subject 
#'          in tuples of sizes `r_y1`, `r_y2`, ..., `r_yn`.
#' @param r_y A vector of length `n` with elements `[r_y1, r_y2,..., r_yn]`, the numbers of replicate
#'            measurements made by the comparator method on each subject.
#' @param orders A vector of length 3 with elements `[p, d_x, d_y]` representing the orders of the bias and 
#'               precision polynomials to be fit. If order_select = TRUE, these entries are taken to be maximum 
#'               orders and all polynomials up to these orders will be fit.
#'                     
#' @return A matrix with `B` rows and one column for every parameter. Row `b` provides the parameter estimates 
#'         for bootstrap sample `b`. Columns are organized as (mu, sigma, beta's, omega_x's, omega_y's).
#'
do_bootstrap <- function(B, x, r_x, y, r_y, orders){
  n <- length(r_x)
  id_x <- rep(1:n, times = r_x)
  id_y <- rep(1:n, times = r_y)
  res_boot <- foreach(b=1:B, .combine=rbind) %dopar% {
    id_boot <- sample(n, size = n, replace = TRUE)
    r_x_star <- r_x[id_boot]
    x_star <- as.vector(unlist(sapply(id_boot, FUN = function(u){x[which(id_x == u)]})))
    r_y_star <- r_y[id_boot]
    y_star <- as.vector(unlist(sapply(id_boot, FUN = function(u){y[which(id_y == u)]})))
    fit_model(x = x_star, r_x = r_x_star, y = y_star, r_y = r_y_star,
              orders = orders, order_select = FALSE)$Estimates
  }
  return(res_boot)
}


#' Calculate bootstrap-based confidence intervals for the individual model parameters.
#' 
#' @param boot A matrix with `B` rows and one column for every parameter. This could be the output of `do_bootstrap()`.
#' @param alpha A scalar corresponding to 1 minus the desired confidence level.
#' @param theta_hat A vector of parameter estimates, ordered as (mu, sigma, beta's, omega_x's, omega_y's).
#' @param type A string indicating the type of bootstrap confidence intervals to calculate. 
#'               Either `"standard"` or `"percentile"`.
#' 
#' @return A data frame with two columns containing lower and upper bounds of (1-alpha)x100% bootstrap CIs for each of the 
#'         model parameters. Rows are organized by as (mu, sigma, beta's, omega_x's, omega_y's).
#'         
calc_param_CIs <- function(boot, alpha, theta_hat, type = "standard"){
if(type == "standard"){
    st_error <- apply(X = boot, MARGIN = 2, FUN = sd)
    low <- theta_hat - qnorm(1-alpha/2)*st_error
    upp <- theta_hat + qnorm(1-alpha/2)*st_error
  }else if(type == "percentile"){
    low <- apply(X = res_boot, MARGIN = 2, FUN = quantile, prob = alpha/2)
    upp <- apply(X = res_boot, MARGIN = 2, FUN = quantile, prob = 1-alpha/2)
  }
  parameter_ci <- data.frame(LB = low, UB = upp)
  return(parameter_ci)
}


#' Calculate, return, and possibly visualize, the bias polynomial and associated confidence interval/band.
#' 
#' @param theta A vector of parameter values (could be true, estimates, or bootstrapped) for the model, 
#'             ordered (mu, sigma, beta's, omega_x's, omega_y's).
#' @param orders A vector of length 3 with elements `[p, d_x, d_y]` representing the orders of the bias and 
#'               precision polynomials.
#' @param s A vector of the underlying trait values over which to calculate bias.
#' @param interval_info A list containing elements which control the confidence interval/band. If `NULL`, no interval/band is calculated.
#'                      `type` A string indicating the type of bootstrap confidence intervals to calculate. Either "standard" or "percentile".
#'                      `boot` A matrix with `B` rows and one column for every parameter. This could be the output of `do_bootstrap()`.
#'                      `alpha` A scalar corresponding to 1 minus the desired confidence level.  
#' @param plot_info A list containing elements which control the plot constructed. If `NULL`, no plot is created.
#'                  `leg_location` A string which specifies the location of the legend.  
#'                  `colour` A string which specifies the `col` of the graphic.
#'                  `ylimits` A vector of length 2 which specifies the `ylim` values for the plot.
#'                  `s_name` A string which specifies the name of the latent trait being measured.
#'                  
#' @return A data frame which includes bias evaluated  across the `s` values as well as corresponding (1-alpha)x100% CIs and CBs, 
#'         if `interval_info` is not `NULL`.
#' @return A plot of bias (and potentially the confidence intervals/band), if `plot_info` is not `NULL`. 
#' 
calculate_bias <- function(theta, orders, s, interval_info, plot_info){
  # Extract relevant parameters
  p <- orders[1]
  beta_indx <- 3:(3+p)
  beta <- theta[beta_indx]

  # Bias estimate
  bias <- (t(sapply(X = s, FUN = function(x){x^(0:p)})) %*% beta) - s
  
  # Interval estimates
  if(!is.null(interval_info)){
    if(interval_info$type == "standard"){
      # Margins of error
      V_bias <- cov(interval_info$boot[, beta_indx])
      partial_d_bias <- t(sapply(X = s, FUN = function(x){x^(0:p)}))
      h_bias_ci <- rep(0, length(s))
      h_bias_cb <- rep(0, length(s))
      for(i in 1:length(s)){
        # Bias ME
        h_bias_ci[i] <- sqrt(qchisq(1-interval_info$alpha, df = 1) * matrix(partial_d_bias[i,], nrow = 1) %*% V_bias %*% matrix(partial_d_bias[i,], ncol = 1))
        h_bias_cb[i] <- sqrt(qchisq(1-interval_info$alpha, df = p+1) * matrix(partial_d_bias[i,], nrow = 1) %*% V_bias %*% matrix(partial_d_bias[i,], ncol = 1))
      }
      # Confidence Intervals
      bias_lo_ci <- bias - h_bias_ci
      bias_up_ci <- bias + h_bias_ci
      # Confidence Bands
      bias_lo_cb <- bias - h_bias_cb
      bias_up_cb <- bias + h_bias_cb
    }else if(interval_info$type == "percentile"){
      B <- nrow(interval_info$boot)
      bias_star <- matrix(0, nrow = B, ncol = length(s))
      for(b in 1:B){
        # Bootstrap estimates
        beta_hat_star <- interval_info$boot[b, beta_indx]
        # Bootstrapped bias calculations
        bias_star[b,] <- (t(sapply(X = s, FUN = function(x){x^(0:p)})) %*% beta_hat_star) - s
      }
      # Confidence Intervals
      bias_lo_ci <- apply(X = bias_star, MARGIN = 2, FUN = quantile, p = interval_info$alpha/2)
      bias_up_ci <- apply(X = bias_star, MARGIN = 2, FUN = quantile, p = 1-interval_info$alpha/2)
      # Confidence Bands
      res <- MBD(x = bias_star, plotting = FALSE)
      middle_C <- which(res$MBD > quantile(res$MBD, p = interval_info$alpha))
      bias_lo_cb <- apply(X = bias_star[middle_C,], MARGIN = 2, FUN = min)
      bias_up_cb <- apply(X = bias_star[middle_C,], MARGIN = 2, FUN = max)
    }
  }
  # Return results -- table
  if(is.null(interval_info)){
    result <- data.frame(Bias = bias)  
  }else{
    result <- data.frame(Bias = bias, 
                         CI_L = bias_lo_ci, 
                         CI_U = bias_up_ci, 
                         CB_L = bias_lo_cb,
                         CB_U = bias_up_cb)
  }
  # Return results -- plot
  if(!is.null(plot_info)){
    if(!is.null(interval_info)){
      plot_agreement(metric = "Bias", s = s, estimate = bias, interval = result, ylimits = plot_info$ylimits, colour = plot_info$colour, leg_location = plot_info$leg_loc, alpha = interval_info$alpha, s_name = plot_info$s_name)
    }else{
      plot_agreement(metric = "Bias", s = s, estimate = bias, interval = NULL, ylimits = plot_info$ylimits, colour = plot_info$colour, s_name = plot_info$s_name)
    }
  }
  return(result)
}


#' Calculate, return, and possibly visualize, the precision polynomial and associated confidence interval/band.
#' 
#' @param theta A vector of parameter values (could be true, estimates, or bootstrapped) for the model, 
#'             ordered (mu, sigma, beta's, omega_x's, omega_y's).
#' @param orders A vector of length 3 with elements `[p, d_x, d_y]` representing the orders of the bias and 
#'               precision polynomials.
#' @param method A string indicating whether precision is calculated for for the reference method ("R"), 
#'               or the comparator method ("C").                
#' @param s A vector of the underlying trait values over which to calculate precision.
#' @param interval_info A list containing elements which control the confidence interval/band. If `NULL`, no interval/band is calculated.
#'                      `type` A string indicating the type of bootstrap confidence intervals to calculate. Either `"standard"` or `"percentile"`.
#'                      `boot` A matrix with `B` rows and one column for every parameter. This could be the output of `do_bootstrap()`.
#'                      `alpha` A scalar corresponding to 1 minus the desired confidence level.  
#' @param plot_info A list containing elements which control the plot constructed. If `NULL`, no plot is created.
#'                  `leg_location` A string which specifies the location of the legend.  
#'                  `colour` A string which specifies the `col` of the graphic.
#'                  `ylimits` A vector of length 2 which specifies the `ylim` values for the plot.
#'                  `s_name` A string which specifies the name of the latent trait being measured.
#'                  `method_name` A string which specifies whether the method is the `reference` or the `comparator`.
#'                  
#' @return A data frame which includes precision evaluated  across the `s` values as well as corresponding (1-alpha)x100% CIs and CBs, if `interval_info` is not `NULL`.
#' @return A plot of precision (and potentially the confidence intervals/band), if `plot_info` is not `NULL`.
#' 
calculate_precision <- function(theta, orders, method, s, interval_info, plot_info){
  # Extract relevant parameters
  p <- orders[1]
  if(method == "R"){
    omega_indx <- (4+p):(4+p+orders[2])  
    d <- orders[2]
  }else if(method == "C"){
    omega_indx <- (5+p+orders[2]):length(theta)  
    d <- orders[3]
  }
  omega <- theta[omega_indx]
  # Precision estimate
  if(d == 0){
    precision <- rep(omega, length(s)) 
  }else{
    precision <- (t(sapply(X = s, FUN = function(x){x^(0:d)})) %*% omega)    
  }
  # Interval estimates
  if(!is.null(interval_info)){
    if(interval_info$type == "standard"){
      # Margins of error
      if(d == 0){
        V_prec <- var(interval_info$boot[, omega_indx])
      }else{
        V_prec <- cov(interval_info$boot[, omega_indx])
      }
      if(d == 0){
        partial_d_prec <- matrix(1, nrow = length(s), ncol = 1)
      }else{
        partial_d_prec <- t(sapply(X = s, FUN = function(x){x^(0:d)}))  
      }
      h_prec_ci <- rep(0, length(s))
      h_prec_cb <- rep(0, length(s))
      for(i in 1:length(s)){
        # Precision ME
        h_prec_ci[i] <- sqrt(qchisq(1-interval_info$alpha, df = 1) * matrix(partial_d_prec[i,], nrow = 1) %*% V_prec %*% matrix(partial_d_prec[i,], ncol = 1))
        h_prec_cb[i] <- sqrt(qchisq(1-interval_info$alpha, df = d+1) * matrix(partial_d_prec[i,], nrow = 1) %*% V_prec %*% matrix(partial_d_prec[i,], ncol = 1))
      }
      # Confidence Intervals
      prec_lo_ci <- precision - h_prec_ci
      prec_up_ci <- precision + h_prec_ci
      # Confidence Bands
      prec_lo_cb <- precision - h_prec_cb
      prec_up_cb <- precision + h_prec_cb
    }else if(interval_info$type == "percentile"){
      B <- nrow(interval_info$boot)
      precision_star <- matrix(0, nrow = B, ncol = length(s))
      for(b in 1:B){
        # Bootstrap estimates
        omega_hat_star <- interval_info$boot[b, omega_indx]
        # Bootstrapped precision calculations
        if(d == 0){
          precision_star[b,] <- rep(omega_hat_star, length(s))  
        }else{
          precision_star[b,] <- (t(sapply(X = s, FUN = function(x){x^(0:d)})) %*% omega_hat_star)
        }
      }
      # Confidence Intervals
      prec_lo_ci <- apply(X = precision_star, MARGIN = 2, FUN = quantile, p = interval_info$alpha/2)
      prec_up_ci <- apply(X = precision_star, MARGIN = 2, FUN = quantile, p = 1-interval_info$alpha/2)
      # Confidence Bands
      res <- MBD(x = precision_star, plotting = FALSE)
      middle_C <- which(res$MBD > quantile(res$MBD, p = interval_info$alpha))
      prec_lo_cb <- apply(X = precision_star[middle_C,], MARGIN = 2, FUN = min)
      prec_up_cb <- apply(X = precision_star[middle_C,], MARGIN = 2, FUN = max)
    }
  }
  # Return results -- table
  if(is.null(interval_info)){
    result <- data.frame(Precision = precision)  
  }else{
    result <- data.frame(Precision = precision, 
                         CI_L = prec_lo_ci,
                         CI_U = prec_up_ci,
                         CB_L = prec_lo_cb,
                         CB_U = prec_up_cb)
  }
  # Return results -- plot
  if(!is.null(plot_info)){
    if(!is.null(interval_info)){
      plot_agreement(metric = "Precision", s = s, estimate = precision, interval = result, ylimits = plot_info$ylimits, colour = plot_info$colour, leg_location = plot_info$leg_loc, alpha = plot_info$alpha, s_name = plot_info$s_name)
    }else{
      plot_agreement(metric = "Precision", s = s, estimate = precision, interval = NULL, ylimits = plot_info$ylimits, colour = plot_info$colour, s_name = plot_info$s_name)
    }
    if(!is.null(plot_info$method_name)){
      title(ylab = paste0("Precision (", plot_info$method_name, ")"))  
    }else{
      if(method == "R"){
        title(ylab = "Precision (reference)")  
      }else{
        title(ylab = "Precision (comparator)")
      }      
    }
  }
  return(result)
}


#' Calculate, return, and possibly visualize, the probability of agreement and associated confidence interval/band.
#' 
#' @param theta A vector of parameter values (could be true, estimates, or bootstrapped) for the model, 
#'             ordered (mu, sigma, beta's, omega_x's, omega_y's).
#' @param delta An input that specifies the indifference region. Can be a scalar, a string, or a matrix. The scalar should be `x` 
#'              such that the indifference region is defined as `[-x,x]`. The string should be used to specify a percent tolerance 
#'              whereby `delta = "x%"`. The matrix should have 2 columns containing the lower bound (first column) and upper bound 
#'              (second column) of the indifference region. The matrix should have many rows as elements of `s`. 
#' @param orders A vector of length 3 with elements `[p, d_x, d_y]` representing the orders of the bias and 
#'               precision polynomials.
#' @param s A vector of the underlying trait values over which to calculate the PoA.
#' @param interval_info A list containing elements which control the confidence interval/band. If `NULL`, no interval/band is calculated.
#'                      `type` A string indicating the type of bootstrap confidence intervals to calculate. Either `"standard"` or `"percentile"`.
#'                      `boot` A matrix with `B` rows and one column for every parameter. This could be the output of `do_bootstrap()`.
#'                      `alpha` A scalar corresponding to 1 minus the desired confidence level.  
#' @param plot_info A list containing elements which control the plot constructed. If `NULL`, no plot is created.
#'                  `leg_location` A string which specifies the location of the legend.  
#'                  `colour` A string which specifies the `col` of the graphic.
#'                  `ylimits` A vector of length 2 which specifies the `ylim` values for the plot.
#'                  `s_name` A string which specifies the name of the latent trait being measured.
#'                  
#' @return A data frame which includes PoA evaluated  across the `s` values as well as corresponding (1-alpha)x100% CIs and CBs, if `interval_info` is not `NULL`.
#' @return A plot of the PoA (and potentially the confidence intervals/band), if `plot_info` is not `NULL`.
#' 
calculate_poa <- function(theta, delta, orders, s, interval_info, plot_info){
  # Extract relevant parameters
  p <- orders[1]
  d_x <- orders[2]
  d_y <- orders[3]
  beta_indx <- 3:(3+p)
  beta <- theta[beta_indx]
  omega_x_indx <- (4+p):(4+p+d_x)
  omega_x <- theta[omega_x_indx]
  omega_y_indx <- (5+p+d_x):length(theta)
  omega_y <- theta[omega_y_indx]
  # Bias
  bias <- (t(sapply(X = s, FUN = function(x){x^(0:p)})) %*% beta) - s
  # Reference precision
  if(d_x == 0){
    precision_x <- rep(omega_x, length(s)) 
  }else{
    precision_x <- (t(sapply(X = s, FUN = function(x){x^(0:d_x)})) %*% omega_x)    
  }
  # Comparator precision
  if(d_y == 0){
    precision_y <- rep(omega_y, length(s)) 
  }else{
    precision_y <- (t(sapply(X = s, FUN = function(x){x^(0:d_y)})) %*% omega_y)    
  }
  # Probability of agreement estimate
  if(is.character(delta)){
    delta <- as.numeric(strsplit(x = delta, split = "%"))/100
    delta_mat <- matrix(0, nrow = length(s), ncol = 2)
    delta_mat[,1] <- -delta*s
    delta_mat[,2] <- delta*s
  }else if(is.matrix(delta)){
    delta_mat <- delta
  }else{
    delta_mat <- matrix(c(-delta,delta), ncol = 2)
  }
  deltaL <- delta_mat[,1]
  deltaU <- delta_mat[,2]
  zU <- (deltaU - bias) / sqrt(precision_x^2 + precision_y^2)
  zL <- (deltaL - bias) / sqrt(precision_x^2 + precision_y^2)
  poa <- pnorm(zU) - pnorm(zL)
  
  # Interval estimates
  if(!is.null(interval_info)){
    if(interval_info$type == "standard"){
      # Margin of error
      V_poa <- cov(interval_info$boot[, c(beta_indx, omega_x_indx, omega_y_indx)])
      h_poa_ci <- rep(0, length(s))
      h_poa_cb <- rep(0, length(s))
      for(i in 1:length(s)){
        # PoA ME
        partial_d_beta <- (-dnorm(zU[i]) + dnorm(zL[i])) * (s[i]^(0:p)) / sqrt(precision_x[i]^2 + precision_y[i]^2)
        if(nrow(delta_mat) == 1){
          partial_d_omega_x <- (-dnorm(zU[i]) * (deltaU - bias[i]) + dnorm(zL[i]) * (deltaL - bias[i])) * precision_x[i] * (s[i]^(0:d_x)) / (precision_x[i]^2 + precision_y[i]^2)^(3/2)
          partial_d_omega_y <- (-dnorm(zU[i]) * (deltaU - bias[i]) + dnorm(zL[i]) * (deltaL - bias[i])) * precision_y[i] * (s[i]^(0:d_y)) / (precision_x[i]^2 + precision_y[i]^2)^(3/2)      
        }else{
          partial_d_omega_x <- (-dnorm(zU[i]) * (deltaU[i] - bias[i]) + dnorm(zL[i]) * (deltaL[i] - bias[i])) * precision_x[i] * (s[i]^(0:d_x)) / (precision_x[i]^2 + precision_y[i]^2)^(3/2)
          partial_d_omega_y <- (-dnorm(zU[i]) * (deltaU[i] - bias[i]) + dnorm(zL[i]) * (deltaL[i] - bias[i])) * precision_y[i] * (s[i]^(0:d_y)) / (precision_x[i]^2 + precision_y[i]^2)^(3/2)      
        }
        partial_d_poa <- matrix(c(partial_d_beta, partial_d_omega_x, partial_d_omega_y)/((poa[i]-1)*log(1-poa[i])), ncol = 1)
        h_poa_ci[i] <- sqrt(qchisq(1-interval_info$alpha, df = 1) * t(partial_d_poa) %*% V_poa %*% partial_d_poa)
        h_poa_cb[i] <- sqrt(qchisq(1-interval_info$alpha, df = p + d_x + d_y + 3) * t(partial_d_poa) %*% V_poa %*% partial_d_poa)
      }
      # Confidence Intervals
      cloglog_poa_lo_ci <- log(-log(1-poa)) - h_poa_ci
      cloglog_poa_up_ci <- log(-log(1-poa)) + h_poa_ci
      poa_lo_ci <- 1-exp(-exp(cloglog_poa_lo_ci))
      poa_lo_ci[poa==1] <- 1
      poa_up_ci <- 1-exp(-exp(cloglog_poa_up_ci))  
      poa_up_ci[poa==1] <- 1  
      # Confidence Bands
      cloglog_poa_lo_cb <- log(-log(1-poa)) - h_poa_cb
      cloglog_poa_up_cb <- log(-log(1-poa)) + h_poa_cb
      poa_lo_cb <- 1-exp(-exp(cloglog_poa_lo_cb))
      poa_lo_cb[poa==1] <- 1
      poa_up_cb <- 1-exp(-exp(cloglog_poa_up_cb))  
      poa_up_cb[poa==1] <- 1  
    }
    else if(interval_info$type == "percentile"){
      B <- nrow(interval_info$boot)
      poa_star <- matrix(0, nrow = B, ncol = length(s))
      for(b in 1:B){
        # Bootstrap estimates
        beta_hat_star <- interval_info$boot[b, beta_indx]
        omega_x_hat_star <- interval_info$boot[b, omega_x_indx]
        omega_y_hat_star <- interval_info$boot[b, omega_y_indx]
        theta_hat_star <- c(beta_hat_star, omega_x_hat_star, omega_y_hat_star)
        # Bootstrapped bias calculations
        bias_star <- (t(sapply(X = s, FUN = function(x){x^(0:p)})) %*% beta_hat_star) - s
        # Bootstrapped reference precision calculations
        if(d_x == 0){
          precision_x_star <- rep(omega_x_hat_star, length(s))  
        }else{
          precision_x_star <- (t(sapply(X = s, FUN = function(x){x^(0:d_x)})) %*% omega_x_hat_star)
        }
        # Bootstrapped comparator precision calculations
        if(d_y == 0){
          precision_y_star <- rep(omega_y_hat_star, length(s))  
        }else{
          precision_y_star <- (t(sapply(X = s, FUN = function(x){x^(0:d_y)})) %*% omega_y_hat_star)  
        }
        # Bootstrapped probability of agreement calculations
        poa_star[b,] <- pnorm((deltaU - bias_star) / sqrt(precision_x_star^2 + precision_y_star^2)) - pnorm((deltaL - bias_star) / sqrt(precision_x_star^2 + precision_y_star^2)) 
      }
      # Confidence Intervals
      poa_lo_ci <- apply(X = poa_star, MARGIN = 2, FUN = quantile, p = interval_info$alpha/2)
      poa_up_ci <- apply(X = poa_star, MARGIN = 2, FUN = quantile, p = 1-interval_info$alpha/2)
      # Confidence Bands
      res <- MBD(x = poa_star, plotting = FALSE)
      middle_C <- which(res$MBD > quantile(res$MBD, p = interval_info$alpha))
      poa_lo_cb <- apply(X = poa_star[middle_C,], MARGIN = 2, FUN = min)
      poa_up_cb <- apply(X = poa_star[middle_C,], MARGIN = 2, FUN = max)
    }
  }
  # Return results -- table
  if(is.null(interval_info)){
    result <- data.frame(PoA = poa)  
  }else{
    result <- data.frame(PoA = poa, 
                         CI_L = poa_lo_ci,
                         CI_U = poa_up_ci,
                         CB_L = poa_lo_cb,
                         CB_U = poa_up_cb)
  }
  # Return results -- plot
  if(!is.null(plot_info)){
    if(!is.null(interval_info)){
      plot_agreement(metric = "PoA", s = s, estimate = poa, interval = result, ylimits = plot_info$ylimits, colour = plot_info$colour, leg_location = plot_info$leg_location, alpha = interval_info$alpha, s_name = plot_info$s_name)
    }else{
      plot_agreement(metric = "PoA", s = s, estimate = poa, interval = NULL, ylimits = plot_info$ylimits, colour = plot_info$colour, s_name = plot_info$s_name)
    }
  }
  return(result)
}


#' Helper function construct bias, precision, and PoA plots.
#' 
#' @param metric A string indicating which type of plot is being constructed ("Bias", or "Precision", or "PoA").
#' @param s A vector of the underlying trait values that define the horizontal axis of the plot.
#' @param estimate A vector with the same length as `s` containing the estimated values of the metric to be plotted.
#' @param interval A data frame with 4 columns titled CI_L, CI_U, CB_L, CB_U, containing the lower and upper bounds 
#'                 of the confidence intervals/bands. This should have as many rows as elements in `s`. If `NULL`, no
#'                 confidence intervals/bands are plotted.
#' @param ylimits A vector of length 2 which specifies the `ylim` values for the plot.
#' @param colour A string which specifies the `col` of the graphic.
#' @param leg_location A string which specifies the location of the legend.  
#' @param s_name A string which specifies the name of the latent trait being measured.
#' @param alpha A scalar corresponding to 1 minus the desired confidence level. Ignored if `interval = NULL`.  
#' 
#' @return A plot of the metric (and potentially the confidence intervals/band).
#' 
plot_agreement <- function(metric, s, estimate, interval, ylimits, colour, leg_location, s_name, alpha){
  if(!is.null(interval)){
    lo_ci <- interval$CI_L
    up_ci <- interval$CI_U
    lo_cb <- interval$CB_L
    up_cb <- interval$CB_U
    if(metric != "PoA"){
      if(is.null(ylimits)){
        ylimits <- c(min(lo_cb), max(up_cb))
      }
      if(metric == "Bias"){
        ylabel <- metric
      }else{
        ylabel = ""
      }
    }else{
      ylimits <- c(0,1)
      ylabel <- "Probability"
    }
    plot(x = s, y = estimate, type = "l", col = colour, ylim = ylimits, ylab = ylabel, main = "", xlab = "")
    polygon(x = c(s, rev(s)), y = c(lo_ci, rev(up_ci)), col = adjustcolor(colour, 0.25), border = NA)
    polygon(x = c(s, rev(s)), y = c(lo_cb, rev(up_cb)), col = adjustcolor(colour, 0.25), border = NA)
    lines(x = s, y = estimate, col = colour, lwd = 2)
    if(metric == "Bias"){abline(h = 0, lty = 2)}
    legend(leg_location, legend = c(paste0(metric, " Estimate"), paste0(100*(1-alpha), "% Conf Interval"), paste0(100*(1-alpha), "% Conf Band")), 
           col = c(colour, adjustcolor(colour, 0.5), adjustcolor(colour, 0.25)), 
           lty = 1, bty = "n", cex = 0.8, lwd = c(2,8,8))
  }else{
    if(metric != "PoA"){
      if(is.null(ylimits)){
        ylimits <- c(min(estimate), max(estimate))
      }
      ylabel <- metric
    }else{
      ylimits <- c(0,1)
      ylabel <- "Probability"
    }
    plot(x = s, y = estimate, type = "l", col = colour, main = "", ylab = ylabel, lwd = 2, ylim = ylimits, xlab = "")
    if(metric == "Bias"){abline(h = 0, lty = 2)}
  }
  if(is.null(s_name)){s_name = "s"} 
  title(xlab = s_name)
}


#' Estimation of the best linear approximation for each subject's underlying latent trait.
#' 
#' @param x A vector of length `N_x` of repeated observations on `n` subjects with the reference method. 
#'          Note that `N_x = r_x1 + r_x2 + ... + r_xn` and the replicate measurements are ordered by subject 
#'          in tuples of sizes `r_x1`, `r_x2`, ..., `r_xn`.
#' @param r_x A vector of length `n` with elements `[r_x1, r_x2,..., r_xn]`, the numbers of replicate
#'            measurements made by the reference method on each subject.
#'            
#' @return A vector of length `n` containing the best linear approximation of the underlying latent trait 
#'         for each of the `n` subjects.
#' 
calculate_bla <- function(x, r_x){
  # Calculate mean of X
  mean_x <- mean(x)
  
  # Calculate SSW for X
  low_x <- (cumsum(r_x)-r_x)+1
  upp_x <- cumsum(r_x)
  n <- length(r_x)
  subject_means_x <- rep(0, n)
  subject_vars_x <- rep(0, n)
  for(i in 1:n){
    subject_means_x[i] <- mean(x[low_x[i]:upp_x[i]])
    subject_vars_x[i] <- var(x[low_x[i]:upp_x[i]])
  }
  subject_means_x_rep <- rep(x = subject_means_x, times = r_x) 
  SSW_x <- sum((x - subject_means_x_rep)^2)
  
  # Calculate SSB for X
  SSB_x <- sum((subject_means_x_rep - mean_x)^2)
  
  # Estimate mu, sigma2_s
  mu_hat <- mean_x
  N_x <- sum(r_x)
  sigma2_hat <- (SSB_x - sum((1-r_x/N_x) * subject_vars_x)) / (N_x - sum(r_x^2)/N_x)
  
  # Calculate BLA of S
  s_hat <- mu_hat + sigma2_hat * (subject_means_x - mu_hat) / (sigma2_hat + subject_vars_x/r_x)
  return(s_hat)
}


#' Calculate, return, and possibly visualize, the conditional probability of agreement and associated confidence interval.
#'
#' @param x A vector of length `N_x` of repeated observations on `n` subjects with the reference method. 
#'          Note that `N_x = r_x1 + r_x2 + ... + r_xn` and the replicate measurements are ordered by subject 
#'          in tuples of sizes `r_x1`, `r_x2`, ..., `r_xn`.
#' @param r_x A vector of length `n` with elements `[r_x1, r_x2,..., r_xn]`, the numbers of replicate
#'            measurements made by the reference method on each subject.
#' @param y A vector of length `N_y` of repeated observations on `n` subjects with the comparator method. 
#'          Note that `N_y = r_y1 + r_y2 + ... + r_yn` and the replicate measurements are ordered by subject 
#'          in tuples of sizes `r_y1`, `r_y2`, ..., `r_yn`.
#' @param r_y A vector of length `n` with elements `[r_y1, r_y2,..., r_yn]`, the numbers of replicate
#'            measurements made by the comparator method on each subject.
#' @param orders A vector of length 3 with elements `[p, d_x, d_y]` representing the orders of the bias and 
#'               precision polynomials to be fit. If order_select = TRUE, these entries are taken to be maximum 
#'               orders and all polynomials up to these orders will be fit.
#' @param delta An input that specifies the indifference region. Can be a scalar, a string, or a matrix. The scalar should be `x` 
#'              such that the indifference region is defined as `[-x,x]`. The string should be used to specify a percent tolerance 
#'              whereby `delta = "x%"`. The matrix should have 2 columns containing the lower bound (first column) and upper bound 
#'              (second column) of the indifference region. The matrix should have many rows as elements of `s`. 
#' @param B A scalar indicating how many bootstrap samples to draw
#' @param alpha A scalar corresponding to 1 minus the desired confidence level.  
#' @param plot_info A list containing elements which control the plot constructed. If `NULL`, no plot is created.
#'                  `colour` A string which specifies the `col` of the graphic.
#'                  `s_name` A string which specifies the name of the latent trait being measured.
#'            
#' @return A list with elements:
#'            `PoA` A vector of length `n` containing the conditional PoA estimates for each subject.
#'            `LCL` A vector of length `n` containing the lower confidence limit for the conditional PoA value for each subject.
#'            `UCL` A vector of length `n` containing the upper confidence limit for the conditional PoA value for each subject.
#'            `BLA` A vector of length `n` containing the best linear approximation for each subject's true underlying trait value.
#'            `BLA_x` A vector of length `n` containing the best linear approximation for each subject's true underlying trait value, based on the reference method's measurements.
#'            `BLA_y` A vector of length `n` containing the best linear approximation for each subject's true underlying trait value, based on the calibrated comparator method's measurements.
#'            `PoA_star` A matrix with `B` rows and `n` columns containing bootstrap estimates of the conditional PoA for each subject.
#' @return A plot of the conditional PoA with associated confidence intervals, if `plot_info` is not `NULL`.
#' 
conditional_poa <- function(x, r_x, y, r_y, orders, delta, B, alpha, plot_info){
  # Extract relevant parameters
  n <- length(r_x)
  p <- orders[1]
  d_x <- orders[2]
  d_y <- orders[3]
  # Sample Estimate
  res <- fit_model(x = x, r_x = r_x, y = y, r_y = r_y, orders = orders, order_select = FALSE)
  theta_hat <- res$Estimates
  s_hat_x <- res$BLA
  beta_hat <- theta_hat[3:(3+p)]
  y_cal <- (y - beta_hat[1])/beta_hat[2]
  s_hat_y <- calculate_bla(y_cal, r_y)
  s_hat <- (r_x*s_hat_x + r_y*s_hat_y) / (r_x + r_y)
  poa <- calculate_poa(theta = theta_hat, delta = delta, orders = orders, s = s_hat, interval_info = NULL, plot_info = NULL)$PoA
  # Conditional bootstrap estimates
  id_x <- rep(1:n, times = r_x)
  id_y <- rep(1:n, times = r_y)
  poa_star <- foreach(i = 1:n, .combine = rbind) %dopar% {
    poa_star_i <- foreach(b = 1:B, .combine = c) %dopar% {
      id_boot <- sample((1:n)[-i], size = n-1, replace = TRUE)
      if(i == 1){
        id_boot <- c(i, id_boot)
      }else if(i == n){
        id_boot <- c(id_boot, i)
      }else{
        id_boot <- c(id_boot[1:(i-1)], i, id_boot[i:(n-1)])
      }
      r_x_star <- r_x[id_boot]
      r_y_star <- r_y[id_boot]
      x_star <- as.vector(unlist(sapply(id_boot, FUN = function(u){x[which(id_x == u)]})))
      y_star <- as.vector(unlist(sapply(id_boot, FUN = function(u){y[which(id_y == u)]})))
      res_star <- fit_model(x = x_star, r_x = r_x_star, y = y_star, r_y = r_y_star,
                            orders = orders, order_select = FALSE)
      theta_hat_star <- res_star$Estimates
      s_hat_x_star <- res_star$BLA
      beta_hat_star <- theta_hat_star[3:(3+p)]
      y_star_cal <- (y_star - beta_hat_star[1])/beta_hat_star[2]
      s_hat_y_star <- calculate_bla(y_star_cal, r_y_star)
      s_hat_star <- (r_x_star*s_hat_x_star + r_y_star*s_hat_y_star) / (r_x_star + r_y_star)
      calculate_poa(theta = theta_hat_star, delta = delta, orders = orders, s = s_hat_star[i], 
                    interval_info = NULL, plot_info = NULL)$PoA
    }
    poa_star_i
  }
  # Confidence intervals
  LCL <- as.numeric(apply(X = poa_star, MARGIN = 1, FUN = quantile, p = alpha/2))
  UCL <- as.numeric(apply(X = poa_star, MARGIN = 1, FUN = quantile, p = 1-alpha/2))
  if(!is.null(plot_info)){
    if(is.null(plot_info$colour)){
      colour <- "red"
    }else{
      colour <- plot_info$colour
    }
    plot(x = s_hat, y = poa, ylim = c(0,1), col = "white", pch = 16, cex = 0.75,
         ylab = "Probability", xlab = "")
    for(i in 1:n){
      segments(x0 = s_hat[i], y0 = LCL[i], x1 = s_hat[i], y1 = UCL[i], col = colour, lwd = 0.5)
    }
    points(x = s_hat, y = poa, col = colour, pch = 16, cex = 0.75)
    points(x = s_hat, y = LCL, pch = "-", , col = colour)
    points(x = s_hat, y = UCL, pch = "-", col = colour)
    if(is.null(plot_info$s_name)){
      title(xlab = expression(hat(s)))
    }else{
      title(xlab = plot_info$s_name)  
    }
  }
  return(list(PoA = poa, LCL = LCL, UCL = UCL, BLA = s_hat, BLA_x = s_hat_x, BLA_y = s_hat_y, PoA_star = poa_star))
}


#' Construct plots of the raw data to informally visualize agreement
#' 
#' @param x A vector of length `N_x` of repeated observations on `n` subjects with the reference method. 
#'          Note that `N_x = r_x1 + r_x2 + ... + r_xn` and the replicate measurements are ordered by subject 
#'          in tuples of sizes `r_x1`, `r_x2`, ..., `r_xn`.
#' @param r_x A vector of length `n` with elements `[r_x1, r_x2,..., r_xn]`, the numbers of replicate
#'            measurements made by the reference method on each subject.
#' @param y A vector of length `N_y` of repeated observations on `n` subjects with the comparator method. 
#'          Note that `N_y = r_y1 + r_y2 + ... + r_yn` and the replicate measurements are ordered by subject 
#'          in tuples of sizes `r_y1`, `r_y2`, ..., `r_yn`.
#' @param r_y A vector of length `n` with elements `[r_y1, r_y2,..., r_yn]`, the numbers of replicate
#'            measurements made by the comparator method on each subject.
#' @param s_name A string which specifies the name of the latent trait being measured.
#' @param ref_name A string which specifies the name of the reference method.
#' @param com_name A string which specifies the name of the comparator method.
#'            
#' @return A 2x2 figure containing various scatter plots and a Bland-Altman plot that visualize the raw data.
#' 
data_vizualization <- function(x, r_x, y, r_y, s_name = "Individual", ref_name = "Reference", com_name = "Comparator"){
  # Visualize agreement
  n <- length(r_x)
  par(mfrow=c(2,2))
  # Method-specific scatter plot (reference)
  subject_x <- rep(x = 1:n, times = r_x)
  x_avg <- aggregate(x, by = list(subject_x), FUN = mean)$x
  x_avgs <- rep(x_avg, times = r_x) # repeat each subject average r times
  plot(x = x_avgs, y = x, col = "black", pch = 19,
       main = paste0(s_name, " Measurements vs. Subject Average (", ref_name, ")"),
       xlab = "Subject Average",
       ylab = "Individual Measurement",
       ylim = c(min(x, y), max(x, y)),
       xlim = c(min(x, y), max(x, y)))
  abline(a=0, b=1, col = "red")
  legend("bottomright", inset = 0.02, legend = "Line of Equality", col = "red", cex = 0.8, bty = "n", lty = 1)
  
  # Method-specific scatter plot (comparator)
  subject_y <- rep(x = 1:n, times = r_y)
  y_avg <- aggregate(y, by = list(subject_y), FUN = mean)$x
  y_avgs <- rep(y_avg, times = r_y) # repeat each subject average r times
  plot(x = y_avgs, y = y, col = "black", pch = 19,
       main = paste0(s_name, " Measurements vs. Subject Average (", com_name, ")"),
       xlab = "Subject Average",
       ylab = "Individual Measurement",
       ylim = c(min(x, y), max(x, y)),
       xlim = c(min(x, y), max(x, y)))
  abline(a=0, b=1, col = "red")
  legend("bottomright", inset = 0.02, legend = "Line of Equality", col = "red", cex = 0.8, bty = "n", lty = 1)
  
  # Reference vs comparator scatter plot 
  plot(x_avg, y_avg, col="black", pch=19, main = "Scatter Plot of Subject Averages",
       xlim = c(min(c(x_avg, y_avg)), max(c(x_avg, y_avg))),
       ylim = c(min(c(x_avg, y_avg)), max(c(x_avg, y_avg))),
       xlab = paste0("Subject Average Measurement (", ref_name, ")"),
       ylab = paste0("Subject Average Measurement (", com_name, ")"))
  abline(a = 0, b = 1, lty = 1, col = "red")
  legend("bottomright", inset = 0.02, legend = "Line of Equality", col = "red", cex = 0.8, bty = "n", lty = 1)
  
  # Bland-Altman plot
  avg <- (x_avg + y_avg)/2 # averages
  diff <-  y_avg - x_avg # differences
  avgdiff <- mean(diff) # average difference
  var_x <- mean(aggregate(x, by = list(subject_x), FUN = var)$x) # repeatability for MS1 (pooled across individuals)
  var_y <- mean(aggregate(y, by = list(subject_y), FUN = var)$x) # repeatability for MS2 (pooled across individuals)
  var_d <- var(diff) + (1-(1/n)*sum(1/r_x))*var_x + (1-(1/n)*sum(1/r_y))*var_y # estimate of stdev of single measurement differences (based on variance of difference of averages)
  LLA <- avgdiff - qnorm(0.975)*sqrt(var_d)
  ULA <- avgdiff + qnorm(0.975)*sqrt(var_d)
  loa_se <- sqrt(var_d/n + ((qnorm(0.975)^2)/(2*var_d)) * ((var_d^2/(n-1)) + (((1-(1/n)*sum(1/r_x))^2)*var_x^2)/sum(r_x-1) + (((1-(1/n)*sum(1/r_y)))*var_y^2)/sum(r_y-1)))  
  LLA_up <- LLA + qnorm(0.975)*loa_se
  LLA_lo <- LLA - qnorm(0.975)*loa_se
  ULA_up <- ULA + qnorm(0.975)*loa_se
  ULA_lo <- ULA - qnorm(0.975)*loa_se
  
  plot(x = avg, y = diff, col="black", pch=19, main = "Bland and Altman Difference Plot",
       xlab = "Average Measurements: (X+Y)/2",
       ylab = "Differences: Y-X", 
       ylim = c(min(LLA_lo, min(diff)), max(ULA_up, max(diff))),
       xlim = range(avg))
  abline(h = avgdiff, col="red", lty = 1)
  
  x_range <- seq(min(avg)-0.1*diff(range(avg)), max(avg)+0.1*diff(range(avg)), length.out = 1000)
  polygon(x = c(x_range, rev(x_range)),
          y = c(rep(ULA_up, length(x_range)), rep(ULA_lo, length(x_range))),
          col = adjustcolor("red", 0.2), border = FALSE)
  abline(h=ULA, col="red", lty=2)
  polygon(x = c(x_range, rev(x_range)),
          y = c(rep(LLA_up, length(x_range)), rep(LLA_lo, length(x_range))),
          col = adjustcolor("red", 0.2), border = FALSE)
  abline(h=LLA, col="red", lty=2)
  legend("bottomright", inset = 0.02, legend = c("Average difference","Limits of agreement (LoA)", "95% CI for LoA"), 
         col = c("red", "red", adjustcolor("red", 0.3)), 
         lty = c(1,2,1), lwd = c(1,1,10), cex = 0.8, bty = "n")
}


#' Construct diagnostic plots from model residuals.
#' 
#' @param x A vector of length `N_x` of repeated observations on `n` subjects with the reference method. 
#'          Note that `N_x = r_x1 + r_x2 + ... + r_xn` and the replicate measurements are ordered by subject 
#'          in tuples of sizes `r_x1`, `r_x2`, ..., `r_xn`.
#' @param r_x A vector of length `n` with elements `[r_x1, r_x2,..., r_xn]`, the numbers of replicate
#'            measurements made by the reference method on each subject.
#' @param y A vector of length `N_y` of repeated observations on `n` subjects with the comparator method. 
#'          Note that `N_y = r_y1 + r_y2 + ... + r_yn` and the replicate measurements are ordered by subject 
#'          in tuples of sizes `r_y1`, `r_y2`, ..., `r_yn`.
#' @param r_y A vector of length `n` with elements `[r_y1, r_y2,..., r_yn]`, the numbers of replicate
#'            measurements made by the comparator method on each subject.
#' @param fitted_model A model object return by `fit_model()`.
#' @param orders A vector of length 3 with elements `[p, d_x, d_y]` representing the orders of the bias and 
#'               precision polynomials to be fit. If order_select = TRUE, these entries are taken to be maximum 
#'               orders and all polynomials up to these orders will be fit.
#'               
#' @return A 2x3 figure containing scatter plot (left), histogram (middle), and QQ-plot (right) of the 
#'         standardized residuals for the reference method (top) and comparator method (bottom).
#'               
model_diagnostics <- function(x, r_x, y, r_y, fitted_model, orders){
  # Extract relevant information
  p <- orders[1]
  d_x <- orders[2]
  d_y <- orders[3]
  beta_indx <- 3:(3+p)
  omega_x_indx <- (4+p):(4+p+d_x)
  omega_y_indx <- (5+p+d_x):(5+p+d_x+d_y)
  theta_hat <- fitted_model$Estimates
  par(mfrow=c(2,3))
  # Calculate reference method residuals
  s_hat_x <- rep(fitted_model$BLA, times = r_x)
  fitted_x <- s_hat_x
  omega_x_hat <- theta_hat[omega_x_indx]
  if(d_x == 0){
    sd_x_hat <- omega_x_hat  
  }else{
    sd_x_hat <- (t(sapply(X = s_hat_x, FUN = function(x){x^(0:d_x)})) %*% omega_x_hat)
  }
  residuals_x <- (x - fitted_x) / sd_x_hat
  # Scatter plot of reference method residuals
  plot(x = fitted_x, y = residuals_x, pch = 16, col = adjustcolor("black", 0.5),
       xlab = "Fitted Values", ylab = "Standardized Residuals", main = "")
  abline(h = 0, col = "red")
  # Histogram of reference method residuals
  hist(residuals_x, xlab = "Standardized Residuals", main = "Reference Method Diagnostics")
  # QQ-plot of reference method residuals
  qqx <- qqPlot(residuals_x, col.lines = "red", xlab = "Normal Quantiles", ylab = "Residual Quantiles", main = "")
  # Calculate comparator method residuals
  s_hat_y <- rep(fitted_model$BLA, times = r_y)
  beta_hat <- theta_hat[beta_indx]
  fitted_y <- (t(sapply(X = s_hat_y, FUN = function(x){x^(0:p)})) %*% beta_hat)
  omega_y_hat <- theta_hat[omega_y_indx]
  if(d_y == 0){
    sd_y_hat <- omega_y_hat  
  }else{
    sd_y_hat <- (t(sapply(X = s_hat_y, FUN = function(x){x^(0:d_y)})) %*% omega_y_hat)
  }
  residuals_y <- (y - fitted_y) / sd_y_hat
  # Scatter plot of comparator method residuals
  plot(x = fitted_y, y = residuals_y, pch = 16, col = adjustcolor("black", 0.5),
       xlab = "Fitted Values", ylab = "Standardized Residuals", main = "")
  abline(h = 0, col = "red")
  # Histogram of comparator method residuals
  hist(residuals_y, xlab = "Standardized Residuals", main = "Comparator Method Diagnostics")
  # QQ-plot of comparator method residuals
  qqy <- qqPlot(residuals_y, col.lines = "red", xlab = "Normal Quantiles", ylab = "Residual Quantiles", main = "")
}



