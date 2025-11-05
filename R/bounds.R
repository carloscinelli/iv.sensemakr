
# iv bounder --------------------------------------------------------------




ovb4iv_partial_r2_bound <- function(fs,
                                    rf,
                                    instrument,
                                    benchmark,
                                    kd = 1,
                                    ky = kd,
                                    alpha = 0.05,
                                    reduce = TRUE){

  # check if lower limit or upper limit of CI
  tau.r      <- coef(rf)[instrument]/coef(fs)[instrument]
  sign.tau.r <- sign(tau.r)

  if ( (sign.tau.r > 0 & reduce) | (sign.tau.r < 0 & !reduce) ) {
    ci.limit <- "lower"
  } else {
    ci.limit <- "upper"
  }


  r2yxj.x <- rep(NA, length(benchmark))
  for(i in seq_along(benchmark)){
    partials <- get_partials(fs = fs, rf = rf, benchmark = benchmark[i])

    # benchmark for Y
    r2yxj.x[i] <- do.call("maxR2", partials)$objective
  }

  # benchmark for d
  z.form <- as.formula(paste(instrument, "~ . -", instrument))
  z.reg   <- update(rf, z.form)
  r2dxj.x  <- sensemakr::partial_r2(z.reg, covariates = benchmark)

  bounds <- vector(mode = "list", length = length(benchmark))

  for(i in seq_along(benchmark)){
    # label
    label <- sensemakr:::label_maker(benchmark_covariate = benchmark[i], kd = kd, ky = ky, digits = 2)

    bounds[[i]] <- sensemakr:::ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i],
                                                            r2yxj.dx = r2yxj.x[i],
                                                            kd = kd,
                                                            ky = ky,
                                                            bound_label = label)
    names(bounds[[i]])[2:3] <- c("r2zw.x", "r2y0w.zx")
    bounds[[i]]$bound_value <- iv_adjusted_limit(fs = fs,
                                                 rf = rf,
                                                 instrument = instrument,
                                                 r2zw.x = bounds[[i]]$r2zw.x,
                                                 r2y0w.zx = bounds[[i]]$r2y0w.zx,
                                                 alpha = alpha,
                                                 ci.limit = ci.limit)
  }
  bounds <- do.call("rbind", bounds)

  return(bounds)
}



# function residualizes  y, d, and xj
get_partials <- function(fs, rf, benchmark){
  form <- as.formula(paste("~ . -", benchmark))
  fsj  <- update(fs, form)
  rfj  <- update(rf, form)

  xj.form <- as.formula(paste(benchmark, "~ . -", benchmark))
  xj.reg   <- update(rf, xj.form)

  y.zx  <- resid(rfj)
  d.zx  <- resid(fsj)
  xj.zx <- resid(xj.reg)

  out <- list(y.zx  = y.zx,
              d.zx  = d.zx,
              xj.zx = xj.zx)
  return(out)
}

# function finds the maximum correlation of y0.zx with xj.zx
# for any tau0
# for now we do it numerically (which is pretty fast).
# change to algebraic solution later.
maxR2 <- function(y.zx, d.zx, xj.zx) {
  cor.y0xj <- function(tau0) {
    ytau0.zx <- y.zx - tau0*d.zx
    cor(ytau0.zx, xj.zx)^2
  }
  max.r2 <- optimize(cor.y0xj, interval = c(-1e5, 1e5), maximum = T)
  return(max.r2)
}

# these functions are for plotting only.
# not exported for now
iv_adjusted_limit <- function(...){
  UseMethod("iv_adjusted_limit")
}


iv_adjusted_limit.lm <- function(fs,
                                 rf,
                                 instrument,
                                 r2zw.x,
                                 r2y0w.zx,
                                 alpha = 0.05,
                                 ci.limit = c("lower", "upper"),
                                 ...){
  call.args   <- list(r2zw.x,
                      r2y0w.zx,
                      alpha = alpha,
                      ci.limit = ci.limit,
                      ...)
  iv.data     <- iv_model_helper(fs = fs, rf = rf, instrument = instrument)
  args        <- c(iv.data, call.args)
  out         <- do.call(iv_adjusted_limit, args)
  return(out)
}


iv_adjusted_limit.numeric <- function(fs.coef,
                                      fs.se,
                                      rf.coef,
                                      rf.se,
                                      rho,
                                      dof,
                                      r2zw.x,
                                      r2y0w.zx,
                                      alpha = 0.05,
                                      ci.limit = c("lower", "upper"),
                                      ...){

  ci.limit <- match.arg(ci.limit)

  fun <- switch(ci.limit,
                lower = quad_lower_limit,
                upper = quad_upper_limit)

  ts  <- critical_t(alpha = alpha,
                    dof = dof - 1)

  CIF <- max_ci_factor(r2dz.x = r2zw.x,
                       r2yz.dx = r2y0w.zx,
                       crit_f = ts/sqrt(dof - 1))

  abc <- abc_builder(fs.coef  = fs.coef,
                     fs.se    = fs.se,
                     rf.coef  = rf.coef,
                     rf.se    = rf.se,
                     rho      = rho,
                     crit.thr = CIF*sqrt(dof))

  limits <- do.call("fun", abc)

  return(limits)
}
