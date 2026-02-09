

# Generic Plot ------------------------------------------------------------
##' Sensitivity analysis plots for IV sensemakr
##'
##' This function provides the contour plots of the sensitivity analysis results obtained with the function \code{\link{sensemakr}} for IV. It is basically a dispatcher to the core plot function \code{\link{ovb_contour_plot}}.
##'
##' @param x an object of class \code{iv.sensemakr} created with the \code{\link{sensemakr}} function.
##' @inheritParams ovb_contour_plot
##'@export
plot.iv.sensemakr = function(x,
                             sensitivity.of = c("ci", "lwr", "upr", "t-value"),
                             parm = "iv",
                             ...) {

  # # Capture additional arguments
  # additional_args <- list(...)
  #
  # # Resolve conflicts: Override x$pars with additional_args if conflicts exist
  # for (arg in names(additional_args)) {
  #   if (arg %in% names(x$pars)) {
  #     x$pars[[arg]] <- additional_args[[arg]]
  #   }
  # }
  #
  # type <- match.arg(type)
  # sensitivity.of <- match.arg(sensitivity.of)
  #
  # # Call the dispatch function of interest
  # if(sensitivity.of == "ci"){
  #   oldpar <- par(mfrow = c(1,2))
  #   ovb_contour_plot(x$unadjusted,
  #                    benchmark_covariates = x$pars$benchmark_covariates,
  #                    kz = x$pars$kz,
  #                    ky = x$pars$ky,
  #                    kd = x$pars$kd,
  #                    sensitivity.of = "lwr",
  #                    main = "Lower Limit CI",
  #                    parm = parm,
  #                    r2zw.x = x$pars$r2zw.x,
  #                    r2y0w.zx = x$pars$r2y0w.zx,
  #                    bound_label = x$pars$bound_label, ...)
  #
  #   ovb_contour_plot(x$unadjusted,
  #                    benchmark_covariates = x$pars$benchmark_covariates,
  #                    kz = x$pars$kz,
  #                    ky = x$pars$ky,
  #                    kd = x$pars$kd,
  #                    main = "Upper Limit CI",
  #                    sensitivity.of = "upr",
  #                    parm = parm,
  #                    r2zw.x = x$pars$r2zw.x,
  #                    r2y0w.zx = x$pars$r2y0w.zx,
  #                    bound_label = x$pars$bound_label, ...)
  #
  #   on.exit(par(oldpar))
  # }
  # else{
  #   ovb_contour_plot(x$unadjusted,
  #                    benchmark_covariates = x$pars$benchmark_covariates,
  #                    kz = x$pars$kz,
  #                    ky = x$pars$ky,
  #                    kd = x$pars$kd,
  #                    sensitivity.of = sensitivity.of,
  #                    parm = parm,
  #                    r2zw.x = x$pars$r2zw.x,
  #                    r2y0w.zx = x$pars$r2y0w.zx,
  #                    bound_label = x$pars$bound_label, ...)
  # }

  additional_args <- list(...)
  pars <- c("benchmark_covariates",
            "kz", "ky", "kd", "r2zw.x",
            "r2y0w.zx", "bound_label")
  combined_args <- modifyList(x$pars[pars], additional_args)

  sensitivity.of <- match.arg(sensitivity.of)

  if (sensitivity.of == "ci") {
    oldpar <- par(mfrow = c(1,2))

    do.call("ovb_contour_plot", c(list(x$unadjusted,
                                       sensitivity.of = "lwr",
                                       main = "Lower Limit CI",
                                       parm = parm),
                                  combined_args))

    do.call("ovb_contour_plot", c(list(x$unadjusted,
                                       sensitivity.of = "upr",
                                       main = "Upper Limit CI",
                                       parm = parm),
                                  combined_args))

    on.exit(par(oldpar))
  } else {
    do.call("ovb_contour_plot", c(list(x$unadjusted,
                                       sensitivity.of = sensitivity.of,
                                       parm = parm),
                                  combined_args))
  }
}


# contour plots -----------------------------------------------------------

##' Contour plots of omitted variable bias for IV
##'
##' @description
##' Contour plots of omitted variable bias for sensitivity analysis of instrumental variable estimates.
##'
##' The main inputs are an \code{\link{iv_fit}} model, and the covariates used for benchmarking the strength of omitted variables.
##'
##' If \code{parm = "iv"} (default) contour plots of the IV estimate are shown. The horizontal axis of the plot shows hypothetical values of the partial R2 of latent variables with the instrument. The vertical axis shows hypothetical values of the partial R2 of latent variables with the (pot.) outcome. The contour levels represent the adjusted lower limit (or upper limit) of the Anderson-Rubin confidence interval of the IV estimate, or the t-value for testing a specific null hypothesis.  The reference points are the bounds on the partial R2 of latent variables if they were k times “as strong” as the observed covariate used for benchmarking (see arguments kz and ky). The dotted red line show the chosen critical threshold (for instance, zero): latent variables with such strength (or stronger) are sufficient to invalidate the research conclusions.
##'
##' if \code{parm = "fs"} or \code{parm = "rf"}, then contour plots of the first-stage and reduced-form regression are shown. See, e.g, \code{\link[sensemakr]{ovb_contour_plot.lm}}.
##'
##' See Cinelli and Hazlett (2020, 2025) for details.
##' @references
##'
##'  Cinelli, C. and Hazlett, C. (2020), "Making Sense of Sensitivity: Extending Omitted Variable Bias." Journal of the Royal Statistical Society, Series B (Statistical Methodology).
##'
##' Cinelli, C. and Hazlett, C. (2025), "An Omitted Variable Bias Framework for Sensitivity Analysis of Instrumental Variables." Biometrika. \doi{10.1093/biomet/asaf004}
##'
##' @export
ovb_contour_plot <- sensemakr::ovb_contour_plot

##' @inheritParams sensemakr
##' @param ... further arguments and graphical parameters.
##' @param sensitivity.of should the contour plot show adjusted lower limits of confidence intervals (\code{"lwr"}), upper limit of confidence intervals (\code{"upr"}) or t-values (\code{"t-value"})?
##' @param parm contour plots of which estimate? Options are \code{iv} for instrumental variable estimates, \code{fs} for first-stage estimates, and \code{rf} for reduced-form estimates.
##' @param xlab label of x axis. If `NULL`, default label is used.
##' @param ylab label of y axis. If `NULL`, default label is used.
##'
##' @details
##' Other parameters include:
##' \describe{
##'  \item{\code{alpha}}{significance level.}
##'  \item{\code{threshold}}{ critical threshold, default is \code{0}.}
##'  \item{\code{lim}}{limits for the axes.}
##'  \item{\code{lim.x}}{limits for the x axis. Default is \code{lim}.}
##'  \item{\code{lim.y}}{limits for the y axis. Default is \code{lim}.}
##'  \item{\code{nlevels}}{number of levels in the contour plot.}
##'  \item{\code{col.contour}}{color of the contour lines.}
##'  \item{\code{col.thr.line}}{color of the threshold line.}
##'  \item{\code{label.text}}{should benchmark label texts be shown? Default is \code{TRUE}.}
##'  \item{\code{cex.label.text}}{character size of label text. Default is \code{.7}.}
##'  \item{\code{label.bump.x}}{bump on the x coordinate of label text.}
##'  \item{\code{label.bump.y}}{bump on the y coordinate of label text.}
##'  \item{\code{cex.lab}}{The magnification to be used for x and y labels relative to the current setting of cex.}
##'  \item{\code{cex.main}}{The magnification to be used for main titles relative to the current setting of cex.}
##'  \item{\code{cex.axis}}{The magnification to be used for axis annotation relative to the current setting of cex.}
##'  \item{\code{asp}}{the y/x aspect ratio. Default is 1.}
##' }
##' If \code{parm = "fs"} or \code{parm = "rf"} the function is simply a wrapper to the sensemakr function \code{\link[sensemakr]{ovb_contour_plot.lm}} on the first-stage or reduced-form \code{\link{lm}} models.
##'
##'
##'
##' @exportS3Method sensemakr::ovb_contour_plot iv_fit
##' @rdname ovb_contour_plot
ovb_contour_plot.iv_fit <- function(model,
                                    benchmark_covariates = NULL,
                                    kz = 1,
                                    ky = kz,
                                    kd = kz,
                                    sensitivity.of = c("lwr", "upr", "t-value"),
                                    parm = "iv",
                                    r2zw.x = NULL,
                                    r2y0w.zx = r2zw.x,
                                    bound_label = "manual bound",
                                    xlab = NULL,
                                    ylab = NULL,
                                    ...) {
  x <- model
  parm <- match.arg(parm,
                    choices = c("iv",  "fs", "rf"),
                    several.ok = FALSE)
  sensitivity.of <- match.arg(sensitivity.of)

  if (parm == "iv") {
    if (sensitivity.of %in% c("lwr", "upr")) {
      ovb4iv_contour_plot.iv_fit(model = x,
                                 benchmark_covariates = benchmark_covariates,
                                 kz = kz, ky = ky, ci.limit = sensitivity.of,
                                 xlab = xlab, ylab = ylab,
                                 r2zw.x = r2zw.x,
                                 r2y0w.zx = r2y0w.zx,
                                 bound_label = bound_label,
                                 ...)
    } else {
      if (is.null(xlab)) {
        xlab <- expression(paste("Partial ", R^2, " of omitted variable(s) with the instrument"))
      }
      if (is.null(ylab)){
        ylab <- expression(paste("Partial ", R^2, " of omitted variable(s) with the pot. outcome"))
      }
      sensemakr::ovb_contour_plot(model = x$models$ar,
                                  sensitivity.of = sensitivity.of,
                                  treatment = "z",
                                  benchmark_covariates = benchmark_covariates,
                                  kd = kz,
                                  ky = kd,
                                  xlab = xlab,
                                  ylab = ylab,
                                  ...)
    }
  }

  if (parm == "fs") {
    if (is.null(xlab)) {
      xlab <- expression(paste("Partial ", R^2, " of omitted variable(s) with the instrument"))
    }
    if (is.null(ylab)){
      ylab <- expression(paste("Partial ", R^2, " of omitted variable(s) with the treatment"))
    }
    sensemakr::ovb_contour_plot(model = x$models$fs,
                                sensitivity.of = sensitivity.of,
                                treatment = "z",
                                benchmark_covariates = benchmark_covariates,
                                kd = kz,
                                ky = kd,
                                xlab = xlab,
                                ylab = ylab,
                                ...)
  }

  if (parm == "rf") {
    if (is.null(xlab)) {
      xlab <- expression(paste("Partial ", R^2, " of omitted variable(s) with the instrument"))
    }
    if (is.null(ylab)){
      ylab <- expression(paste("Partial ", R^2, " of omitted variable(s) with the outcome"))
    }
    sensemakr::ovb_contour_plot(model = x$models$rf,
                                sensitivity.of = sensitivity.of,
                                treatment = "z",
                                benchmark_covariates = benchmark_covariates,
                                kd = kz,
                                ky = ky,
                                xlab = xlab,
                                ylab = ylab,
                                ...)
  }
}



ovb4iv_contour_plot.iv_fit <- function(model,
                                       benchmark_covariates = NULL,
                                       kz = 1,
                                       ky = kz,
                                       alpha = 0.05,
                                       ci.limit = c("lwr", "upr"),
                                       round = 3,
                                       threshold = model$pars$h0,
                                       r2zw.x = NULL,
                                       r2y0w.zx = r2zw.x,
                                       bound_label = "Manual Bound",
                                       ...){

  ci.limit  <- match.arg(ci.limit)

  call.args   <- list(alpha = alpha,
                      round = round,
                      threshold = threshold,
                      ci.limit = ci.limit,
                      ...)
  fs <- model$models$fs
  rf <- model$models$rf
  iv.data     <- iv_model_helper(fs = fs,
                                 rf = rf,
                                 instrument = "z")
  args        <- c(iv.data, call.args)
  bounds <- NULL
  if (!is.null(benchmark_covariates)) {
    bounds       <- ovb4iv_partial_r2_bound(model,
                                            benchmark = benchmark_covariates,
                                            kz = kz,
                                            ky = ky,
                                            alpha = alpha,
                                            ci.limit = ci.limit)
  }

  if (!is.null(r2zw.x)) {

    manual_bounds <- data.frame(bound_label = bound_label,
                                r2zw.x = r2zw.x,
                                r2y0w.zx = r2y0w.zx)

    manual_bounds$bound_value <- iv_adjusted_limit(fs = model$models$fs,
                                                   rf = model$models$rf,
                                                   instrument = "z",
                                                   ci.limit = ci.limit,
                                                   r2zw.x = r2zw.x,
                                                   r2y0w.zx = r2y0w.zx,
                                                   alpha = alpha)

    bounds <- rbind(bounds, manual_bounds)
  }

  args        <- c(args, bounds)
  out         <- do.call(ovb4iv_contour_plot.numeric, args)
  return(invisible(out))
}




ovb4iv_contour_plot.numeric <- function(fs.coef,
                                        rf.coef,
                                        fs.se,
                                        rf.se,
                                        rho,
                                        dof,
                                        ci.limit = c("lwr", "upr"),
                                        alpha = 0.05,
                                        threshold = 0,
                                        r2zw.x = NULL,
                                        r2y0w.zx = r2zw.x,
                                        bound_label = "",
                                        bound_value = NULL,
                                        lim = max(c(0.4, r2zw.x + 0.1, r2y0w.zx + 0.1)),
                                        lim.x = lim,
                                        lim.y = lim,
                                        nlevels = 5,
                                        col.contour = "black",
                                        col.thr.line = "red",
                                        label.text = TRUE,
                                        cex.label.text = .7,
                                        label.bump.x = lim.x*(1/15),
                                        label.bump.y = lim.y*(1/15),
                                        xlab = NULL,
                                        ylab = NULL,
                                        cex.lab = .8,
                                        cex.axis = .8,
                                        cex.main = 1,
                                        asp = lim.x/lim.y,
                                        list.par = list(mar = c(4,4,1,1), pty = "s"),
                                        round = 3,
                                        ...) {


  # Set up the grid for the contour plot
  grid_values.x = seq(0, lim.x, by = lim.x/400)
  grid_values.y = seq(0, lim.y, by = lim.y/400)

  # check if lwr limit or upr limit of CI
  tau.r      <- rf.coef/fs.coef
  sign.tau.r <- sign(tau.r)

  ci.limit <- match.arg(ci.limit)



  FUN <- function(X, Y){
    iv_adjusted_limit( fs.coef = fs.coef,
                       rf.coef = rf.coef,
                       fs.se   = fs.se,
                       rf.se   = rf.se,
                       rho     = rho,
                       dof     = dof,
                       r2zw.x = X,
                       r2y0w.zx = Y,
                       alpha = alpha,
                       ci.limit = ci.limit)
  }

  z_axis = outer(X = grid_values.x,  Y = grid_values.y, FUN = FUN)


  out = list(r2zw.x = grid_values.x,
             r2y0w.zx = grid_values.y,
             value = z_axis)


  default_levels    <- quantile(z_axis[is.finite(z_axis)], seq(0, 1, length.out = max(2, nlevels)))

  # set minimum to 1% percentile, this will be the -Infinite if there are infinities
  default_levels[1] <- quantile(z_axis[is.finite(z_axis)], 0.01)

  # set maximum to 99% percentile, this will be the +Infinite if there are infinities
  default_levels[length(default_levels)] <- quantile(z_axis[is.finite(z_axis)], 0.99)

  # unique, sort and round for prettier labels
  contour.lines <- list()
  while (length(contour.lines) == 0) {
    try.levels <- round(sort(unique(c(default_levels, threshold))), round)

    # compute actual levels that contour will use
    contour.lines <- contourLines(x = grid_values.x, y = grid_values.y, z = z_axis, levels = try.levels)
    contour.levels <- unique(sapply(contour.lines, function(x) x$level))
    round <- round + 1
  }
  lab.th2 <- th2 <- NULL

  if ( ci.limit == "lwr") {
    min <- min(z_axis)
    if (!is.finite(min)) {
      th2 <- min(contour.levels)
      lab.th2 <- "-Inf"
    }

  } else {

    max <- max(z_axis)
    if (!is.finite(max)) {
      th2 <- max(contour.levels)
      lab.th2 <- "Inf"

    }
  }


  # Plot contour plot:
  if (is.null(list.par)) {
    oldpar <- par(mar = c(5, 5, 4, 1) + .1)
    on.exit(par(oldpar))
  } else {
    if (!is.list(list.par)) stop("list.par needs to be a named list")
    oldpar <- do.call("par", list.par)
    on.exit(par(oldpar))
  }


  if (is.null(xlab)) {
    xlab <- expression(paste("Partial ", R^2, " of omitted variable(s) with the instrument"))
  }

  if (is.null(ylab)) {
    ylab <- expression(paste("Partial ", R^2, " of omitted variable(s) with the pot. outcome"))
  }

  too_close <- abs(contour.levels - threshold) < min(diff(contour.levels[contour.levels!=round(threshold, round)]))*.25
  line_color <- rep("black", length(contour.levels))
  line_color[contour.levels == th2] <- "transparent"
  line_color[too_close] <- "transparent"
  line_color[contour.levels == round(threshold, round)] <- "red"


  line_type <- rep(1, length(contour.levels))
  line_type[contour.levels == round(threshold, round)] <- 2

  line_width <- rep(1, length(contour.levels))
  line_width[contour.levels == round(threshold, round)] <- 2


  labels <- contour.levels
  labels[contour.levels  %in% c(th2)] <- lab.th2

  contour(x = grid_values.x,
          y = grid_values.y,
          z = z_axis,
          # nlevels = nlevels,
          levels = contour.levels,
          labels = labels,
          xlab = xlab,
          ylab = ylab,
          cex.main = cex.main,
          cex.lab = cex.lab,
          cex.axis = cex.axis,
          asp = asp,
          col = line_color,
          lty = line_type,
          lwd = line_width,
          ...)

  if (!is.null(th2)) {
    contour(grid_values.x, grid_values.y,
            z_axis,
            levels = c(th2),
            label = c(lab.th2),
            add = TRUE,
            labcex = 1.2*cex.lab,
            col = col.thr.line,
            lwd = 2,
            lty = 2, ...)
  }

  # Add the point of the initial estimate.
  points(0, 0, pch = 17, col = "black", cex = 1)

  text(0.0 + label.bump.x,
       0.00 + label.bump.y,
       paste0("Observed\n(",
              signif(z_axis[1], 2),
              ")"),

       cex = cex.label.text)

  # add bounds
  if (!is.null(r2zw.x)) {

    check_r2(r2dz.x = r2zw.x, r2yz.dx = r2y0w.zx)


    sensemakr::add_bound_to_contour(r2dz.x = r2zw.x,
                                             r2yz.dx = r2y0w.zx,
                                             bound_value = bound_value,
                                             bound_label = bound_label,
                                             label.text = label.text,
                                             label.bump.x = label.bump.x,
                                             label.bump.y = label.bump.y,
                                             cex.label.text = cex.label.text,
                                             round = round)
    out$bounds = data.frame(bound_label = bound_label,
                            r2zw.x = r2zw.x,
                            r2y0w.zx = r2y0w.zx,
                            bound_value = bound_value,
                            stringsAsFactors = FALSE)
  }

  invisible(out)
}



check_r2 <- function(r2dz.x, r2yz.dx) {

  if (!is.null(r2dz.x)) {
    if (any(!is.numeric(r2dz.x) | r2dz.x < 0 | r2dz.x > 1)) {
      stop("partial R2 must be a number between zero and one")
    }
  }

  if (!is.null(r2yz.dx) &&
      any(!is.numeric(r2yz.dx) | r2yz.dx < 0 | r2yz.dx > 1)) {
    stop("partial R2 must be a number between zero and one, if provided")
  }
}



# Functions for fast computation of lwr and upr limit -----------------
# these functions are for fast plotting only, they are not exported.

iv_adjusted_limit <- function(...){
  UseMethod("iv_adjusted_limit")
}


##' @noRd
##' @exportS3Method iv_adjusted_limit lm
iv_adjusted_limit.lm <- function(fs,
                                 rf,
                                 instrument,
                                 r2zw.x,
                                 r2y0w.zx,
                                 alpha = 0.05,
                                 ci.limit = c("lwr", "upr"),
                                 max = TRUE,
                                 ...){
  call.args   <- list(r2zw.x,
                      r2y0w.zx,
                      alpha = alpha,
                      ci.limit = ci.limit,
                      max = max,
                      ...)
  iv.data     <- iv_model_helper(fs = fs, rf = rf, instrument = instrument)
  args        <- c(iv.data, call.args)
  out         <- do.call(iv_adjusted_limit, args)
  return(out)
}


##' @noRd
##' @exportS3Method iv_adjusted_limit numeric
iv_adjusted_limit.numeric <- function(fs.coef,
                                      fs.se,
                                      rf.coef,
                                      rf.se,
                                      rho,
                                      dof,
                                      r2zw.x,
                                      r2y0w.zx,
                                      alpha = 0.05,
                                      ci.limit = c("lwr", "upr"),
                                      max = TRUE,
                                      ...){

  ci.limit <- match.arg(ci.limit)

  fun <- switch(ci.limit,
                lwr = quad_lwr_limit,
                upr = quad_upr_limit)


  t.dagger <- sensemakr::adjusted_critical_value(r2dz.x = r2zw.x,
                                                 r2yz.dx = r2y0w.zx,
                                                 dof = dof,
                                                 alpha = alpha,
                                                 max = max)
  abc <- abc_builder(fs.coef  = fs.coef,
                     fs.se    = fs.se,
                     rf.coef  = rf.coef,
                     rf.se    = rf.se,
                     rho      = rho,
                     crit.thr = t.dagger)

  limits <- do.call("fun", abc)

  return(limits)
}




# iv bounder --------------------------------------------------------------



ovb4iv_partial_r2_bound <- function(model,
                                    benchmark,
                                    kz = 1,
                                    ky = kz,
                                    alpha = 0.05,
                                    ci.limit = c("lwr", "upr"),
                                    value = T
                                    ){

  # check if lwr limit or upr limit of CI
  tau.r      <- model$estimates$iv$estimate
  sign.tau.r <- sign(tau.r)

  ci.limit <- match.arg(ci.limit)


  r2yxj.x <- rep(NA, length(benchmark))
  for(i in seq_along(benchmark)){
    partials <- get_partials(model, benchmark = benchmark[i])

    # benchmark for Y
    r2yxj.x[i] <- do.call("maxR2", partials)$objective
  }

  # benchmark for d
  z  <- model$data$z
  x  <- model$data[, !colnames(model$data) %in% c("y", "d", "z") , drop = FALSE]
  z.reg <- lm(formula(as.data.frame(cbind(z, x))), data = model$data)
  r2dxj.x  <- sensemakr::partial_r2(z.reg, covariates = benchmark)

  bounds <- vector(mode = "list", length = length(benchmark))

  for(i in seq_along(benchmark)){
    # label
    label <- label_maker(benchmark_covariate = benchmark[i], kd = kz, ky = ky, digits = 2)

    bounds[[i]] <- sensemakr::ovb_partial_r2_bound(r2dxj.x = r2dxj.x[i],
                                                            r2yxj.dx = r2yxj.x[i],
                                                            kd = kz,
                                                            ky = ky,
                                                            bound_label = label)
    names(bounds[[i]])[2:3] <- c("r2zw.x", "r2y0w.zx")
    if(value){
      bounds[[i]]$bound_value <- iv_adjusted_limit(fs = model$models$fs,
                                                   rf = model$models$rf,
                                                   instrument = "z",
                                                   r2zw.x = bounds[[i]]$r2zw.x,
                                                   r2y0w.zx = bounds[[i]]$r2y0w.zx,
                                                   alpha = alpha,
                                                   ci.limit = ci.limit)
    }
  }
  bounds <- do.call("rbind", bounds)

  return(bounds)
}

label_maker <- function(benchmark_covariate, kd, ky, digits = 2) {
  # Generate the label text
  variable_text = ifelse(
    is.null(benchmark_covariate),
    "\n",
    paste0(" ", benchmark_covariate)
  )

  multiplier_text = ifelse(
    ky == kd,
    paste0(round(ky, digits = digits)),
    paste0(round(kd, digits = digits), "/", round(ky, digits = digits))
  )

  bound_label = paste0(paste0(multiplier_text, "x"), variable_text)

  return(bound_label)

}

# function residualizes  y, d, and xj
get_partials <- function(model, benchmark){
  x  <- model$data[,!(names(model$data) %in% c("y", "d", "z")), drop = FALSE]
  y  <- model$data$y
  z  <- model$data$z
  d  <- model$data$d
  xj <- x[,benchmark, drop=F]
  xx <- x[,!(colnames(x) %in% benchmark), drop = FALSE]

  y.zx  <- resid(lm(formula(as.data.frame(cbind(y, z, xx))), data = model$data))
  d.zx  <- resid(lm(formula(as.data.frame(cbind(d, z, xx))), data = model$data))
  xj.zx <- resid(lm(formula(as.data.frame(cbind(xj, z, xx))), data = model$data))

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
  max.r2 <- optimize(cor.y0xj, interval = c(-1e5, 1e5), maximum = TRUE)
  return(max.r2)
}


# Quadratic equation solvers ----------------------------------------------

## vectorized function to compute lwr limit
## this works because if a < 0 then it is ifinite
## if a > 0 then subtracting is the lwr limit, regardless of the sign of b (delta is always positive)
quad_lwr_limit <- function(a, b, c){
  n <- check_abc_length(a = a, b = b, c = c)
  lwr_limit <- rep(NA, n)
  delta       <- delta(a = a, b = b, c = c)
  lwr_limit[a <= 0] <- -Inf
  finite <- a > 0 & delta >= 0
  lwr_limit[finite] <- (-b[finite] - sqrt(delta[finite])) / (2*a[finite])
  return(lwr_limit)
}

## vectorized function to compute upr limit
## this works because if a < 0 then it is ifinite
## if a > 0 then summing is the upr limit, regardless of the sign of b (delta is always positive)
quad_upr_limit <- function(a, b, c){
  n <- check_abc_length(a = a, b = b, c = c)
  upr_limit <- rep(NA, n)
  delta       <- delta(a = a, b = b, c = c)
  upr_limit[a < 0] <- Inf
  finite <- a > 0 & delta > 0
  upr_limit[finite] <- (-b[finite] + sqrt(delta[finite])) / (2*a[finite])
  return(upr_limit)
}


## Finds roots using Bhaskara's formula
quad_roots <- function(a, b, c){

  n <- check_abc_length(a = a, b = b, c = c)

  delta <- delta(a = a, b = b, c = c)

  roots <- matrix(nrow = n, ncol = 2)
  root1 <- (-b - sqrt(delta))/(2*a)
  root2 <- (-b + sqrt(delta))/(2*a)
  roots[,1] <- pmin(root1, root2)
  roots[,2] <- pmax(root1, root2)

  return(roots)
}

## delta
delta <- function(a, b, c){
  b^2 - 4*a*c
}



## Finds all values below zero
# quad_ineq <- function(a, b, c){
#   n <- check_abc_length(a = a, b = b, c = c)
#
#   # computes delta
#   delta       <- delta(a = a, b = b, c = c)
#
#   # initialize all NA
#   lwr_limit <- rep(NA, n)
#   upr_limit <- rep(NA, n)
#
#   # if a < 0, infinities
#   lwr_limit[a < 0] <- -Inf
#   upr_limit[a < 0] <- Inf
#
#   finite <- a > 0 & delta > 0
#   lwr_limit[finite] <- (-b[finite] - sqrt(delta[finite])) / (2*a[finite])
#   upr_limit[finite] <- (-b[finite] - sqrt(delta[finite])) / (2*a[finite])
#
#   # check finite cases
#
#
#   # initialize CI
#   CI <- vector(mode = "list", length = n)
#   CI[1:n] <- c()
#
# }

quad_ineq <- function(a, b, c){

  delta <- delta(a, b, c)


  if (a < 0) {

    if (delta < 0) {
      # if a < 0 and delta < 0, the whole curve is below zero

      CI <- c(-Inf, Inf)

    } else {
      # if a < 0 and delta > 0, the curve crosses zero, and we have a disjoint set
      roots <- quad_roots(a, b, c)
      CI.1 <- c(-Inf, roots[1])
      CI.2 <- c(roots[2], Inf)
      CI   <- cbind(CI.1, CI.2)

    }
  } else {

    if (delta < 0) {
      # if a > 0 and delta < 0, the whole curve is above zero

      CI <- NULL

    } else {
      # if a > 0 and delta > 0, only parts of the curve is below zero, and boundaries are given by the roots

      CI <- quad_roots(a, b, c)

    }
  }

  return(CI)
}



# creates quadratic coefficients
abc_builder <- function(fs.coef,
                        fs.se,
                        rf.coef,
                        rf.se,
                        rho,
                        crit.thr) {

  # quadratic coefficients
  a  <- (fs.coef^2 - (fs.se^2)*(crit.thr^2))
  b  <- (2*rho*rf.se*fs.se*(crit.thr^2) - 2*rf.coef*fs.coef)
  c  <- (rf.coef^2 - (rf.se^2)*(crit.thr^2))

  out <- list(a = a,
              b = b,
              c = c)
  return(out)

}


check_abc_length <- function(a,b,c){
  na <- length(a)
  nb <- length(b)
  nc <- length(c)
  if (na != nb | na != nc) {
    stop("Length of a, b and c must be the same.")
  }
  return(invisible(na))
}
