
##' @export
##' @rdname robustness_value
robustness_value <- sensemakr::robustness_value

##' @export
##' @rdname robustness_value
extreme_robustness_value <- sensemakr::extreme_robustness_value

##'@export
##'@rdname robustness_value
xrv <- sensemakr::xrv

##'@export
##'@rdname robustness_value
rv <- sensemakr::rv

##' Computes the (extreme) robustness value for IV
##'
##' @description
##' Computes robustness values for \code{\link{iv_fit}} objects, adapting the
##' robustness value definitions of \pkg{sensemakr} to instrumental variables as
##' described in Cinelli and Hazlett (2025).
##' For \code{parm = "iv"}, returns robustness values for the IV estimate; for
##' \code{parm = "fs"} or \code{parm = "rf"}, dispatches to \pkg{sensemakr}
##' methods on the corresponding \code{\link{lm}} models.
##'
##'
##' @param model an \code{\link{iv_fit}} model
##' @param parm  parameter for which the robustness value is computed. Default is \code{iv}, meaning that the robustness value of the IV estimate is computed. Other options are to compute the robustness value of auxiliary estimates, such as the first stage (\code{fs}) or the reduced form (\code{rf}).
##' @param q percent change of the effect estimate that would be deemed problematic. Default is 1, which means a reduction (increase) of 100\% of the current effect estimate (bring estimate to zero). It has to be greater than zero.
##' @param alpha significance level.
##' @param min  in many cases, researchers are interested in biases as large or larger than a certain amount (for instance, the strength of confounding to bring a positive estimate to zero or below). Setting \code{min = TRUE} (default) computes the robustness value for such cases. Setting \code{min = FALSE} computes the robustness value for a bias of exactly \code{q}.
##' @param ... further arguments passed to or from other methods.
##'
##' @returns A numeric value with the (extreme) robustness value.
##'
##' @examples
##' data("card")
##' y <- card$lwage
##' d <- card$educ
##' z <- card$nearc4
##' x <- model.matrix( ~ exper + expersq + black + south + smsa + reg661 + reg662 +
##'                      reg663 + reg664 + reg665+ reg666 + reg667 + reg668 + smsa66,
##'                    data = card)
##' card.fit <- iv_fit(y, d, z, x)
##'
##' # robustness value of the IV estimate
##' rv(card.fit)
##'
##' # extreme robustness value
##' xrv(card.fit)
##'
##' # robustness values for first-stage and reduced-form
##' rv(card.fit, parm = "fs")
##' rv(card.fit, parm = "rf")
##'
##' @references
##'
##' Cinelli, C. and Hazlett, C. (2025), "An Omitted Variable Bias Framework for Sensitivity Analysis of Instrumental Variables." Biometrika. \doi{10.1093/biomet/asaf004}
##'
##' @exportS3Method sensemakr::extreme_robustness_value iv_fit
##' @rdname robustness_value
extreme_robustness_value.iv_fit <- function(model,
                                            parm = "iv",
                                            q=1,
                                            alpha=0.05,
                                            min = TRUE, ...){

  check_alpha(alpha)
  # sensemakr:::check_q(q)
  parm <- match.arg(parm,
                    choices = c("iv",  "fs", "rf"),
                    several.ok = FALSE)
  if(parm == "iv"){
    args <- iv_model_helper(fs = model$models$fs, rf = model$models$rf, instrument = "z")
    xrv <- do.call("rv_from_summary_stats", c(args, which = "xrv", q = q, alpha = alpha, min = min))
    return(xrv)
  } else {
    xrv <- sensemakr::xrv(model$models[[parm]],covariates = "z", q = abs(q), alpha = alpha)
    names(xrv) <- parm
    return(xrv)
  }
}



##' @exportS3Method sensemakr::robustness_value iv_fit
##' @rdname robustness_value
robustness_value.iv_fit <- function(model, parm = "iv", q=1, alpha=0.05, min = TRUE, ...){
  check_alpha(alpha)
  # sensemakr:::check_q(q)
  parm <- match.arg(parm,
                    choices = c("iv",  "fs", "rf"),
                    several.ok = FALSE)
  if(parm == "iv"){
    args <- iv_model_helper(fs = model$models$fs, rf = model$models$rf, instrument = "z")
    rv <- do.call("rv_from_summary_stats", c(args, which = "rv", q = q, alpha = alpha, min = min))
    return(rv)
  } else {
    rv <- sensemakr::rv(model$models[[parm]],covariates = "z", q = abs(q), alpha = alpha)
    names(rv) <- parm
    return(rv)
  }
}



check_alpha <- function(alpha) {
  # Error: alpha, if provided, was non-numeric or out of bounds
  if ((!is.numeric(alpha) || length(alpha) > 1 ||
       alpha < 0 || alpha > 1)) {
    stop("`alpha` must be between 0 and 1.")
  }
}


# internal functions ------------------------------------------------------


rv_from_summary_stats <- function(fs.coef,
                                  fs.se,
                                  rf.coef,
                                  rf.se,
                                  rho,
                                  dof,
                                  which = c("rv", "xrv"),
                                  q = 1,
                                  alpha = 0.05,
                                  min = TRUE){
  which <- match.arg(which)
  fun <- switch(which,
                rv  = sensemakr::rv,
                xrv = sensemakr::xrv)
  tau.r    <- rf.coef/fs.coef
  tau0     <- (1 - q)*tau.r
  phi.r    <- rf.coef - tau0*fs.coef
  se.phi.r <- sqrt(rf.se^2 + tau0^2*fs.se^2 - 2*tau0*rho*fs.se*rf.se)
  t.phi.r  <- abs(phi.r/se.phi.r)
  t.fs.r   <- abs(fs.coef/fs.se)

  rv.phi <- fun(t_statistic = t.phi.r, dof = dof, q = 1, alpha = alpha)

  rv.fs  <- fun(t_statistic = t.fs.r, dof = dof, q = 1, alpha = alpha)
  if(min){
    rv.iv <- min(rv.phi, rv.fs)
  } else {
    rv.iv <- rv.phi
  }
  attr(rv.iv,"q") <- q
  attr(rv.iv,"min") <- min
  attr(rv.iv,"alpha") <- alpha
  class(rv.iv) <- c("numeric", "rv.iv")
  names(rv.iv) <- "iv"
  return(rv.iv)
}


# sensitivity_stats -------------------------------------------------------

##' Sensitivity statistics for instrumental variable estimates
##'
##' @description
##' Convenience function that computes robustness values for IV estimates as well as auxiliary first stage and reduced form regressions.
##'
##' @returns A \code{\link{data.frame}} with columns for the estimate, confidence interval
##' bounds (lower and upper), t-value, extreme robustness value (\code{xrv_qa}),
##' robustness value (\code{rv_qa}), and the parameters used (\code{q}, \code{min},
##' \code{alpha}, \code{dof}).
##'
##' @examples
##' data("card")
##' y <- card$lwage
##' d <- card$educ
##' z <- card$nearc4
##' x <- model.matrix( ~ exper + expersq + black + south + smsa + reg661 + reg662 +
##'                      reg663 + reg664 + reg665+ reg666 + reg667 + reg668 + smsa66,
##'                    data = card)
##' card.fit <- iv_fit(y, d, z, x)
##'
##' # sensitivity statistics for the IV estimate
##' sensitivity_stats(card.fit)
##'
##' # sensitivity statistics for the first-stage
##' sensitivity_stats(card.fit, parm = "fs")
##'
##' @export
sensitivity_stats <- sensemakr::sensitivity_stats

##' @inheritParams sensemakr
##' @inheritParams ovb_contour_plot
##' @param ... further arguments passed to or from other methods.
##' @exportS3Method sensemakr::sensitivity_stats iv_fit
##' @export
##' @rdname sensitivity_stats
sensitivity_stats.iv_fit <- function(model, parm = "iv", q = 1, alpha = 0.05, min = TRUE, ...){

  parm <- match.arg(parm,
                    choices = c("iv", "fs","rf"),
                    several.ok = FALSE)

  check_alpha(alpha)
  # sensemakr:::check_q(q)

  fit2 <- iv_fit(y = model$data$y,
                       d = model$data$d,
                       z = model$data$z,
                       x = as.matrix(model$data[,-c(1:3)]),
                       h0 = (1-q)*model$estimates$iv$estimate,
                       alpha = alpha)
  ests <- fit2$estimates[parm]
  estimate <- get(ests, "estimate")
  conf.int <- t(get(ests, "conf.int", fun = range))
  colnames(conf.int) <- c("lwr","upr")
  t.value  <- get(ests, "t.value")
  rv_qa    <- sapply(names(ests),  function(x) rv(fit2, parm = x, q=q, alpha = alpha))
  xrv_qa   <- sapply(names(ests), function(x) xrv(fit2, parm = x, q=q, alpha = alpha))
  dof      <- get(ests, "dof")
  alpha    <- get(ests, "alpha")
  # h0       <- get(ests, "h0")

  sensitivity_stats      <- as.data.frame(cbind(estimate, conf.int, t.value, xrv_qa, rv_qa, #h0,
                                                q= q, min = min, alpha, dof))
  return(sensitivity_stats)
}

##' @inheritParams sensemakr
##' @inheritParams ovb_contour_plot
##' @exportS3Method sensemakr::sensitivity_stats iv.sensemakr
##' @export
##' @rdname sensitivity_stats
sensitivity_stats.iv.sensemakr <- function(model, parm = "iv", q = 1, alpha = 0.05, min = TRUE, ...){
  do.call("rbind", model$sensitivity_stats[parm])
}

get <- function(ests, what, fun = function(x) x, names = NULL){
  out <- sapply(ests, function(x) fun(x[[what]]))
  if (!is.null(names)) {
    rownames(out) <- names
  }
  out
}

# print -------------------------------------------------------------------

##' @export
print.rv.iv <- function(x, digits = 3,...){
  value <- x
  attributes(value) <- list(names = names(value))
  class(value) <- "numeric"
  print(value, digits = digits)
  q <- attr(x, "q")
  alpha <- attr(x, "alpha")
  min <- attr(x, "min")
  if(min){
    cat("Parameters: q >=", q)
  } else {
    cat("Parameters: q =", q)
  }
  if (!is.null(alpha)) cat(", alpha =", alpha,"\n")
}

printSenStats <- function(x, digits = 3, note = TRUE){
  cat("\n")
  cat("Sensitivity Statistics:")
  cat("\n")
  q <- x$q
  min <- x$min
  alpha <- x$alpha
  rv <- x$rv_qa
  xrv <- x$xrv_qa
  cat("  Extreme Robustness Value:", format(xrv, digits = digits))
  cat("\n")
  cat("  Robustness Value:",format(rv, digits = digits))
  cat("\n")
  if(note){
    cat("Note:",
        paste0("q", ifelse(min, " >= ", " = "), q, ","),
        paste0("alpha = ", alpha,"."))
    cat("\n")
  }
  cat("\n")
}
