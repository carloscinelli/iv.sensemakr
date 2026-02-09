# Fitting -----------------------------------------------------------------

##' Instrumental Variable Estimation using the Anderson-Rubin approach
##'
##' @description
##'
##' \code{iv_fit} computes instrumental variable estimates and confidence intervals using the Anderson-Rubin (AR) approach (Anderson and Rubin, 1949).  This approach is numerically identical to Fieller's theorem (Fieller, 1954). See Cinelli and Hazlett (2025) for further discussion.
##'
##' The AR point estimate is numerically identical to the point estimate of two-stage least squares (2SLS) and it is given by the ratio of the reduced-form to the first-stage regression coefficient. Confidence intervals, however, are constructed differently. 2SLS is equivalent to using the delta-method to obtain the variance of the ratio estimator, and then proceeding by assuming the ratio is asymptotically normal. This approximation can fail when instruments are "weak." The Anderson-Rubin approach instead uses a test inversion procedure to construct confidence intervals. This procedure has correct coverage regardless of instrument strength, at the (inevitable) cost of eventually obtaining unbounded confidence intervals.
##'
##' @param y \code{\link{numeric}} \code{\link{vector}} with the outcome.
##' @param d \code{\link{numeric}} \code{\link{vector}} with the treatment.
##' @param z \code{\link{numeric}} \code{\link{vector}} with the instrument.
##' @param x (optional) \code{\link{numeric}} \code{\link{matrix}} with observed covariates.
##' @param h0 null hypothesis for the target parameter (the IV estimate).
##' @param alpha significance test for hypothesis tests and confidence intervals.
##'
##' @return
##' An object of class \code{iv_fit}, containing:
##' \describe{
##'  \item{ \code{data} }{A \code{\link{data.frame}} with the data used for fitting the models.}
##'  \item{ \code{models} }{A \code{\link{list}} with the \code{\link{lm}} models used for obtaining the IV estimates. This includes the first-stage (FS), reduced-form (RF), and Anderson-Rubin (AR) regressions.}
##'  \item{ \code{estimates} }{A \code{\link{list}} with the summary information of IV estimates, as well as summary information of the auxiliary estimates of the FS, RF, and AR regression. }
##'  \item{\code{pars}}{A \code{\link{list}} with the parameters of the call, such as the null hypothesis \code{h0} and the significance level \code{alpha}.}
##'  }
##'
##'  @returns
##'  The function returns an object of class \code{iv_fit}.
##'
##' @references
##'
##' Anderson, T.W. and Rubin, H. (1949), Estimation of the parameters of a single equation in a complete system of stochastic equations, Annals of Mathematical Statistics, 20, 46-63.
##'
##' Fieller, E. C. (1954). Some problems in interval estimation. Journal of the Royal Statistical Society: Series B (Methodological), 16(2), 175-185.
##'
##' Cinelli, C. and Hazlett, C. (2025), "An Omitted Variable Bias Framework for Sensitivity Analysis of Instrumental Variables." Biometrika. \doi{10.1093/biomet/asaf004}
##' @export
iv_fit <- function(y, d, z, x = NULL, h0 = 0, alpha = 0.05){

  # capture original variable names
  y_name <- deparse(substitute(y))
  d_name <- deparse(substitute(d))
  z_name <- deparse(substitute(z))

  y <- check_num(y, "y")
  d <- check_num(d, "d")
  z <- check_num(z, "z")
  if (!is.null(x)) {
    x <- x[,apply(x, 2, var) != 0, drop = F]
    x <- cbind("constant" = rep(1,nrow(x)), x)
  } else {
    x <- rep(1, length(y))
  }
  x <- check_mat(x, "x")


  n <- c(ny = length(y), nd = length(d), nz = length(z), nx = nrow(x))

  if(length(unique(n)) != 1){
    stop(paste("Variables y, d, z, x must have the same sample size.\n",
                paste0(" Current sample sizes are ", paste(n, collapse = ","), "."))
         )
  }



  mat <- na.omit(cbind(y,d,z,x))

  if (!is.null(attr(mat, which = "na.action"))) {
    message("Your variables have missing values, and these were omitted.\n",
             "Proceeding with a complete case analysis.")
  }

  # y <- mat[,1]
  # d <- mat[,2]
  # z <- mat[,3]
  # x <- mat[, -c(1:3)]
  #
  # data <- list(y = y,d = d,z = z, x = x)
  data <- as.data.frame(mat)
  # colnames(data)

  # first stage
  fs.formula       <- as.formula(paste0("d ~ -1 + z + ",
                                        paste(sapply(colnames(x), function(x) paste0("`",x,"`")), collapse = " + ")))
  fs               <- lm(fs.formula, data = data)
  fs.estimates     <- suppressWarnings(model.helper(fs, name = "FS Estimates", covariate = "z", alpha = alpha))


  # reduced form
  rf.formula       <- as.formula(paste0("y ~ -1 + z + ",
                                        paste(sapply(colnames(x), function(x) paste0("`",x,"`")), collapse = " + ")))
  rf               <- lm(rf.formula, data = data)
  rf.estimates     <- suppressWarnings(model.helper(rf, name = "RF Estimates", covariate = "z", alpha = alpha))

  # anderson-rubin (h0)
  tau              <- h0
  ar.formula       <- as.formula(paste0("y-tau*d ~ -1 + z + ",
                                        paste(sapply(colnames(x), function(x) paste0("`",x,"`")), collapse = " + ")))
  ar               <- lm(ar.formula, data = data)
  ar.estimates     <- suppressWarnings(model.helper(ar, name = "FS Estimates", covariate = "z", alpha = alpha))

  # IV estimates
  # Hyp. Test
  iv.t.value       <- ar.estimates$t.value
  iv.p.value       <- ar.estimates$p.value

  # point estimate
  iv.estimate      <- unname(coef(rf)["z"]/coef(fs)["z"])

  # confidence interval
  iv.data          <- suppressWarnings(iv_model_helper(fs = fs, rf = rf, instrument = "z"))
  args             <- c(iv.data, alpha = alpha)
  iv.ci            <- do.call(what = "ar_confint", args = args)

  iv.estimates <- list(name = "IV Estimates",
                       estimate = iv.estimate,
                       # se = NA,
                       conf.int = iv.ci,
                       t.value  = iv.t.value,
                       p.value  = iv.p.value,
                       h0  = h0,
                       alpha = alpha,
                       dof = ar$df.residual)

  out <- list()
  out$data         <- data
  out$models       <- list(ar = ar, fs = fs, rf = rf)
  out$estimates    <- list(iv = iv.estimates,
                           ar = ar.estimates,
                           fs = fs.estimates,
                           rf = rf.estimates)
  out$pars <- list(alpha = alpha, h0 = h0,
                   y_name = y_name, d_name = d_name, z_name = z_name)
  class(out) <- "iv_fit"
  return(out)

}

# ar confint
ar_confint <- function(fs.coef,  fs.se, rf.coef, rf.se, rho, dof, alpha = 0.05) {

  # critical t
  # ts <- qt(p = 1 - (alpha)/2, df = dof)
  ts <- sqrt(qf(p = log(alpha),
                df1 = 1, df2 = dof, log.p = T, lower.tail = F))
  # names <- paste0(c(1-))

  # quadratic coefficients
  a  <- (fs.coef^2 - (fs.se^2)*(ts^2))
  b  <- (2*rho*rf.se*fs.se*(ts^2) - 2*rf.coef*fs.coef)
  c  <- (rf.coef^2 - (rf.se^2)*(ts^2))
  delta <- b^2 - 4*a*c

  # checks cases
  if ( a < 0 & delta < 0) {
    return(c(-Inf, Inf))
  }
  if (a < 0 & delta > 0) {
    roots <- range(c(-b + sqrt(delta))/(2*a), (-b - sqrt(delta))/(2*a))
    ci <- c(-Inf, min(roots), max(roots), Inf)
    return(ci)
  }

  if ( a >  0 & delta < 0) {
    return(NULL)
  }

  if ( a > 0 & delta >= 0) {
    roots <- range(c(-b + sqrt(delta))/(2*a), (-b - sqrt(delta))/(2*a))
    return(roots)
  }
}





# Helpers and checker -----------------------------------------------------

check_num <- function(x, name = "y"){
  if(!is.numeric(x)) stop(paste(name, "must be a numeric vector."))
  if(is.matrix(x) && ncol(x) > 1) stop(paste(name, "has", ncol(x), "columns. It must be a numeric vector."))
  if(is.matrix(x) && ncol(x) == 1) x <- c(x)
  return(x)
}

check_mat <- function(x, name = "x"){
  if(!is.numeric(x)) stop(paste(name, "must be a numeric matrix."))
  if(!is.matrix(x)) {
    x <- cbind(x)
  }
  if(is.null(colnames(x))) {
    colnames(x) <- paste0("x", seq_along(x))
  }
  return(x)
}


model.helper <- function(model, name, covariate = "z", alpha = 0.05){
  summ <- suppressWarnings(coef(summary(model))[covariate,])
  out <- list()
  out$name      <- name
  out$estimate  <- unname(summ["Estimate"])
  out$se        <- unname(summ["Std. Error"])
  out$conf.int  <- unname(confint(model, level = 1-alpha)[covariate,])
  out$t.value   <- unname(summ["t value"])
  out$p.value   <- unname(summ["Pr(>|t|)"])
  out$h0        <- 0
  out$alpha     <- alpha
  out$dof       <- model$df.residual
  return(out)
}


iv_model_helper <- function(fs, rf, instrument, ...){

  summ.fs <- suppressWarnings(sensemakr::model_helper(model = fs, covariates = instrument))
  summ.rf <- suppressWarnings(sensemakr::model_helper(model = rf, covariates = instrument))
  rho     <- rho(fs, rf)

  if (summ.fs$dof != summ.rf$dof) stop("Degrees of freedom of first-stage and reduced form regressions differ.")

  iv.data <- list(fs.coef   = summ.fs$estimate,
                  fs.se     = summ.fs$se,
                  rf.coef   = summ.rf$estimate,
                  rf.se     = summ.rf$se,
                  rho       = rho,
                  dof       = summ.rf$dof)
  return(iv.data)
}


# coef, confint, summary --------------------------------------------------

##' Extracts point estimates and confidence intervals of an \code{iv_fit} model.
##'
##' @description
##'
##'
##' The function \code{coef} extracts point estimates of an \code{\link{iv_fit}} model.
##'
##' The function \code{confint} extracts confidence intervals of an \code{\link{iv_fit}} model.
##'
##'
##' @param object an object of class \code{\link{iv_fit}}.
##' @param parm which estimate to return. Options are \code{"iv"} for instrumental variable estimate, \code{"fs"} for the first-stage estimate and \code{"rf"} for the reduced-form estimate.
##' @param ... arguments passed to other methods.
##'
##' @returns \code{coef} returns a numeric vector with the estimates of interest.
##'
##' @examples
##' # prepare data
##' data("card")
##' y <- card$lwage
##' d <- card$educ
##' z <- card$nearc4
##' x <- model.matrix( ~ exper + expersq + black + south + smsa + reg661 + reg662 +
##'                      reg663 + reg664 + reg665+ reg666 + reg667 + reg668 + smsa66,
##'                    data = card)
##'
##' # fit iv model
##' card.fit <- iv_fit(y, d, z, x)
##'
##' # extract coefficients
##' coef(card.fit)
##' coef(card.fit, parm = "fs")
##' coef(card.fit, parm = "rf")
##'
##' # extract confidence intervals
##' confint(card.fit)
##' confint(card.fit, parm = "fs")
##' confint(card.fit, parm = "rf")
##'
##' @export
coef.iv_fit <- function(object, parm = "iv", ...) {
  parm <- match.arg(parm,
                    choices = c("iv",  "fs", "rf"),
                    several.ok = T)
  # which <- which[which %in% names(object$estimates)]
  sapply(object$estimates[parm], function(x) x$estimate)
}

# coef.summary.iv_fit <- function(object, ...) {
#   sapply(object, function(x)x$estimate)
# }



##'
##' @param level coverage level (i.e, 1-alpha). If not provided, it uses the same level as the one provided in \code{\link{iv_fit}}.
##' @returns \code{confint} returns a numeric vector with the confidence interval of interest.
##' @rdname coef.iv_fit
##' @export
confint.iv_fit <- function(object, parm =c("iv", "fs", "rf"), level, ...){
  parm <- match.arg(parm, several.ok = F)
  alpha <- object$pars$alpha
  if (!missing(level)) {
    alpha <- 1 - level
  }

  if (parm == "iv") {
    iv.data        <- with(object$models, iv_model_helper(fs = fs, rf = rf, instrument = "z"))
    args           <- c(iv.data, alpha = alpha)
    ci             <- do.call(what = "ar_confint", args = args)
    attr(ci, "level") <- 1 - alpha
    class(ci) <- c("numeric", "iv.ci")
    return(ci)
  } else {
    ci <- confint(object$models[[parm]], parm = "z", level = 1-alpha)
    rownames(ci) <- parm
    return(ci)
  }

}

##' @param object an object of class \code{\link{iv_fit}}.
##' @param ... arguments passed to other methods
##' @rdname print.iv_fit
##' @export
summary.iv_fit <- function(object,...){
  out <- object$estimates
  class(out) <- "summary.iv_fit"
  return(out)
}





# Printing ----------------------------------------------------------------

##' print and summary methods for \code{iv_fit}
##'
##' The print and summary methods provide verbal descriptions of the results obtained with the function \code{\link{iv_fit}}.
##'
##' @param x an object of class \code{\link{iv_fit}}.
##' @param digits minimal number of significant digits
##' @export
print.iv_fit <- function(x, digits = 3, ...){
  # skip first line
  cat("\n")

  # body
  cat("Instrumental Variable Estimation")
  cat("\n")
  cat("(Anderson-Rubin Approach)")
  cat("\n")
  cat("=============================================")
  printEst(x$estimates$iv, digits = digits)
  cat("=============================================")
  # cat("\n")
  cat("\n")
  cat("See summary for first stage and reduced form.")
  cat("\n")

  # skip last line
  cat("\n")
}


##' @export
print.summary.iv_fit <- function(x, digits = 3, ...){
  # skip first line
  cat("\n")

  # body
  cat("Instrumental Variable Estimation")
  cat("\n")
  cat("(Anderson-Rubin Approach)")
  cat("\n")
  cat("=============================================")
  printEst(x$iv, digits = digits)
  cat("---------------------------------------------")
  printEst(x$fs, digits = digits)
  cat("---------------------------------------------")
  printEst(x$rf, digits = digits)
  cat("=============================================")
  cat("\n")

  # skip last line
  cat("\n")

}

##' @export
print.iv.ci <- function(x, ...) {
  cat("IV conf. interval:\n")
  cat(printCI(x, ...))
  level <- attr(x, "level")
  cat("\nlevel =", level)
}


if.print <- function(x, label, digits, fun = format){
  if (!is.null(x) & !all(is.na(x))) {
    cat("","", paste0(label,":"), fun(x, digits = digits))
    cat("\n")
  }
}


printEst <- function(x, digits = 3, note = T){
  cat("\n")
  cat(paste0(x$name, ":"))
  cat("\n")
  if.print(x$estimate, "Coef. Estimate", digits = digits)
  if.print(x$se, "Standard Error", digits = digits)
  if.print(x$t.value, "t-value", digits = digits)
  if.print(x$p.value, "p-value", digits = digits)
  if.print(x$conf.int, "Conf. Interval", digits = digits, fun = printCI)
  if(note){
    cat("Note:", paste0("H0 = ", format(x$h0, digits = digits), ","),
        paste0("alpha = ", format(x$alpha,digits = digits), ","),
        paste0("df = ", format(x$dof, digits = digits), "."))
    cat("\n")
  }
}


printCI <- function(ci, digits = 3){
  if (length(ci) > 2) {
    out <- paste0("[",paste(round(ci[1:2], digits = digits), collapse = " , "),"]")
    out <- paste(out, "U", paste0("[",paste(round(ci[3:4], digits = digits), collapse = ", "),"]"))
  } else {
    out <- paste0("[",paste(round(ci,digits = digits), collapse = ", "),"]")
  }
  return(out)
}

rho <- function(fs, rf) cor(resid(fs), resid(rf))





