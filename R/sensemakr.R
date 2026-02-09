
##' Sensitivity Analysis of Instrumental Variable Estimates
##' @description
##' This function performs sensitivity analysis of instrumental variable estimates, as discussed in Cinelli and Hazlett (2025).
##' The main input is an object of class \code{\link{iv_fit}}. It returns an object of class \code{iv.sensemakr} with several pre-computed sensitivity statistics for reporting. After running \code{sensemakr} you may directly use the \code{plot}, \code{print} and \code{summary} methods in the returned object.
##'
##' @returns An object of class \code{iv.sensemakr}, containing:
##' \describe{
##'  \item{ \code{pars} }{A \code{list} with the general parameters used when calling sensemakr.}
##'  \item{ \code{unadjusted} }{A \code{list} with the original, unadjusted results.}
##'  \item{ \code{sensitivity_stats} }{A \code{list} with the sensitivity statistics of the IV, First-Stage, and Reduced-Form regressions.}
##'  \item{ \code{bounds} }{A \code{list} with bounds on the strength of latent variables if they were "k times" as strong as the benchmark covariates.}
##'  }
##'
##' @export
sensemakr <- sensemakr::sensemakr

##' @param model a model created with the function \code{\link{iv_fit}}.
##' @param benchmark_covariates  character vector of the names of covariates that will be used to bound the plausible strength of the latent variables.
##' @param kz numeric vector. Parameterizes how many times stronger the latent variables are related to the instrument in comparison to the observed benchmark covariates.
##' Default value is \code{1} (latent variable is as strong as benchmark covariate).
##' @param ky numeric vector. Parameterizes how many times stronger the latent variables are related to the (pot.) outcome in comparison to the observed benchmark covariates.
##' @param kd numeric vector. Parameterizes how many times stronger the latent variables are related to the treatment in comparison to the observed benchmark covariates. Default value is the same as \code{kz}.
##' @param r2zw.x (optional) hypothetical partial R2 of latent variables
##' W with the instrument Z, given observed covariates X.
##' @param r2y0w.zx (optional) hypothetical partial R2 of latent variables W with the (pot.) outcome Y(0) given Z and X. Default is the same as \code{r2zw.x}.
##' @param bound_label label to bounds provided manually in \code{r2zw.x} and \code{r2y0w.zx}.
##' @param q 	 percent change of the effect estimate that would be deemed problematic. Default is 1, which means a reduction of 100\% of the current effect estimate (bring estimate to zero).
##' @param alpha significance level.
##' @param min should we consider biases as large or larger than a certain amount? Default is \code{TRUE}.
##' @param ... arguments passed to other methods.
##' @importFrom stats as.formula coef confint cor formula lm na.omit optimize qf quantile resid setNames var
##' @importFrom grDevices contourLines
##' @importFrom graphics contour par points text
##' @importFrom utils modifyList
##' @exportS3Method sensemakr::sensemakr iv_fit
##' @examples
##'
##' # loads package
##' library(iv.sensemakr)
##'
##' # loads dataset
##' data("card")
##'
##' # prepares data
##' y <- card$lwage  # outcome
##' d <- card$educ   # treatment
##' z <- card$nearc4 # instrument
##' x <- model.matrix( ~ exper + expersq + black + south + smsa + reg661 + reg662 +
##'                      reg663 + reg664 + reg665+ reg666 + reg667 + reg668 + smsa66,
##'                    data = card) # covariates
##' # fits IV model
##' card.fit <- iv_fit(y,d,z,x)
##'
##' # see results
##' card.fit
##'
##' # runs sensitivity analysis
##' card.sens <- sensemakr(card.fit, benchmark_covariates = c("black", "smsa"))
##'
##' # see results
##' card.sens
##'
##' # sensitivity contour plot
##' plot(card.sens, lim = 0.09)
##' @references
##'
##' Cinelli, C. and Hazlett, C. (2025), "An Omitted Variable Bias Framework for Sensitivity Analysis of Instrumental Variables." Biometrika. \doi{10.1093/biomet/asaf004}
##' @rdname sensemakr
sensemakr.iv_fit <- function(model,
                             benchmark_covariates = NULL,
                             kz = 1,
                             ky = kz,
                             kd = kz,
                             r2zw.x = NULL,
                             r2y0w.zx = r2zw.x,
                             bound_label = "Manual Bound",
                             q = 1,
                             alpha = 0.05,
                             min = TRUE,
                             ...
                             ){
  out <- list()

  out$pars <- list(q = q,
                   alpha = alpha,
                   h0 = (1-q)*model$estimates$iv$estimate,
                   min = min,
                   benchmark_covariates = benchmark_covariates,
                   kz = kz,
                   ky = ky,
                   kd = kd,
                   r2zw.x = r2zw.x,
                   r2y0w.zx = r2y0w.zx,
                   bound_label = bound_label)

  out$unadjusted <- iv_fit(y = model$data$y,
                                 d = model$data$d,
                                 z = model$data$z,
                                 x = as.matrix(model$data[,-c(1:3)]),
                                 h0 = (1-q)*model$estimates$iv$estimate,
                                 alpha = alpha)
  model <- out$unadjusted

  out$sensitivity_stats$iv <- suppressWarnings(sensitivity_stats(model = model,
                                                                 parm = "iv",
                                                                 q = q,
                                                                 alpha = alpha,
                                                                 min = min))

  out$sensitivity_stats$fs <- suppressWarnings(sensitivity_stats(model = model,
                                                                 parm = "fs",
                                                                 q = 1,
                                                                 alpha = alpha,
                                                                 min = min))

  out$sensitivity_stats$rf <- suppressWarnings(sensitivity_stats(model = model,
                                                                 parm = "rf",
                                                                 q = 1,
                                                                 alpha = alpha,
                                                                 min = min))

  if (!is.null(benchmark_covariates)) {
    out$bounds$iv <- suppressWarnings(ovb4iv_bounds(model, kz = kz, ky = ky, benchmarks = benchmark_covariates,alpha = alpha))
    out$bounds$iv$t.dagger <- sensemakr::adjusted_critical_value(r2dz.x = out$bounds$iv$r2zw.x,
                                                                 r2yz.dx = out$bounds$iv$r2y0w.zx,
                                                                 dof = model$estimates$iv$dof)

    out$bounds$fs <- suppressWarnings(sensemakr::ovb_bounds(model$models$fs, treatment = "z",
                                           kd = kz, ky= kd,
                                           benchmark_covariates = benchmark_covariates,
                                           alpha = alpha))
    out$bounds$fs$t.dagger <- sensemakr::adjusted_critical_value(r2dz.x = out$bounds$fs$r2dz.x,
                                                                 r2yz.dx = out$bounds$fs$r2yz.dx,
                                                                 dof = model$estimates$fs$dof)
    out$bounds$fs$lwr <- with(model$estimates$fs, estimate - out$bounds$fs$t.dagger*se)
    out$bounds$fs$upr <-  with(model$estimates$fs, estimate + out$bounds$fs$t.dagger*se)
    out$bounds$fs <- out$bounds$fs[, c("bound_label", "r2dz.x", "r2yz.dx","lwr", "upr","t.dagger")]
    names(out$bounds$fs) <-  c("bound_label", "r2zw.x", "r2dw.zx", "lwr",  "upr", "t.dagger")

    out$bounds$rf <- suppressWarnings(sensemakr::ovb_bounds(model$models$rf, treatment = "z",
                                           kd = kz, ky= ky, benchmark_covariates = benchmark_covariates,
                                           alpha = alpha))
    out$bounds$rf$t.dagger <- sensemakr::adjusted_critical_value(r2dz.x = out$bounds$rf$r2dz.x,
                                                                 r2yz.dx = out$bounds$rf$r2yz.dx,
                                                                 dof = model$estimates$rf$dof)
    out$bounds$rf$lwr <- with(model$estimates$rf, estimate - out$bounds$rf$t.dagger*se)
    out$bounds$rf$upr <-  with(model$estimates$rf, estimate + out$bounds$rf$t.dagger*se)
    out$bounds$rf <- out$bounds$rf[, c("bound_label", "r2dz.x", "r2yz.dx","lwr", "upr", "t.dagger" )]
    names(out$bounds$rf) <-  c("bound_label", "r2zw.x", "r2yw.zx", "lwr", "upr", "t.dagger")

  }

  if(!is.null(r2zw.x)){
    manual_bounds <- data.frame(bound_label = bound_label, r2zw.x = r2zw.x, r2y0w.zx = r2y0w.zx)
    manual_bounds$lwr <- iv_adjusted_limit(fs = model$models$fs,
                                           rf = model$models$rf,
                                           instrument = "z",
                                           ci.limit = "lwr",
                                           r2zw.x = r2zw.x,
                                           r2y0w.zx = r2y0w.zx,
                                           alpha = alpha)
    manual_bounds$upr <- iv_adjusted_limit(fs = model$models$fs,
                                           rf = model$models$rf,
                                           instrument = "z",
                                           ci.limit = "upr",
                                           r2zw.x = r2zw.x,
                                           r2y0w.zx = r2y0w.zx,
                                           alpha = alpha)
    manual_bounds$t.dagger <- sensemakr::adjusted_critical_value(r2dz.x = r2zw.x,
                                                                 r2yz.dx = r2y0w.zx,
                                                                 dof = model$estimates$rf$dof)
    out$bounds$iv <- rbind(manual_bounds, out$bounds$iv)

  }

  class(out) <- "iv.sensemakr"
  out
}


##' Extract estimates of an \code{iv.sensemakr} object
##'
##' This function extracts the estimate, lower limit, upper limit, t-value, and (extreme) robustness values of an \code{iv.sensemakr} object, created with the function \code{\link{sensemakr}}.
##'
##' @inheritParams coef.iv_fit
##' @export
coef.iv.sensemakr <- function(object, parm = "iv", ...){

  parm <- match.arg(parm,
                    several.ok = T,
                    choices = c("iv", "fs", "rf"))

  do.call("rbind", object$sensitivity_stats)[parm, ]
}



##' @export
print.iv.sensemakr <- function(x,digits = 3,...){
  cat("\n")
  cat("Sensitivity Analysis for Instrumental Variables\n")
  cat("(Anderson-Rubin Approach)")
  cat("\n")
  cat("=============================================================")
  printEst(x$unadjusted$estimates$iv, digits = digits, note = F)
  # cat("-------------------------------------------------------------")
  printSenStats(x$sensitivity_stats$iv, digits = digits, note = F)
  # cat("-------------------------------------------------------------")
  if(!is.null(x$bounds$iv)){
    cat("Bounds on Omitted Variable Bias:")
    cat("\n")
    print(setNames(x$bounds$iv, c("Bound Label", "R2zw.x", "R2y0w.zx", "Lower CI", "Upper CI", "Crit. Thr.")), digits = digits, row.names = F)
    cat("\n")
  }
  cat("Note:",
      paste0("H0 = ", format(x$pars$h0, digits = digits), ","),
      paste0("q", ifelse(x$pars$min, " >= ", " = "), format(x$pars$q, digits = digits), ","),
      paste0("alpha = ", format(x$pars$alpha, digits = digits), ","),
      paste0("df = ", format(x$unadjusted$estimates$iv$dof,digits = digits), "."))
  cat("\n")
  cat("=============================================================")
  cat("\n")
  cat("See summary for first stage and reduced form.")
  cat("\n")

  # skip last line
  cat("\n")
}

#' Sensitivity analysis print and summary methods for \code{iv.sensemakr}
#'
#' @description
#' The \code{print} and \code{summary} methods provide verbal descriptions of the sensitivity analysis results
#' obtained with the function \code{\link{sensemakr}}.
#'
#' @param ... arguments passed to other methods.
#' @param object an object of class \code{\link{sensemakr}}.
#' @param x an object of class \code{\link{sensemakr}}.
#' @param digits minimal number of \emph{significant} digits.
#' @export
#' @rdname print.iv.sensemakr
summary.iv.sensemakr <- function(object, ...){
  class(object) <- "summary.iv.sensemakr"
  object
}

##' @export
##' @rdname print.iv.sensemakr
print.summary.iv.sensemakr <- function(x, digits = 3, ...){
  cat("\n")
  cat("Sensitivity Analysis for Instrumental Variables\n")
  cat("(Anderson-Rubin Approach)")
  cat("\n")
  cat("=============================================================")
  printEst(x$unadjusted$estimates$iv, digits = digits, note = F)
  # cat("-------------------------------------------------------------")
  printSenStats(x$sensitivity_stats$iv, digits = digits, note = F)
  if(!is.null(x$bounds$iv)){
    cat("Bounds on Omitted Variable Bias:")
    cat("\n")
    print(setNames(x$bounds$iv, c("Bound Label", "R2zw.x", "r2y0w.zx", "Lower CI", "Upper CI", "Crit. Thr.")), digits = digits, row.names = F)
    cat("\n")
  }
  cat("Note:",
      paste0("H0 = ", format(x$pars$h0, digits = digits), ","),
      paste0("q", ifelse(x$pars$min, " >= ", " = "), format(x$pars$q, digits = digits), ","),
      paste0("alpha = ", format(x$pars$alpha, digits = digits), ","),
      paste0("df = ", format(x$unadjusted$estimates$iv$dof,digits = digits), "."))
  cat("\n")
  cat("-------------------------------------------------------------")
  printEst(x$unadjusted$estimates$fs, digits = digits, note = F)
  # cat("-------------------------------------------------------------")
  printSenStats(x$sensitivity_stats$fs, digits = digits, note = F)
  if(!is.null(x$bounds$fs)){
    cat("Bounds on Omitted Variable Bias:")
    cat("\n")
    print(setNames(x$bounds$fs, c("Bound Label", "R2zw.x", "R2dw.zx", "Lower CI", "Upper CI", "Crit. Thr.")), digits = digits, row.names = F)
    cat("\n")
  }
  cat("Note:",
      paste0("H0 = 0,"),
      paste0("q = 1,"),
      paste0("alpha = ", format(x$pars$alpha, digits = digits), ","),
      paste0("df = ", format(x$unadjusted$estimates$fs$dof,digits = digits), "."))
  cat("\n")
  cat("-------------------------------------------------------------")
  printEst(x$unadjusted$estimates$rf, digits = digits, note = F)
  # cat("-------------------------------------------------------------")
  printSenStats(x$sensitivity_stats$rf, digits = digits, note = F)
  if(!is.null(x$bounds$rf)){
    cat("Bounds on Omitted Variable Bias:")
    cat("\n")
    print(setNames(x$bounds$rf, c("Bound Label", "R2zw.x", "R2yw.zx", "Lower CI", "Upper CI", "Crit. Thr.")), digits = digits,row.names = F)
    cat("\n")
  }
  cat("Note:",
      paste0("H0 = 0,"),
      paste0("q = 1,"),
      paste0("alpha = ", format(x$pars$alpha, digits = digits), ","),
      paste0("df = ", format(x$unadjusted$estimates$fs$dof,digits = digits), "."))
  cat("\n")
  cat("=============================================================")
  cat("\n")
  cat("\n")
}

# ovb_minimal_reporting -------------------------------------------------------

##' Minimal sensitivity reporting for IV estimates
##'
##' @description
##' This function produces LaTeX or HTML code for a minimal sensitivity analysis table
##' for instrumental variable estimates, as suggested in Cinelli and Hazlett (2025).
##' For objects of class \code{\link{sensemakr}} (from the \pkg{sensemakr} package), it
##' dispatches to \code{sensemakr::\link[sensemakr]{ovb_minimal_reporting}}.
##'
##' @param x an object of class \code{iv.sensemakr} or \code{sensemakr}.
##' @param digits minimal number of \emph{significant} digits.
##' @param verbose if \code{TRUE}, the function prints the code with \code{\link{cat}}.
##' @param format code format to print: \code{"latex"}, \code{"html"} (requires mathjax), or \code{"pure_html"}.
##' @param ... further arguments passed to the table-building functions.
##'   Optional overrides include \code{outcome_label} and \code{treatment_label}.
##'
##' @return The function returns the LaTeX or HTML code invisibly as a character string and also prints it
##' with \code{\link{cat}} when \code{verbose = TRUE}.
##'
##' @references
##' Cinelli, C. and Hazlett, C. (2025), "An Omitted Variable Bias Framework for Sensitivity Analysis of Instrumental Variables." Biometrika. \doi{10.1093/biomet/asaf004}
##'
##' @examples
##' # loads package
##' library(iv.sensemakr)
##'
##' # loads dataset
##' data("card")
##'
##' # prepares data
##' y <- card$lwage
##' d <- card$educ
##' z <- card$nearc4
##' x <- model.matrix( ~ exper + expersq + black + south + smsa + reg661 + reg662 +
##'                      reg663 + reg664 + reg665+ reg666 + reg667 + reg668 + smsa66,
##'                    data = card)
##'
##' # fits IV model and runs sensitivity analysis
##' card.fit <- iv_fit(y, d, z, x)
##' card.sens <- sensemakr(card.fit, benchmark_covariates = c("black", "smsa"))
##'
##' # latex code
##' ovb_minimal_reporting(card.sens)
##'
##' # html code (pure html, no mathjax needed)
##' ovb_minimal_reporting(card.sens, format = "pure_html")
##'
##' @export
ovb_minimal_reporting <- function(x, digits = 3, verbose = TRUE, format = c("latex", "html", "pure_html"), ...) {
  if (inherits(x, "iv.sensemakr")) {
    format <- match.arg(format)
    fun <- switch(format,
                  latex     = iv_latex_table,
                  html      = iv_html_table,
                  pure_html = iv_html_table_no_mathjax)
    fun(x = x, digits = digits, verbose = verbose, ...)
  } else {
    sensemakr::ovb_minimal_reporting(x = x, digits = digits, verbose = verbose, format = format, ...)
  }
}


# internal table functions -------------------------------------------------

iv_latex_table <- function(x, digits = 3, verbose = TRUE, ...) {
  table_settings <- list(...)

  # labels: default to outcome/treatment names from internal models
  outcome_label   <- if (!is.null(table_settings[["outcome_label"]])) table_settings[["outcome_label"]] else all.vars(formula(x$unadjusted$models$rf))[1]
  treatment_label <- if (!is.null(table_settings[["treatment_label"]])) table_settings[["treatment_label"]] else all.vars(formula(x$unadjusted$models$fs))[1]

  q     <- x$pars$q
  alpha <- x$pars$alpha

  # extract stats
  ss  <- x$sensitivity_stats$iv
  est <- x$unadjusted$estimates$iv

  table_begin <- paste0(
    "\\begin{table}[!h]\n",
    "\\centering\n",
    "\\begin{tabular}{lrrrrrr}\n")

  outcome_header <- sprintf(
    "\\multicolumn{7}{c}{Outcome: \\textit{%s}} \\\\\n",
    outcome_label)

  coeff_header <- paste0(
    "\\hline \\hline \n",
    "Treatment: & Est. & Lower CI & Upper CI & t-value ",
    "& $XRV_{q = ", q, ", \\alpha = ", alpha, "}$ ",
    "& $RV_{q = ", q, ", \\alpha = ", alpha, "}$ ",
    " \\\\ \n",
    "\\hline \n")

  coeff_results <- paste0(
    "\\textit{", treatment_label, "} & ",
    round(ss$estimate, digits), " & ",
    round(ss$lwr, digits), " & ",
    round(ss$upr, digits), " & ",
    round(ss$t.value, digits), " & ",
    100 * round(ss$xrv_qa, digits), "\\% & ",
    100 * round(ss$rv_qa, digits), "\\% \\\\ \n")

  # footnote
  dof <- est$dof
  if (!is.null(x$bounds$iv)) {
    row <- x$bounds$iv[1, , drop = FALSE]
    footnote <- paste0(
      "\\hline \n",
      "df = ", dof, " & & ",
      "\\multicolumn{5}{r}{ ",
      "\\small ",
      "\\textit{Bound (", row$bound_label, ")}: ",
      "$R^2_{Z\\sim W| {\\bf X}}$ = ",
      100 * round(row$r2zw.x, digits), "\\%, ",
      "$R^2_{Y(0)\\sim W| Z, {\\bf X}}$ = ",
      100 * round(row$r2y0w.zx, digits), "\\%",
      "} \\\\\n")
  } else {
    footnote <- paste0(
      "\\hline \n",
      "df = ", dof, " & & ",
      "\\multicolumn{5}{r}{ }")
  }

  tabular_end <- "\\end{tabular}\n"
  caption <- if (!is.null(table_settings[["caption"]])) paste0("\\caption{", table_settings[["caption"]], "} \n") else ""
  label   <- if (!is.null(table_settings[["label"]]))   paste0("\\label{",   table_settings[["label"]],   "} \n") else ""
  table_end <- "\\end{table}"

  table <- paste0(table_begin, outcome_header, coeff_header,
                  coeff_results, footnote, tabular_end,
                  caption, label, table_end)

  if (verbose) cat(table)
  return(invisible(table))
}


iv_html_table <- function(x, digits = 3, verbose = TRUE, ...) {
  table_settings <- list(...)

  outcome_label   <- if (!is.null(table_settings[["outcome_label"]])) table_settings[["outcome_label"]] else all.vars(formula(x$unadjusted$models$rf))[1]
  treatment_label <- if (!is.null(table_settings[["treatment_label"]])) table_settings[["treatment_label"]] else all.vars(formula(x$unadjusted$models$fs))[1]

  q     <- x$pars$q
  alpha <- x$pars$alpha

  ss  <- x$sensitivity_stats$iv
  est <- x$unadjusted$estimates$iv

  table_begin <- "<table>\n"

  coeff_header <- paste0(
    "<thead>\n",
    "<tr>\n",
    '\t<th style="text-align:left;border-bottom: 1px solid transparent;border-top: 1px solid black"> </th>\n',
    '\t<th colspan = 6 style="text-align:center;border-bottom: 1px solid black;border-top: 1px solid black"> Outcome: ',
    outcome_label, '</th>\n',
    "</tr>\n",
    "<tr>\n",
    '\t<th style="text-align:left;border-top: 1px solid black"> Treatment </th>\n',
    '\t<th style="text-align:right;border-top: 1px solid black"> Est. </th>\n',
    '\t<th style="text-align:right;border-top: 1px solid black"> Lower CI </th>\n',
    '\t<th style="text-align:right;border-top: 1px solid black"> Upper CI </th>\n',
    '\t<th style="text-align:right;border-top: 1px solid black"> t-value </th>\n',
    '\t<th style="text-align:right;border-top: 1px solid black"> $XRV_{q = ', q, ', \\alpha = ', alpha, '}$ </th>\n',
    '\t<th style="text-align:right;border-top: 1px solid black"> $RV_{q = ', q, ', \\alpha = ', alpha, '}$ </th>\n',
    "</tr>\n",
    "</thead>\n")

  coeff_results <- paste0(
    "<tbody>\n <tr>\n",
    '\t<td style="text-align:left; border-bottom: 1px solid black"><i>',
    treatment_label, "</i></td>\n",
    '\t<td style="text-align:right;border-bottom: 1px solid black">',
    round(ss$estimate, digits), " </td>\n",
    '\t<td style="text-align:right;border-bottom: 1px solid black">',
    round(ss$lwr, digits), " </td>\n",
    '\t<td style="text-align:right;border-bottom: 1px solid black">',
    round(ss$upr, digits), " </td>\n",
    '\t<td style="text-align:right;border-bottom: 1px solid black">',
    round(ss$t.value, digits), " </td>\n",
    '\t<td style="text-align:right;border-bottom: 1px solid black">',
    100 * round(ss$xrv_qa, digits), "\\% </td>\n",
    '\t<td style="text-align:right;border-bottom: 1px solid black">',
    100 * round(ss$rv_qa, digits), "\\% </td>\n",
    "</tr>\n</tbody>\n")

  table_end <- "</table>"

  dof <- est$dof
  if (!is.null(x$bounds$iv)) {
    row <- x$bounds$iv[1, , drop = FALSE]
    footnote <- paste0(
      '<tr>\n',
      "<td colspan = 7 style='text-align:right;border-top: 1px solid black;border-bottom: 1px solid transparent;font-size:11px'>",
      "Note: df = ", dof, "; ",
      "Bound ( ", row$bound_label, " ):  ",
      "$R^2_{Z\\sim W| {\\bf X}}$ = ",
      100 * round(row$r2zw.x, digits), "\\%, ",
      "$R^2_{Y(0)\\sim W| Z, {\\bf X}}$ = ",
      100 * round(row$r2y0w.zx, digits), "\\%",
      "</td>\n</tr>\n")
  } else {
    footnote <- paste0(
      '<tr>\n',
      "<td colspan = 7 style='text-align:right;border-top: 1px solid black;border-bottom: 1px solid transparent;font-size:11px'>",
      "Note: df = ", dof, "; ",
      "</td>\n</tr>\n")
  }

  table <- paste0(table_begin, coeff_header, coeff_results, footnote, table_end)

  if (verbose) cat(table)
  return(invisible(table))
}


iv_html_table_no_mathjax <- function(x, digits = 3, verbose = TRUE, ...) {
  table_settings <- list(...)

  outcome_label   <- if (!is.null(table_settings[["outcome_label"]])) table_settings[["outcome_label"]] else all.vars(formula(x$unadjusted$models$rf))[1]
  treatment_label <- if (!is.null(table_settings[["treatment_label"]])) table_settings[["treatment_label"]] else all.vars(formula(x$unadjusted$models$fs))[1]

  q     <- x$pars$q
  alpha <- x$pars$alpha

  ss  <- x$sensitivity_stats$iv
  est <- x$unadjusted$estimates$iv

  table_begin <- "<table style='align:center'>\n"

  coeff_header <- paste0(
    "<thead>\n",
    "<tr>\n",
    '\t<th style="text-align:left;border-bottom: 1px solid transparent;border-top: 1px solid black"> </th>\n',
    '\t<th colspan = 6 style="text-align:center;border-bottom: 1px solid black;border-top: 1px solid black"> Outcome: ',
    outcome_label, '</th>\n',
    "</tr>\n",
    "<tr>\n",
    '\t<th style="text-align:left;border-top: 1px solid black"> Treatment </th>\n',
    '\t<th style="text-align:right;border-top: 1px solid black"> Est. </th>\n',
    '\t<th style="text-align:right;border-top: 1px solid black"> Lower CI </th>\n',
    '\t<th style="text-align:right;border-top: 1px solid black"> Upper CI </th>\n',
    '\t<th style="text-align:right;border-top: 1px solid black"> t-value </th>\n',
    '\t<th style="text-align:right;border-top: 1px solid black"> XRV<sub>q = ', q, ', &alpha; = ', alpha, '</sub> </th>\n',
    '\t<th style="text-align:right;border-top: 1px solid black"> RV<sub>q = ', q, ', &alpha; = ', alpha, '</sub> </th>\n',
    "</tr>\n",
    "</thead>\n")

  coeff_results <- paste0(
    "<tbody>\n <tr>\n",
    '\t<td style="text-align:left; border-bottom: 1px solid black"><i>',
    treatment_label, "</i></td>\n",
    '\t<td style="text-align:right;border-bottom: 1px solid black">',
    round(ss$estimate, digits), " </td>\n",
    '\t<td style="text-align:right;border-bottom: 1px solid black">',
    round(ss$lwr, digits), " </td>\n",
    '\t<td style="text-align:right;border-bottom: 1px solid black">',
    round(ss$upr, digits), " </td>\n",
    '\t<td style="text-align:right;border-bottom: 1px solid black">',
    round(ss$t.value, digits), " </td>\n",
    '\t<td style="text-align:right;border-bottom: 1px solid black">',
    100 * round(ss$xrv_qa, digits), "% </td>\n",
    '\t<td style="text-align:right;border-bottom: 1px solid black">',
    100 * round(ss$rv_qa, digits), "% </td>\n",
    "</tr>\n</tbody>\n")

  table_end <- "</table>"

  dof <- est$dof
  if (!is.null(x$bounds$iv)) {
    row <- x$bounds$iv[1, , drop = FALSE]
    footnote <- paste0(
      '<tr>\n',
      "<td colspan = 7 style='text-align:right;border-bottom: 1px solid transparent;font-size:11px'>",
      "Note: df = ", dof, "; ",
      "Bound ( ", row$bound_label, " ):  ",
      "R<sup>2</sup><sub>Z~W|X</sub> = ",
      100 * round(row$r2zw.x, digits),
      "%, R<sup>2</sup><sub>Y(0)~W|Z,X</sub> = ",
      100 * round(row$r2y0w.zx, digits), "%",
      "</td>\n</tr>\n")
  } else {
    footnote <- paste0(
      '<tr>\n',
      "<td colspan = 7 style='text-align:right;border-bottom: 1px solid transparent;font-size:11px'>",
      "Note: df = ", dof, "; ",
      "</td>\n</tr>\n")
  }

  table <- paste0(table_begin, coeff_header, coeff_results, footnote, table_end)

  if (verbose) cat(table)
  return(invisible(table))
}


ovb4iv_bounds <- function(model, kz, ky,
                          benchmarks,
                          alpha = 0.05){
  out <- ovb4iv_partial_r2_bound(model, benchmark = benchmarks, kz = kz, ky= ky, alpha = alpha, value = F)
  out$lwr <- iv_adjusted_limit(fs = model$models$fs,
                               rf = model$models$rf,
                               instrument = "z",
                               ci.limit = "lwr",
                               r2zw.x = out$r2zw.x,
                               r2y0w.zx = out$r2y0w.zx,
                               alpha = alpha)
  out$upr <- iv_adjusted_limit(fs = model$models$fs,
                               rf = model$models$rf,
                               instrument = "z",
                               ci.limit = "upr",
                               r2zw.x = out$r2zw.x,
                               r2y0w.zx = out$r2y0w.zx,
                               alpha = alpha)
  # out$t.dagger <- adjusted_critical_value(r2dz.x = out$bounds$iv$r2zw.x,
  #                                         r2yz.dx = out$bounds$iv$r2y0w.zx,
  #                                         dof = card.fit$estimates$iv$dof,
  #                                         alpha = card.fit$pars$alpha)
  rownames(out) <- NULL
  out
}
