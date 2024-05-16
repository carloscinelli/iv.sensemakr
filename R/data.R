#' National Longitudinal Survey of Young Men (NLSYN)
#'
#' @description
#'
#' Data used in Card (1995). Consists of a sample of 3,010 individuals from the National Longitudinal Survey of Young Men (NLSYM).
#'
#'The treatment is \code{educ}, the outcome is \code{lwage} and the instrument is \code{nearc4}.
#'
#' @docType data
#'
#' @usage data('card')
#'
#' @format A data.frame with 3010 observations and 34 variables:
#' \itemize{
#'  \item \strong{id:} person identifier
#'  \item \strong{nearc2:} =1 if near 2 yr college, 1966
#'  \item \strong{nearc4:} =1 if near 4 yr college, 1966
#'  \item \strong{educ:} years of schooling, 1976
#'  \item \strong{age:} in years
#'  \item \strong{fatheduc:} father's schooling
#'  \item \strong{motheduc:} mother's schooling
#'  \item \strong{weight:} NLS sampling weight, 1976
#'  \item \strong{momdad14:} =1 if live with mom, dad at 14
#'  \item \strong{sinmom14:} =1 if with single mom at 14
#'  \item \strong{step14:} =1 if with step parent at 14
#'  \item \strong{reg661:} =1 for region 1, 1966
#'  \item \strong{reg662:} =1 for region 2, 1966
#'  \item \strong{reg663:} =1 for region 3, 1966
#'  \item \strong{reg664:} =1 for region 4, 1966
#'  \item \strong{reg665:} =1 for region 5, 1966
#'  \item \strong{reg666:} =1 for region 6, 1966
#'  \item \strong{reg667:} =1 for region 7, 1966
#'  \item \strong{reg668:} =1 for region 8, 1966
#'  \item \strong{reg669:} =1 for region 9, 1966
#'  \item \strong{south66:} =1 if in south in 1966
#'  \item \strong{black:} =1 if black
#'  \item \strong{smsa:} =1 in in SMSA, 1976
#'  \item \strong{south:} =1 if in south, 1976
#'  \item \strong{smsa66:} =1 if in SMSA, 1966
#'  \item \strong{wage:} hourly wage in cents, 1976
#'  \item \strong{enroll:} =1 if enrolled in school, 1976
#'  \item \strong{KWW:} knowledge world of work score
#'  \item \strong{IQ:} IQ score
#'  \item \strong{married:} =1 if married, 1976
#'  \item \strong{libcrd14:} =1 if lib. card in home at 14
#'  \item \strong{exper:} age - educ - 6
#'  \item \strong{lwage:} log(wage)
#'  \item \strong{expersq:} exper^2
#' }
#'
#' @references
#'
#' Card, D. "Using Geographic Variation in College Proximity to Estimate the Return to Schooling". In L.N. Christofides, E.K. Grant, and R. Swidinsky, editors, Aspects of Labor Market Behaviour: Essays in Honour of John Vanderkamp. Toronto: University of Toronto Press, 1995.
#'
#' @examples
#' data('card')
#' head(card)
"card"
