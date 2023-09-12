#' Maximum Likelihood ratio for H1 versus H0, given t-statistic or p-value
#'
#' Given the \code{t}-statistic for a difference in means,
#' or for a mean difference, and degrees of freedom,
#' determine the maximum likelihood under the alternative
#' H1, and the $t$-statistic for the difference in means
#' that makes the likelihood under H1 a maximum. Also
#' available is the likelihood that corresponds to  a
#' particular value of a particular effect size (mean
#' divided by standard deviation) \code{delta}.
#'
#' @details The function returns the maximum likelihood estimate
#' of the maximum likelihood on the scale of the $t$-statistic,
#' for the likelihood under the alternative, when the
#' when the $t$-statistic is used as non-centrality parameter.
#' This results in a value for the likelihood ratio that differs
#' from (and is smaller than) the standard likelihood ratio
#' statistic. Additionally, return the likelihoods under H0 and H1.
#'
#' @note The likelihood estimate for H1 versus H0 is unchanged if
#' the roles of H0 and H1 are reversed.
#'
#' @param t \code{t}-statistic.  If \code{NULL}, this is calculated
#' from the \code{p}-value.
#' @param df Degrees of freedom.
#' @param nsamp Sample size.  For a two-sample test, this should be
#' a vector of length 2.
#' @param pval \code{p}-value.  If \code{NULL}, this is calculated
#' from the \code{t}-statistic and degrees of freedom.
#' @param delta If not \code{NULL}, this specifies the $t$-statistic
#'  for the difference from H0 that is of interest, allowing the
#'  calculation of the corresponding likelihood and likelihood ratio.
#' @param sd Standard deviation.
#' @param twoSided Set either to \code{TRUE} for a two-sided test,
#'  or \code{FALSE} for a one-sided test.
#' @param showMax Set to \code{TRUE} if the maximum of the likelihood
#'  and the likelihood ratio is required.
#' @return List, with elements
#'  \itemize{
#'  \item t - \code{t}-statistic
#'  \item df - Degrees of freedom
#'  \item pval - P-value
#'  \item likDelta - Likelihood, given difference delta under H0
#'  \item lrDelta - Likelihood ratio, given difference delta under H0
#'  \item maxlik - Maximum likelihood, under allowed alternatives H1
#'  \item lrmax - Maximum likelihood under H1, on the scale of the
#'  $t$-statistic
#'  \item tmax - \code{t}-statistic for difference in means
#'  that makes likelihood under H1 a maximum
#'  }
#'
#' @importFrom stats pt qt
#' @importFrom lattice xyplot
#' @importFrom latticeExtra axis.grid
#'
#' @export
#'
#' @examples
#' likStats <- tTOlr(pval=0.02, nsamp=c(9,9), twoSided=TRUE,
#'                 delta=1.4, sd=1.2)
#' print(likStats,digits=2)
#' likStats <- tTOlr(t=2.58, df=16, nsamp=c(9,9), twoSided=TRUE,
#'                   delta=1.4, sd=1.2)
#' print(likStats,digits=2)
#' likStats <- tTOlr(pval=0.02, nsamp=9, twoSided=FALSE,
#'                 delta=1.4, sd=1.2)
#' print(likStats,digits=2)
#' likStats <- tTOlr(t=2.45, df=8, nsamp=9, twoSided=FALSE,
#'                   delta=1.4, sd=1.2)
#' print(likStats,digits=2)
#'
tTOlr =  function (t = NULL, df = NULL, nsamp = NULL, pval = NULL, delta = NULL,
                   sd = 1, twoSided = TRUE, showMax = TRUE)
{
  if (!is.null(nsamp) & !is.null(df))
    if (df != sum(nsamp) - length(nsamp)) {
      print(paste0("Sample size(s) ", paste(nsamp, collapse = ","),
                   " require df=", sum(nsamp) - length(nsamp), ", not df=",
                   df))
      return()
    }
  if (is.null(df))
    df <- sum(nsamp) - length(nsamp)
  if (!is.null(nsamp)) {
    if (length(nsamp) == 1)
      sdiff <- sqrt(sd^2/nsamp)
    else if (length(nsamp) == 2)
      sdiff = sqrt((1/nsamp[1] + 1/nsamp[2]) * sd^2)
  }
  else sdiff <- NULL
  if (is.null(pval))
    pval <- pt(-abs(t), df = df) * (twoSided + 1)
  if (is.null(t))
    if (twoSided)
      t <- qt((1 - pval/2), df, ncp = 0)
  else t <- qt((1 - pval), df, ncp = 0)
  lik0 = dt(t, df, 0)
  if (twoSided)
    ref0 <- 2 * lik0
  else ref0 <- lik0
  if (!is.null(delta) & !is.null(sdiff))
    likDelta <- dt(t, df, ncp = delta/sdiff)
  else likDelta <- NULL
  if (showMax)
    maxStats <- tTOlr::tTOmaxlik(t, df)
  else maxStats <- NULL
  likStats <- c(c(t = t, pval = pval, lik0 = lik0),
                if (!is.null(likDelta)) c(likDelta = likDelta,  lrDelta = likDelta/ref0) else NULL,
                if (!is.null(maxStats)) with(as.list(maxStats),
                c(maxlik = maxlik, lrmax = maxlik/ref0, tmax = tmax)) else NULL)
  likStats
}
