#' Maximum Likelihood Under H1, Given T-statistic
#'
#' Given the \code{t}-statistic for a difference in means,
#' or for a mean difference, and degrees of freedom,
#' determine the maximum likelihood under the alternative
#' H1, and the $t$-statistic for the difference in means
#' that makes the likelihood under H1 a maximum.
#' Additionally, return the likelihood under H0.
#'
#' @details Because the \code{t}-distribution mean under H1
#' is a random variable, one has a non-central \code{t},
#' and the mode (which gives the maximum) differs somewhat
#' from the mean.
#'
#' @param t \code{t}-statistic.
#' @param df Degrees of freedom.
#'
#' @return List, with elements
#'  \itemize{
#'  \item maxlik - Maximum likelihood under H1
#'  \item tmax - \code{t}-statistic for difference in means that makes
#'  likelihood a maximum
#'  \item lik0 - Likelihood under H0
#'  }
#'
#' @references  van Aubel, A; Gawronski, W (2003).
#' Analytic properties of noncentral distributions.
#' Applied Mathematics and Computation. 141: 3–12.
#' doi:10.1016/S0096-3003(02)00316-8.
#'
#' @importFrom stats dt optimize setNames
#'
#' @export
#'
#' @examples
#' stats <- tTOmaxlik(t=2, df=5)
#' likrat <- stats[['maxlik']]/stats[['lik0']]
#' c("Maximum likelihood ratio"=likrat)
#' ## Likelihood ratio, 1-sided test and 2-sided test, p=0.05
#' tvals1 <- qt(0.05, df=c(2,5,20), lower.tail=FALSE)
#' tvals2 <- qt(0.025, df=c(2,5,20), lower.tail=FALSE)
#' likrat1 <- likrat2 <- numeric(3)
#' for(i in 1:3){
#' stats1 <- tTOmaxlik(t=tvals1[i], df=c(2,5,20)[i])
#' likrat1[i] <- stats1[['maxlik']]/stats1[['lik0']]
#' stats2 <- tTOmaxlik(t=tvals2[i], df=c(2,5,20)[i])
#' likrat2[i] <- stats2[['maxlik']]/(2*stats2[['lik0']])
#' # NB: 2*stats2[['lik0']] in denominator.
#' }
#' likrat <- rbind('One-sided'=likrat1, 'Two-sided'=likrat2)
#' colnames(likrat) <- paste0('df=',c(2,5,20))
#' likrat
#'
tTOmaxlik =  function(t, df)
{
  #
  # under H1 use non-central t distribution
  # ncp is non-centrality paramater
  lik0 <- dt(t,df)
  int <- t*sqrt(c(df/(df+2.5), df/(df+1)))
  opt <- optimize(f=function(x)dt(x,df,ncp=t),maximum=TRUE, interval=int)
  maxlik=opt[["objective"]]
  #
  ## Calculate `fpr` as `(1-prior)/(1-prior+prior*lr)`
  return(list(maxlik=maxlik, tmax=opt[["maximum"]], lik0=lik0))
}
