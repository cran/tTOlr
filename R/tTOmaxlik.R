#' Maximum Likelihood Under H1, Given T-statistic
#'
#' Given the \code{t}-statistic for a difference in means,
#' or for a mean difference, and degrees of freedom,
#' determine the maximum likelihood under the alternative
#' H1, and the $t$-statistic for the difference in means
#' that makes the likelihood under H1 a maximum.
#' Additionally, return the likelihood under H0.
#'
#' #' @details The function returns the maximum likelihood estimate
#' of the maximum likelihood on the scale of the $t$-statistic,
#' for the likelihood under the alternative, when the
#' when the $t$-statistic is used as non-centrality parameter.
#' This results in a value for the likelihood ratio that differs
#' from (and is smaller than) the standard likelihood ratio
#' statistic. Additionally, return the likelihoods under H0 and H1.
#'
#' @param t \code{t}-statistic.
#' @param df Degrees of freedom.
#'
#' @return List, with elements
#'  \itemize{
#'  \item maxlik - Maximum likelihood under H1
#'  \item tmax - \code{t}-statistic for difference in means that makes
#'  likelihood a maximum under H1
#'  \item lik0 - Density (one-sided) under H0
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
  maxlik=opt[["maximum"]]
  #
  ## Calculate `fpr` as `(1-prior)/(1-prior+prior*lr)`
  return(c(maxlik=maxlik, tmax=opt[["maximum"]], lik0=lik0))
}
