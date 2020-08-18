#' Maximum Likelihood Under H1, Given P-value.
#'
#' Given the \code{t}-statistic for a difference in means,
#' or for a mean difference, and degrees of freedom,
#' determine the maximum likelihood under the alternative
#' H1, and the $t$-statistic for the difference in means
#' that makes the likelihood under H1 a maximum.
#' Additionally, return the likelihood under H0.
#'
#' @details As a minimum, supply a $t$-statistic, either degrees
#'
#' @param t \code{t}-statistic.  If \code{NULL}, this is calculated from
#'  the \code{p}-value.
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
#'  \item lik0 - Likelihood under H0
#'  \item likDelta - Likelihood, given difference delta under H0
#'  \item lrDelta - Likelihood ratio, given difference delta under H0
#'  \item maxlik - Maximum likelihood, under allowed alternatives H1
#'  \item lrmax - Maximum of likelihood ratio,
#'  under allowed alternatives H1
#'  \item tmax - \code{t}-statistic for difference in means
#'  that makes likelihood under H1 a maximum
#'  }
#'
#' @importFrom stats pt qt
#' @importFrom lattice xyplot
#'
#' @export
#'
#' @examples
#' likStats <- tTOlr(pval=0.02, nsamp=c(9,9), twoSided=TRUE,
#'                 delta=1.4, sd=1.2)
#' print(unlist(likStats),digits=2)
#' likStats <- tTOlr(t=2.58, df=16, nsamp=c(9,9), twoSided=TRUE,
#'                   delta=1.4, sd=1.2)
#' print(unlist(likStats),digits=2)
#' likStats <- tTOlr(pval=0.02, nsamp=9, twoSided=FALSE,
#'                 delta=1.4, sd=1.2)
#' print(unlist(likStats),digits=2)
#' likStats <- tTOlr(t=2.45, df=8, nsamp=9, twoSided=FALSE,
#'                   delta=1.4, sd=1.2)
#' print(unlist(likStats),digits=2)
#'
tTOlr =  function(t=NULL, df=NULL, nsamp=NULL, pval=NULL, delta=NULL,
                  sd=1, twoSided=TRUE, showMax=TRUE)
{
  if(!is.null(nsamp)&!is.null(df))
    if(df!=sum(nsamp)-length(nsamp)){
    print(paste0('Sample size(s) ',paste(nsamp,collapse=','),
                 ' require df=',sum(nsamp)-length(nsamp), ', not df=',df))
    return()
  }
  if(is.null(df))df <- sum(nsamp)-length(nsamp)
  if(!is.null(nsamp)){
    if(length(nsamp)==1)sdiff <- sqrt(sd^2/nsamp) else
      if(length(nsamp)==2)sdiff=sqrt((1/nsamp[1]+1/nsamp[2])*sd^2)
  } else sdiff <- NULL
  #under H0, use central t distribution
  if(is.null(pval))pval <- pt(-abs(t),df=df)*(twoSided+1)
  if(is.null(t))
    if(twoSided)t<-qt((1-pval/2),df,ncp=0) else
      t<-qt((1-pval),df,ncp=0)
    lik0=dt(t,df,0)
    if(twoSided)ref0 <- 2*lik0 else ref0 <- lik0
    if(!is.null(delta)&!is.null(sdiff))
      likDelta <- dt(t, df, ncp=delta/sdiff) else likDelta <- NULL
    ##
    if(showMax)maxStats <- tTOlr::tTOmaxlik(t, df) else
      maxStats <- NULL
    ## maxlik, lrmax, tmax
    ## Calculate `fpr` as `(1-prior)/(1-prior+prior*lr)`
    likStats <- c(list(t=t, pval=pval, lik0=lik0),
                  if(!is.null(likDelta))list(likDelta=likDelta, lrDelta=likDelta/ref0)else NULL,
                  if(!is.null(maxStats))with(maxStats,
                                             list(maxlik=maxlik, lrmax=maxlik/ref0, tmax=tmax))else NULL)
    likStats
}
