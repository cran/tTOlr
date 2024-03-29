% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tTOmaxlik.R
\name{tTOmaxlik}
\alias{tTOmaxlik}
\title{Maximum Likelihood Under H1, Given T-statistic}
\usage{
tTOmaxlik(t, df)
}
\arguments{
\item{t}{\code{t}-statistic.}

\item{df}{Degrees of freedom.}
}
\value{
List, with elements
 \itemize{
 \item maxlik - Maximum likelihood under H1
 \item tmax - \code{t}-statistic for difference in means that makes
 likelihood a maximum under H1
 \item lik0 - Density (one-sided) under H0
 }
}
\description{
Given the \code{t}-statistic for a difference in means,
or for a mean difference, and degrees of freedom,
determine the maximum likelihood under the alternative
H1, and the $t$-statistic for the difference in means
that makes the likelihood under H1 a maximum.
Additionally, return the likelihood under H0.
}
\details{
#' @details The function returns the maximum likelihood estimate
of the maximum likelihood on the scale of the $t$-statistic,
for the likelihood under the alternative, when the
when the $t$-statistic is used as non-centrality parameter.
This results in a value for the likelihood ratio that differs
from (and is smaller than) the standard likelihood ratio
statistic. Additionally, return the likelihoods under H0 and H1.
}
\examples{
stats <- tTOmaxlik(t=2, df=5)
## Likelihood ratio, 1-sided test and 2-sided test, p=0.05
tvals1 <- qt(0.05, df=c(2,5,20), lower.tail=FALSE)
tvals2 <- qt(0.025, df=c(2,5,20), lower.tail=FALSE)
likrat1 <- likrat2 <- numeric(3)
for(i in 1:3){
stats1 <- tTOmaxlik(t=tvals1[i], df=c(2,5,20)[i])
likrat1[i] <- stats1[['maxlik']]/stats1[['lik0']]
stats2 <- tTOmaxlik(t=tvals2[i], df=c(2,5,20)[i])
likrat2[i] <- stats2[['maxlik']]/(2*stats2[['lik0']])
# NB: 2*stats2[['lik0']] in denominator.
}
likrat <- rbind('One-sided'=likrat1, 'Two-sided'=likrat2)
colnames(likrat) <- paste0('df=',c(2,5,20))
likrat

}
\references{
van Aubel, A; Gawronski, W (2003).
Analytic properties of noncentral distributions.
Applied Mathematics and Computation. 141: 3–12.
doi:10.1016/S0096-3003(02)00316-8.
}
