## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
library(tTOlr)
knitr::opts_chunk$set(echo = TRUE, comment=NA)

## ----pkgs, echo=FALSE---------------------------------------------------------
suppressPackageStartupMessages(library(lattice))

## ----wear---------------------------------------------------------------------
wear <- with(MASS::shoes, rbind(A,B,d=B-A))
colnames(wear) <- rep("",10)
wear

## ----shoes, echo=FALSE--------------------------------------------------------
wear <- with(MASS::shoes, rbind(A, B, d=B-A))
library(magrittr)
shoeStats <- wear['d',] %>% 
  {list(Mean=mean(.), SD=sd(.), n=length(.))} %$%
  {c(., SEM=SD/sqrt(n), t=Mean*sqrt(n)/SD)} %$%
  {c(., pval=2*pt(t, df=n-1, lower.tail=FALSE), df=n-1)}
print(prettyNum(unlist(shoeStats), digits=2, dropTrailing=TRUE),
      quote=FALSE)

## ----shoelik, echo=FALSE------------------------------------------------------
maxlrShoes <- with(shoeStats, 
               tTOmaxlik(t=t, df=df)[['maxlik']]/(2*dt(t,df)))

## ----maxlrGPH, echo=FALSE-----------------------------------------------------
maxlrGPH <- function(sdev=1, df=9, pval=0.05, offset=0.5, twoSided=TRUE, pretitle=""){
#
# under H1 use non-central t distribution
# ncp is non-centrality paramater
if(twoSided)t<-qt((1-pval/2),df,ncp=0) else
  t<-qt((1-pval),df,ncp=0)
lik0 <- dt(t, df)
maxStats <- tTOmaxlik(t, df)
maxlik <- maxStats[['maxlik']]
tmax <- maxStats[["tmax"]]
titl <- paste0(pretitle,"p = ", signif(pval,2),"; LR = ", signif(maxlik,2), "/(2*",signif(lik0,2),") = ", signif(maxlik/(2*lik0),3))
vals <- pretty(c(-4, 4),50)
vals2 <- pretty(c(-4, 5),54)
dframe <- data.frame(x=c(vals,vals2+t),
                 y=c(dt(vals, df=df),
                     dt(vals2+t, df=df,ncp=t)),
                 hyp=rep(c("H0","H1"),c(length(vals),length(vals2))))
#
# under H1 use non-central t distribution
# ncp is non-centrality paramater
atx <- seq(from=-4, to=7, by=1)
# labx <- as.expression( sapply(atx, 
#             function(x) bquote(italic(.(x)*italic(s)))) )
labx <- paste(atx)
col <- trellis.par.get()$superpose.symbol$col 
xyplot(y ~ x, groups=hyp, data=dframe, type='l',
       xlab=expression(italic(t)*"-statistic for difference from NULL (H0)"),
       ylab="Probability density",
       xlim=range(dframe$x), ylim=c(0,max(dframe$y)*1.04),
       scales=list(x=list(at=atx,labels=labx)),
       par.settings=list(clip=list(panel='off',                             layout.heights=list(axis.xlab.padding=0.5))),
       panel=function(x,y,...){
         panel.xyplot(x,y,...)
         panel.abline(v=t, col='gray')
         xlim <- current.panel.limits()$xlim
         ylim <- current.panel.limits()$ylim
         if(twoSided){
          panel.axis("left",at=lik0, outside=TRUE, line.col=col[1],
                    labels='H0', text.col=col[1], text.cex=1.0)
          # panel.arrows(x0=xlim[1], y0=lik0, x1=-t, line.alpha=0.25, 
          #             y1=lik0, length=0.06, col=col[1])
         panel.abline(v=-t, col='gray',lty=1)
         xx <- with(subset(dframe,hyp=='H0'),
                    c(min(x),x[x< -t],-t,-t))
         yy <- with(subset(dframe,hyp=='H0'), c(0,y[x< -t],lik0,0))
         panel.polygon(xx,yy,...,col=col[1],alpha=0.4)
         }
                  xx <- with(subset(dframe,hyp=='H0'),
                    c(t,t,x[x>t],max(x)))
         yy <- with(subset(dframe,hyp=='H0'), c(0,lik0,y[x>t],0))
         panel.polygon(xx,yy,...,col=col[1],alpha=0.4)
         panel.lines(x=c(tmax,0.5*(tmax+xlim[2])),y=c(maxlik,maxlik),
                     lty=2, alpha=0.75)
         panel.axis("right",at=maxlik,labels='H1', outside=TRUE,
                    line.col=col[2],text.col=col[2],text.cex=1.0)
         panel.axis("right",at=maxlik,text.col=col[2], text.cex=1.0,
                    labels=round(maxlik,2),line.col=col[2])
         panel.axis("right",at=lik0, line.col=col[1], labels='H0',
                    outside=TRUE, text.cex=1.0, text.col=col[1])
         panel.axis("right",at=lik0,text.col=col[1],text.cex=1.0,
                    labels=round(lik0,4), line.col=col[1])
       },
       main=list(titl, just="left", font=1, x=0, y=-0.6)
)
}

## ----cap1, echo=FALSE---------------------------------------------------------
cap1 <- "Panel A shows density curves for NULL and for 
the alternative, for a two-sided test with $t$ = 3.35, 
on 9 degrees of freedom, for the comparison of shoe
materials (B versus A) in the dataset MASS::shoes.
Vertical lines are placed at the positions that give 
the $p$-value.  Panel B shows the normal probability 
plot for the \\texttt{B-A} differences in the dataset.
"

## ----DOshoeDens, fig.width=6, fig.asp=0.6, fig.align='center', fig.cap=cap1, out.width='80%', echo=FALSE----
maxlrGPH(pval=shoeStats[['pval']], twoSided=TRUE)

## ----sleep--------------------------------------------------------------------
sleep2 <-with(sleep, rbind(Drug1=extra[group==1], Drug2=extra[group==2],
              d=extra[group==2]-extra[group==1]))
colnames(sleep2) <- rep("",10)
sleep2 
t <- t.test(-extra ~ group, data = sleep, paired=TRUE)

## ----sleeplik, echo=FALSE-----------------------------------------------------
maxlrSleep <- with(t, tTOmaxlik(statistic, df=parameter)[['maxlik']]/
                     (2*dt(statistic,df=parameter)))

## ----cutoff, echo=FALSE-------------------------------------------------------
t <- t.test(-extra ~ group, data = sleep, mu=0.75, 
            paired=TRUE, alternative = 'greater')

## ----cap2, echo=FALSE---------------------------------------------------------
cap2 <- "Panel A shows density curves for NULL and for 
the alternative, for a one-sided test with $t$ = 
2.01 on 9 degrees of freedom.  This is the $t$-statistic 
for the data on the effect of soporofic drugs when 
differences are \\texttt{B-A-0.75}, i.e., interest is in
the strength of evidence that differences are at least
0.75 hours. A vertical line is placed at the position that 
gives the $p$-value.  Panel B shows the normal probability 
plot for the differences.
"

## ----DOmaxlr, fig.width=7, fig.height=3.25, fig.pos='ht', fig.cap=cap2, out.width='100%', echo=FALSE----
gph1 <- maxlrGPH(pval=t$p.value, pretitle=" A: ")
print(gph1, position=c(0,0,0.55,1))
gph2 <- qqmath(sleep2['d',]-0.75, ylab="B-A-0.75",
               main=list(' B: Normal probability plot', 
                       just='left', font=1, x=0, y=-0.6))
print(gph2, position=c(0.55,0,1,1),new=FALSE)

## ----cutoff, echo=TRUE, eval=FALSE--------------------------------------------
#  t <- t.test(-extra ~ group, data = sleep, mu=0.75,
#              paired=TRUE, alternative = 'greater')

## ----cap3, eval=TRUE, echo=FALSE----------------------------------------------
cap3 <- 'Ratio of the maximum likelihood under the alternative
to the likelihood under the NULL, for three different 
choices of $p$-value, for a range of sample sizes, and for a 
range of degrees of freedom.'

## ----pTOmaxlrGPH, echo=FALSE, fig.width=7, fig.asp=0.6, fig.cap=cap3, fig.pos='ht', out.width='100%'----
lrP <- data.frame(df=rep(seq(from=5, to=50,by=5),3), 
                  pval=rep(c(.05,.01,.001),c(10,10,10)),
                  tstat=numeric(30), lr=numeric(30))
for(i in 1:30){
  t <- qt(1-lrP$pval[i]/2,df=lrP$df[i])
  lik0 <- dt(t, df = lrP$df[i])
  lrP[i,c('tstat','lr')] <- c(t, tTOmaxlik(t = t, df = lrP$df[i])[["maxlik"]]/(2*lik0))  ## NB: Two-sided
}
aty <- list(c(2, 4,10,25,50,100,200),c(2,2.5,3,4,5,6))
laby <- list(paste(aty[[1]]),paste(aty[[2]]))
gph <- xyplot(lr+tstat ~ df, groups=pval, data=lrP, auto.key=list(columns=3), ylab="",
                strip=strip.custom(factor.levels=c("Maximum likelihood ratio","t-statistic")),
                scales=list(y=list(relation='free', log=TRUE,
                                   at=aty, labels=laby)), 
                par.settings=lattice::simpleTheme(pch=16))
update(gph,axis=latticeExtra::axis.grid)

## ----cap4, eval=TRUE, echo=FALSE----------------------------------------------
cap4 <- 'False positive risk, for three different 
choices of $p$-value, for a range of sample sizes, and for a 
range of degrees of freedom.'

## ----pTOfprGPH, echo=FALSE, fig.width=7, fig.asp=0.6, fig.cap=cap4, fig.pos='ht', out.width='100%'----
lrP[,"fpr0.1"] <- (1-0.1)/(1-0.1+0.1*lrP[,"lr"])
lrP[,"fpr0.5"] <- (1-0.5)/(1-0.5+0.5*lrP[,"lr"])
aty <- list(c(0.025,0.05,0.1,0.2,0.4),c(0.004,0.01, 0.04,0.1,0.4))
laby <- list(paste(aty[[1]]),paste(aty[[2]]))
gph <- xyplot(fpr0.1+fpr0.5 ~ df, groups=pval, data=lrP, auto.key=list(columns=3), 
                xlab="Degrees of freedom",
                ylab="False positive risk",
                strip=strip.custom(factor.levels=c("Prior = 0.1","Prior = 0.5")),
                scales=list(y=list(relation='free', log=TRUE,
                                   at=aty, labels=laby)), 
                par.settings=lattice::simpleTheme(pch=16))
update(gph,axis=latticeExtra::axis.grid)

## ----explain-power, echo=FALSE------------------------------------------------
cfdenspwr <- function(sdev=1, nsamp=19, sd=1.5, pval=0.05, power=0.8, sig.level=0.05, type="two.sample", alternative = 'one.sided'){
  delta <- power.t.test(n=nsamp, sig.level=sig.level, 
                        type=type,
                        alternative = alternative,
                        power=power)[['delta']]
sed <- sdev*sqrt(2/nsamp)
df <- 2*(nsamp-1)
vals <- pretty(c(-3.25, 3.25),50)
dframe <- data.frame(x=c(vals,vals+delta/sed),
                 y=c(dt(vals, df=(nsamp-1)*2),
                     dt(vals+delta/sed, df=(nsamp-1)*2,ncp=delta/sed)),
                 hyp=rep(c('H0','H1'),rep(length(vals),2)))
tcrit<-qt((1-sig.level),df,ncp=0)
lik0=dt(tcrit,df,0)
#
# under H1 use non-central t distribution
# ncp is non-centrality paramater
lik1=dt(tcrit,df,ncp=delta/sed)
atx <- seq(from=-3, to=6, by=1)
# labx <- as.expression( sapply(atx, 
#             function(x) bquote(italic(.(x)*italic(s)))) )
labx <- paste(atx)
col <- trellis.par.get()$superpose.symbol$col 
xyplot(y ~ x, groups=hyp, data=dframe, type='l',
       xlab=expression(phantom(0)),
       sub=expression(italic(t)*'-statistic'),
       ylab='Probability density',
       xlim=range(dframe$x), ylim=c(0,max(dframe$y)*1.08),
       scales=list(x=list(at=atx,labels=labx)),
       par.settings=list(clip=list(panel='off',
                                  layout.heights=list(axis.xlab.padding=0.5))),
       panel=function(x,y,...){
         panel.xyplot(x,y,...)
         panel.abline(v=tcrit, col='gray')
         xx <- with(subset(dframe,hyp=='H1'),
                    c(tcrit,tcrit,x[x>tcrit],max(x)))
         yy <- with(subset(dframe,hyp=='H1'), c(0,lik1,y[x>tcrit],0)) 
         panel.polygon(xx,yy,...,col=col[2], alpha=0.4)
         panel.text(x=2.82, y=0.16, 
                    expression(atop('Area = 0.8','(= Power)')))
         maxy <- dt(0, df)
         ncp <- delta/sed
         panel.arrows(x0=0, y0=maxy, x1=ncp,
                      y1=maxy, length=0.1, col='gray40', code=3)
         panel.text(x=ncp/2, y=maxy,
                    expression(frac(delta,'SED')))
         xlim <- current.panel.limits()$xlim
         ylim <- current.panel.limits()$ylim
         panel.arrows(x0=xlim[2], y0=lik0, x1=tcrit,
                      y1=lik0, length=0.1, col=col[1])
         panel.arrows(x0=xlim[2], y0=lik1, x1=tcrit, y1=lik1, 
                      length=0.1, col=col[2])
         panel.arrows(1.35,0.035,2.04,0.025, length=0.1)
         panel.text(1.35,0.04, 'Area = 0.05', pos=2, offset=0.25, col=col[1])
         panel.abline(v=tcrit, col='gray',lty=2)
         panel.arrows(x0=tcrit, y0=-0.04, x1=tcrit, y1=0, 
                      length=0.1, col=col[1])
         panel.text(x=tcrit, y=-0.04, 
                substitute(t[crit]==tc*' for '*italic(p)*phantom(0)*'= 0.05',
                       list(tc=round(tcrit,2))), font=1, pos=1, col=col[1])
      
         panel.axis('right',at=lik0,labels='lik0', text.cex=1.0, line.col=col[1])
         panel.axis('right',at=lik0,line.col=col[1],outside=TRUE)
         panel.axis('right',at=lik1,labels='lik1',text.cex=1.0, line.col=col[2])
         panel.axis('right',at=lik1, line.col=col[2],outside=TRUE)
       },
       main=list(expression('One-sided 2-sample '*t*'-test, power=0.8 with '*alpha*'=0.05'), 
                 just='left',font=1, x=0, y=0)
)
}

## ----cap5, echo=FALSE---------------------------------------------------------
cap5 <- 'This illustrates graphically, for a one-sided $t$-test,
the $t$-statistic for the difference in means required to 
achieve a given power.  For this graph, the $t$-statistic
is calculated with 18 degrees of freedom.  The two density 
curves are separated by the amount that gives \\textit{power} 
= 0.8 for $\\alpha$ = 0.05 .'

## ----pwr-gph, fig.width=7, fig.asp=0.6, fig.cap=cap5, out.width='80%', echo=FALSE, fig.align="center"----
cfdenspwr(nsamp=19)

## ----power, eval=TRUE, echo=FALSE---------------------------------------------
delta2 <- power.t.test(n=19, sd=1.5, sig.level=0.05, power=.8,
                      type='two.sample', 
                      alternative='two.sided')[['delta']]
n <- 19
df2 <- 2*(n-2)
## delta2 is the separation between means
dSTD2 <- delta2/sqrt(2/n)   ## Difference/(SE of difference)
tcrit2 <- qt(.95,df=df2)   
delta1 <- power.t.test(n=19, sd=1.5, sig.level=0.05, power=.8,
                       type='one.sample', 
                       alternative='two.sided')[['delta']]
n1 <- df2+1
dSTD1 <- delta1/sqrt(1/n1)   ## Difference/(SE of difference)
tcrit1 <- qt(.95,df=df2)  

## ----power, eval=FALSE--------------------------------------------------------
#  delta2 <- power.t.test(n=19, sd=1.5, sig.level=0.05, power=.8,
#                        type='two.sample',
#                        alternative='two.sided')[['delta']]
#  n <- 19
#  df2 <- 2*(n-2)
#  ## delta2 is the separation between means
#  dSTD2 <- delta2/sqrt(2/n)   ## Difference/(SE of difference)
#  tcrit2 <- qt(.95,df=df2)
#  delta1 <- power.t.test(n=19, sd=1.5, sig.level=0.05, power=.8,
#                         type='one.sample',
#                         alternative='two.sided')[['delta']]
#  n1 <- df2+1
#  dSTD1 <- delta1/sqrt(1/n1)   ## Difference/(SE of difference)
#  tcrit1 <- qt(.95,df=df2)

## ----power12, echo=FALSE, eval=FALSE------------------------------------------
#  cat("Two-sample, n=9, one-sided:\n")
#  round(c(delta=delta2, tcrit=tcrit2,
#            Power=1 - pt(tcrit2, ncp=dSTD2, df=18)), 3)
#  cat("One-sample, n=19, one-sided\n")
#  round(c(delta=delta1, tcrit=tcrit1,
#    Power=1 - pt(tcrit1, ncp=dSTD1, df=18)), 3)

## ----powercalc----------------------------------------------------------------
power.t.test(type='two.sample', alternative='two.sided', power=0.8,
             sig.level=0.05, sd=1.5, delta=1.4)[['n']]

## ----lrVSpval, eval=TRUE, echo=FALSE------------------------------------------
lrVSpval <- function(pval = 10^(seq(from=-3, to=-0.9, by=.1)), type='two.sample', nsamp=9, twoSided=TRUE, delta=c(1.0,1.4), sd=1.2){
i<-0
lrmat <- matrix(nrow=length(pval), ncol=4)
colnames(lrmat) <- c("pval", "lrDelta1", "lrDelta2", "lrmax")
for(p in pval){
  i <- i+1
  lrmat[i,1:2] <- unlist(tTOlr(pval=p,nsamp=nsamp,delta=delta[1],sd=sd,
                     twoSided=twoSided, showMax=FALSE))[c("pval","lrDelta")]
  lrmat[i,3:4] <- unlist(tTOlr(pval=p,nsamp=nsamp,delta=delta[2],sd=sd,
                     twoSided=twoSided, showMax=TRUE))[c("lrDelta","lrmax")]
}
return(as.data.frame(lrmat))
}

## ----cap6, eval=TRUE, echo=FALSE----------------------------------------------
cap6 <- 'Ratio of likelihood under the alternative to the likelihood 
under the NULL, as a function of the calculated $p$-value, with
$n$ = 9 in each sample in a two-sample test, and 
with $\\delta$ = 0.6$s$ set as the minimum difference of interest.
The graph may, alternatively, be interpreted as for $n$ = 19 in
a one-sample test, now with $\\delta$ = 1.225$s$. The left panel
is for one-sided tests, while the right panel is for two-sided
tests.
'

## ----lrVSpGPH, fig.width=7, fig.asp=0.5, out.width="100%", fig.pos='h', fig.cap=cap6, eval=TRUE, echo=FALSE, warning=FALSE----
suppressPackageStartupMessages(library(latticeExtra))
atx <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)
aty <- c(.5,1,2,5,10,20,50,100, 200)
atxlab <- paste(atx)
atylab <- paste(aty)
delta <- c(1.0,1.4)
df <- lrVSpval(nsamp=9,twoSided=FALSE, delta=delta, sd=1.2,
               type='two.sample')
df2 <- lrVSpval(nsamp=9,twoSided=TRUE, delta=delta, sd=1.2,
                type='two.sample')
df2 <- rbind(df,df2)
df2$type <- rep(c("One-sided","Two-sided"), c(nrow(df),nrow(df)))
lab3 <- vector("expression", 3) 
lab3[1] <- 'Maximum'
for(i in 2:3)lab3[i] <- substitute(expression(delta==del),list(del=delta[i-1]))[2]
gph <- xyplot(lrmax+lrDelta1+lrDelta2 ~ pval | type,data=df2,
              par.settings=simpleTheme(lty=c(1,2,6)),
              auto.key=list(text=lab3, columns=3, points=FALSE, lines=TRUE),   
              scales=list(x=list(log=T, at=atx, labels=atxlab),
                          y=list(log=T, at=aty, labels=atylab), 
                          tck=0.5),
              xlab=expression(italic(p)*'-value'),
              ylab="likelihood ratio",type="l")
update(gph,axis=latticeExtra::axis.grid)

