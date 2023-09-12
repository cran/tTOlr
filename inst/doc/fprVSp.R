## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  fig.align='center',
  width=55,
  collapse = TRUE,
  comment = "#>"
)

## ----set-theme, echo=FALSE----------------------------------------------------
sides <- list(tck = 0.6, pad1 = 0.75, pad2 = 0.75)
theme <- list(axis.line = list(alpha = 1, col = 'gray40',
                               fill = "transparent", lty = 1, lwd = 0.5),
              strip.border = list(alpha = 1, col = rep('gray40', 6),
                                  lty = rep(1, 6), lwd = rep(0.5, 6)),
              strip.shingle = list(alpha = 1, col = rep("gray80", 7)),
              box.3d = list(col = 'gray40'),
              axis.components = list(left = sides, top = sides,
                                     right = sides, bottom = sides),
              fontsize = list(text = 10, points = 6))
library(lattice)
trellis.par.set(theme)
size12 <- list(fontsize=list(text=12, points=8))
size10 <- list(fontsize=list(text=10, points=6))

## ----setup, include=FALSE-----------------------------------------------------
library(tTOlr)
knitr::opts_chunk$set(echo = TRUE, comment=NA)

## ----pkgs, echo=FALSE---------------------------------------------------------
suppressPackageStartupMessages(library(lattice))

## ----cap1, echo=FALSE---------------------------------------------------------
cap1 <- "Under the NULL hypothesis, and assuming
distributional assumptions are correct, $p$-values are uniformly
distributed on the unit interval. The strip plots each show the
distribution of values in a random sample of $p$-values, under
NULL hypothesis assumptions.  The ordering from 1 to 0 reflects
the decrease in the $p$-value, for commonly used test statistics, 
as the absolute value of the test statistic increases. Values 
less than the commonly used 0.05 threshold are shown in red."

## ----H0, fig.width=5.75, fig.height=2.25, out.width='98%', fig.show='hold', fig.cap=cap1, echo=FALSE----
set.seed(17)
rr <- data.frame(x=runif(500), gp=rep(rev(letters[1:5]),rep(100,5)))
panelfun <- function(x,y,...){
                  x1 <- x[x>0.05]
                  y1 <- y[x>0.05]
                  x2 <- x[x<=0.05]
                  y2 <- y[x<=0.05]
                  panel.dotplot(x1,y1,...);
                  panel.dotplot(x2,y2,..., col='red');
                  panel.axis(side = "bottom", at=0.05, outside=TRUE, 
                             ticks=TRUE);
                  panel.abline(v=0.05, col='gray')
}
atx <- c(seq(from=1,to=.2, by=-.2),0.05)
gph <- dotplot(gp~x, data=rr, pch='|', alpha=1, panel=panelfun, 
               xlim=c(1,0), xlab=expression(italic(p)*'-value'), 
               ylab='',scales=list(x=list(at=atx), y=list(draw=F)))
update(gph, main=list("Five random samples of 100 from the uniform", 
                 font=1, y=-1, cex=0.8))

## ----lt5pc, echo=FALSE--------------------------------------------------------
cat(rev(sort(with(rr, round(subset(x, x<=0.05 & gp=='a'), 3)))))
cat(rev(sort(with(rr, round(subset(x, x<=0.05 & gp=='b'), 3)))))

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
print(prettyNum(unlist(shoeStats), digits=3, dropTrailing=TRUE),
      quote=FALSE)

## ----CIs, echo=FALSE----------------------------------------------------------
CI95 <- shoeStats[['Mean']] + c(-1,1)*qt(.975,9)*shoeStats[['SEM']]
CI99 <- shoeStats[['Mean']] + c(-1,1)*qt(.995,9)*shoeStats[['SEM']]

## ----sleep--------------------------------------------------------------------
sleep2 <-with(sleep, Pair(extra[group==2], extra[group==1]))
t <- t.test(sleep2 ~ 1, data = sleep)

## ----sleeplik, echo=FALSE-----------------------------------------------------
maxlrSleep <- with(t, tTOmaxlik(statistic, df=parameter)[['maxlik']]/
                     (2*dt(statistic,df=parameter)))

## ----cutoff, echo=FALSE-------------------------------------------------------
tinfo <- t.test(sleep2 ~ 1, mu=0.8, alternative = 'greater')
t <- tinfo[['statistic']]; df <- tinfo[['parameter']]
maxlrSleep.8 <- 
  with(tinfo, tTOmaxlik(t, df))

## ----maxlrGPH, echo=FALSE, fig.width=6, fig.asp=0.6---------------------------
library(lattice)
maxlrGPH <- function(t, df, tmax, maxlik, offset=0.5, twoSided=TRUE, pretitle=""){
#
# under H1 use non-central t distribution
# ncp is non-centrality paramater
lik0 <- dt(t,df)
if(twoSided){
  maxlrtxt <- paste0(signif(maxlik,2), "/(2*",signif(lik0,3),") = ", signif(maxlik/(2*lik0),2)) 
  pval <- pt(abs(t), df, lower.tail = FALSE)*2
} else {maxlrtxt <- paste0(signif(maxlik,2), "/",signif(lik0,2),") = ", signif(maxlik/(lik0),3))
 pval <- pt(abs(t), df, lower.tail = FALSE)
}
titl <- paste0(pretitle, maxlrtxt)
vals <- pretty(c(-4, 4),50)
vals2 <- pretty(c(-4, 4.16),51)
dframe <- data.frame(x=c(vals,vals2+t), y=c(dt(vals, df=df),
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
          panel.arrows(t-1.2,0.035,t+.25,lik0/2.5, length=0.1,
                       alpha=0.5)
          panel.arrows(-t+1.02,0.035,-t-.25,lik0/2.5, length=0.1,
                       alpha=0.5)
          panel.text(0,0.04, paste0('Area=', round(pval,3),"/2"), alpha=0.75,col=col[1])
          panel.axis("left",at=lik0, outside=TRUE, line.col=col[1],
                    labels='H0', text.col=col[1], text.cex=1.0)
          # panel.arrows(x0=xlim[1], y0=lik0, x1=-t, line.alpha=0.25, 
          #             y1=lik0, length=0.06, col=col[1])
         panel.abline(v=-t, col='gray',lty=1)
         xx <- with(subset(dframe,hyp=='H0'),
                    c(min(x),x[x< -t],-t,-t))
         yy <- with(subset(dframe,hyp=='H0'), c(0,y[x< -t],lik0,0))
         panel.polygon(xx,yy,...,col=col[1],alpha=0.4)
         } else 
           {panel.arrows(t-.5,0.035,t+.25,lik0/2.5, length=0.1,
                         alpha=0.5)
             panel.text(t-.5,0.04, paste('Area =',round(pval,3)), pos=2, offset=0.25, alpha=0.75, col=col[1])
           }
                  xx <- with(subset(dframe,hyp=='H0'),
                    c(t,t,x[x>t],max(x)))
         yy <- with(subset(dframe,hyp=='H0'), c(0,lik0,y[x>t],0))
         panel.polygon(xx,yy,...,col=col[1],alpha=0.4)
         panel.lines(x=c(t,0.54*(t+xlim[2])),y=c(lik0,lik0),
                     lty=2, alpha=0.75,col=col[1])
         panel.lines(x=c(tmax,0.54*(t+xlim[2])),y=c(maxlik,maxlik),
                     lty=2, alpha=0.75,col=col[2])
         panel.axis("right",at=maxlik,labels='H1', outside=TRUE,
                    line.col=col[2],text.col=col[2],text.cex=1.0)
         panel.axis("right",at=maxlik,text.col=col[2], text.cex=1.0,
                    labels=round(maxlik,2),line.col=col[2])
         panel.axis("right",at=lik0, line.col=col[1], labels='H0',
                    outside=TRUE, text.cex=1.0, text.col=col[1])
         panel.axis("right",at=lik0,text.col=col[1],text.cex=1.0,
                    labels=round(lik0,3), line.col=col[1])
       },
       main=list(titl, just="left", font=1, x=0, y=-0.6)
)
}

## ----cap2, echo=FALSE---------------------------------------------------------
cap2 <- paste0("Panel A shows density curves for NULL and for 
the alternative, for a one-sided test with ",quote("$t$")," = ", 
round(t,2)," on ",df," degrees of freedom.  This is the ",
quote("$t$"),"-statistic  for the data on the effect of soporofic
drugs when differences are `B-A-0.8`, i.e., interest is
in the strength of evidence that differences are at least
0.8 hours. A vertical line is placed at the position that 
gives the ",quote("$p$"), "-value, here equal to ", round(tinfo$p.value,3),
".  Panel B shows the normal probability 
plot for the differences.")

## ----DOmaxlr, fig.width=7, fig.asp=0.45, fig.pos='ht', fig.cap=cap2, out.width='100%', echo=FALSE----
gph1 <- maxlrGPH(t, df,
         tmax=maxlrSleep.8[['tmax']], 
         maxlik=maxlrSleep.8[['maxlik']],
         twoSided=FALSE, pretitle=" A: Maximum LR = ")
print(gph1, position=c(0,0,0.55,1))
gph2 <- qqmath(-apply(sleep2, 1, diff)-0.8, ylab="B-A-0.8",
               main=list(' B: Normal probability plot', 
                       just='left', font=1, x=0, y=-0.6))
print(gph2, position=c(0.55,0,1,1),new=FALSE)

## ----cutoff, echo=TRUE, eval=FALSE--------------------------------------------
#  tinfo <- t.test(sleep2 ~ 1, mu=0.8, alternative = 'greater')
#  t <- tinfo[['statistic']]; df <- tinfo[['parameter']]
#  maxlrSleep.8 <-
#    with(tinfo, tTOmaxlik(t, df))

## ----all10, echo=FALSE--------------------------------------------------------
print(prettyNum(unlist(shoeStats), digits=2, dropTrailing=TRUE),
      quote=FALSE)

## ----shoes7, echo=FALSE-------------------------------------------------------
wear7 <- with(MASS::shoes, rbind(A, B, d=B-A))[,1:7]
library(magrittr)
shoeStat7 <- wear7['d',] %>% 
  {list(Mean=mean(.), SD=sd(.), n=length(.))} %$%
  {c(., SEM=SD/sqrt(n), t=Mean*sqrt(n)/SD)} %$%
  {c(., pval=2*pt(t, df=n-1, lower.tail=FALSE), df=n-1)}
colnames(wear7) <- rep("",7)
print(wear7[3,,drop=FALSE])
cat("\n")
print(prettyNum(unlist(shoeStat7), digits=2, dropTrailing=TRUE),
      quote=FALSE)

## ----cap5, echo=FALSE---------------------------------------------------------
tmath <- "$t$"
pmath <- "$p$"
cap5 <- paste0("Density curves for NULL and for the alternative,
for a two-sided test with ",tmath," = ", round(shoeStat7[['t']],2),
" on ", shoeStat7[['df']]," degrees of freedom.
Vertical lines are placed at the positions that give 
the ",pmath,"-value, here equal to ",round(shoeStat7[['pval']],3),
".  Panel B shows the normal probability 
plot for the `B-A` differences in the dataset.")

## ----DOshoeDens, fig.width=5.5, fig.asp=0.6, fig.align='center', fig.cap=cap5, out.width='80%', echo=FALSE----
maxlr7 <- with(shoeStat7, tTOlr::tTOmaxlik(t, df))
maxlrGPH(t=shoeStat7[['t']], df=shoeStat7[['df']],
         tmax=maxlr7[['tmax']], maxlik=maxlr7[['maxlik']],
         twoSided=TRUE, pretitle="Maximum Likelihood Ratio = ")

## ----maxlr7, echo=2:3, eval=TRUE----------------------------------------------
t <- shoeStat7[['t']]
stats7 <- list(t=2.256, df=6)  # t is rounded to 2dp
maxlr7 <- with(stats7, tTOlr::tTOmaxlik(t, df))

## ----sig3, echo=FALSE---------------------------------------------------------
print(setNames(as.character(c(round(unlist(maxlr7)[1:2], 3),
                signif(unlist(maxlr7)[3], 3))), 
               names(maxlr7)), quote=FALSE)

## ----cap3, eval=TRUE, echo=FALSE----------------------------------------------
cap3 <- 'Ratio of the maximum likelihood under the alternative
to the likelihood under the NULL, for three different 
choices of $p$-value, for a range of sample sizes, and for a 
range of degrees of freedom.'

## ----pTOmaxlrGPH, echo=FALSE, fig.width=6, fig.asp=0.6, fig.cap=cap3, fig.pos='ht', out.width='100%'----
lrP <- data.frame(df=rep(c(2,seq(from=5, to=50,by=5)),4), 
                  pval=rep(c(.05,.025,.01,.005),c(11,11,11,11)),
                  tstat=numeric(44), lr=numeric(44))
for(i in 1:nrow(lrP)){
  t <- qt(1-lrP$pval[i]/2,df=lrP$df[i])
  dens0 <- dt(t, df = lrP$df[i])
  lrP[i,c('tstat','lr')] <- c(t, tTOmaxlik(t = t, df = lrP$df[i])[["maxlik"]]/(2*dens0))  ## NB: Two-sided
}
atx <- c(2,(1:5)*10)
aty <- list(c(3,4,5,6,8,10,25,50,100),c(2,2.5,3,4,5))
laby <- list(paste(aty[[1]]),paste(aty[[2]]))
gph <- xyplot(lr+tstat ~ df, groups=pval, data=lrP, type=c("p","l"),
              auto.key=list(columns=4), 
              xlab="Degrees of freedom", ylab="",
                strip=strip.custom(factor.levels=c("Maximum likelihood ratio","t-statistic")),
                scales=list(x=list(at=atx),
                            y=list(relation='free', log=TRUE,
                                   at=aty, labels=laby)), 
                par.settings=lattice::simpleTheme(pch=c(1,2,3,16), alpha.line=0.4,lwd=0.5, lty=2))
update(gph,axis=latticeExtra::axis.grid)

## ----cap4, eval=TRUE, echo=FALSE----------------------------------------------
cap4 <- 'False positive risk, for three different choices
of $p$-value, for a range of sample sizes, and for a 
range of degrees of freedom.'

## ----pTOfprGPH, echo=FALSE, fig.width=6, fig.asp=0.6, fig.cap=cap4, fig.pos='ht', out.width='100%'----
lrP[,"fpr0.1"] <- (1-0.1)/(1-0.1+0.1*lrP[,"lr"])
lrP[,"fpr0.5"] <- (1-0.5)/(1-0.5+0.5*lrP[,"lr"])
aty <- list(c(0.025,0.05,0.1,0.2,0.4),c(0.004,0.01, 0.04,0.1,0.4))
laby <- list(paste(aty[[1]]),paste(aty[[2]]))
gph <- xyplot(fpr0.1+fpr0.5 ~ df, groups=pval, data=lrP, auto.key=list(columns=4), 
                xlab="Degrees of freedom",
                ylab="False positive risk",
                strip=strip.custom(factor.levels=c("Prior = 0.1","Prior = 0.5")),
                scales=list(y=list(relation='free', log=TRUE,
                                   at=aty, labels=laby)), 
                par.settings=lattice::simpleTheme(pch=16))
update(gph,axis=latticeExtra::axis.grid)

## ----explain-power, echo=FALSE, fig.width=5.5, fig.asp=0.8,out.width='80%'----
cfdenspwr <- function(df=18, sig.level=0.05, alternative = 'one.sided', ncp=NA){
vals <- pretty(c(-3.25, 3.25),50)
dframe <- data.frame(x=c(vals,vals+ncp),
                 y=c(dt(vals, df=df, ncp=0),
                     dt(vals+ncp, df=df,ncp=ncp)),
                 hyp=rep(c('H0','H1'),rep(length(vals),2)))
tcrit<-qt((1-sig.level),df,ncp=0)
dens0=dt(tcrit,df,0)
#
# under H1 use non-central t distribution
# ncp is non-centrality paramater
lik1=dt(tcrit,df,ncp=ncp)
atx <- seq(from=-3, to=6, by=1)
# labx <- as.expression( sapply(atx, 
#             function(x) bquote(italic(.(x)*italic(s)))) )
labx <- paste(atx)
col <- trellis.par.get()$superpose.symbol$col 
xyplot(y ~ x, groups=hyp, data=dframe, type='l',
       xlab=expression(phantom(0)),
       sub=expression(italic(t)*'-statistic'),
       ylab='Probability density',
       xlim=range(dframe$x), ylim=c(0,max(dframe$y)*1.125),
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
         panel.text(x=tcrit, y=-0.04, 
                substitute(t[crit]==tc*' for '*italic(p)*phantom(0)*
                             '= 0.05',
                       list(tc=round(tcrit,2))), font=1, pos=1, col=col[1])
         maxy <- dt(0, df)
         ncp <- delta/sed
         panel.arrows(x0=0, y0=maxy, x1=ncp,
                      y1=maxy, length=0.1, col='gray40', code=3)
         panel.text(x=ncp/2, y=maxy,
                    expression(frac(delta,'SED')))
         panel.text(1.35,0.04, 'Area = 0.05', pos=2, offset=0.25,
                    col=col[1])
         xlim <- current.panel.limits()$xlim
         ylim <- current.panel.limits()$ylim
         
         panel.arrows(x0=xlim[2], y0=dens0, x1=tcrit,
                      y1=dens0, length=0.1, col=col[1])
         panel.arrows(x0=xlim[2], y0=lik1, x1=tcrit, y1=lik1, 
                      length=0.1, col=col[2])
         panel.arrows(1.35,0.035,2.04,0.025, length=0.1)
         panel.abline(v=tcrit, col='gray',lty=2)
         panel.arrows(x0=tcrit, y0=-0.04, x1=tcrit, y1=0, 
                      length=0.1, col=col[1])
         panel.axis('right',at=dens0,
                    labels=paste('dens0=',round(dens0,3)), 
                    text.cex=1.0, line.col=col[1])
#         panel.axis('right',at=dens0,line.col=col[1],outside=TRUE)
         panel.axis('right',at=lik1,
                    labels=paste('lik1=',round(lik1,3)),
                    text.cex=1.0, line.col=col[2])
#         panel.axis('right',at=lik1, line.col=col[2],outside=TRUE)
       },
       main=list(expression('One-sided 2-sample '*t*'-test, power=0.8 with '*alpha*'=0.05'), 
                 just='left',font=1, x=0, y=0)
)
}

## ----cap6, echo=FALSE---------------------------------------------------------
cap6 <- 'This illustrates graphically, for a one-sided $t$-test,
the $t$-statistic for the difference in means required to 
achieve a given power.  For this graph, the $t$-statistic
is calculated with 18 degrees of freedom.  The two density 
curves are separated by the amount that gives \\textit{power} 
= 0.8 for $\\alpha$ = 0.05 .'

## ----pwr-gph, fig.width=6, fig.asp=0.6, fig.cap=cap6, out.width='80%', echo=FALSE, fig.align="center"----
n <- 19; df <- 2*(n-1); sd <- 1.5; sed <- sd*sqrt(2/n)
## Calculate difference delta between means that gives power=0.8
delta <- power.t.test(n=19, sd=sd, sig.level=0.05, 
                      power=0.8, type="two.sample", 
                      alternative = "one.sided")[['delta']]
## Calculate the non-centrality parameter ncp
ncp <- delta/sed # sed is Standard Error of Difference
cfdenspwr(df=2*(n-1), ncp=ncp)

## ----power, echo=1:7, eval=FALSE,ref.label='pwr-gph'--------------------------
#  n <- 19; df <- 2*(n-1); sd <- 1.5; sed <- sd*sqrt(2/n)
#  ## Calculate difference delta between means that gives power=0.8
#  delta <- power.t.test(n=19, sd=sd, sig.level=0.05,
#                        power=0.8, type="two.sample",
#                        alternative = "one.sided")[['delta']]
#  ## Calculate the non-centrality parameter ncp
#  ncp <- delta/sed # sed is Standard Error of Difference
#  cfdenspwr(df=2*(n-1), ncp=ncp)

## ----power, eval=FALSE--------------------------------------------------------
#  NA

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
  lrmat[i,1:2] <- tTOlr(pval=p,nsamp=nsamp,delta=delta[1],sd=sd,
                     twoSided=twoSided, showMax=FALSE)[c("pval","lrDelta")]
  lrmat[i,3:4] <- tTOlr(pval=p,nsamp=nsamp,delta=delta[2],sd=sd,
                     twoSided=twoSided, showMax=TRUE)[c("lrDelta","lrmax")]
}
return(as.data.frame(lrmat))
}

## ----cap7, eval=TRUE, echo=FALSE----------------------------------------------
cap7 <- 'Ratio of likelihood under the alternative to the likelihood 
under the NULL, as a function of the calculated $p$-value, with
$n$ = 9 in each sample in a two-sample test, and 
with $\\delta$ = 0.6$s$ set as the minimum difference of interest.
The graph may, alternatively, be interpreted as for $n$ = 19 in
a one-sample test, now with $\\delta$ = 1.225$s$. The left panel
is for one-sided tests, while the right panel is for two-sided
tests.
'

## ----lrVSpGPH, fig.width=6, fig.asp=0.5, out.width="100%", fig.pos='ht', fig.cap=cap7, eval=TRUE, echo=FALSE, warning=FALSE----
# suppressPackageStartupMessages(library(latticeExtra))
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

