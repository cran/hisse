## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=TRUE, echo=FALSE----------------------------------------------
suppressPackageStartupMessages(library(deSolve))
suppressPackageStartupMessages(library(hisse))
load("cladeSamplingVariables.Rsave")

## ---- eval=TRUE, tidy=TRUE, tidy.opts=list(width.cutoff=60)--------------
GetProbs <- function(yini, times){
        times=times
        prob.subtree.cal.full <- lsoda(yini, times, func = "maddison_DE_bisse", padded.pars, initfunc="initmod_bisse", dllname = "hisse", rtol=1e-8, atol=1e-8)
        probs.out <- prob.subtree.cal.full[-1,-1]
    return(probs.out)
}

## ---- eval=TRUE, tidy=TRUE, tidy.opts=list(width.cutoff=60)--------------
times=c(0, 20)
yini <-c(E0=1-0.5, E1=1-0.5, D0=0.5, D1=0)
left.branch.probs <- GetProbs(yini, times)
left.branch.probs[1:2]

## ---- eval=TRUE, tidy=TRUE, tidy.opts=list(width.cutoff=60)--------------
times=c(0, 10)
yini <-c(E0=1-0.5, E1=1-0.5, D0=0.5, D1=0)
right.subA.probs <- GetProbs(yini, times)
right.subB.probs <- GetProbs(yini, times)

nodeR <- c(right.subA.probs[3:4] *  right.subB.probs[3:4] * c(0.1, 0.2))
#Arbritarily using the extinction probabilities from side A:
phi_A <- right.subA.probs[1:2]

## ---- eval=TRUE, tidy=TRUE, tidy.opts=list(width.cutoff=60)--------------
times=c(10, 20)
yini <-c(E0=phi_A[1], E1=phi_A[2], D0=nodeR[1], D1=nodeR[2])
right.branch.probs <- GetProbs(yini, times)
right.branch.probs[1:2]

## ---- eval=TRUE, tidy=TRUE, tidy.opts=list(width.cutoff=60)--------------
round(left.branch.probs[1:2],4) == round(right.branch.probs[1:2],4)

## ---- eval=TRUE, tidy=TRUE, tidy.opts=list(width.cutoff=60)--------------
times=c(0, 10)
yini <-c(E0=1-0.25, E1=1-0.25, D0=0.25, D1=0)
right.subA.probs <- GetProbs(yini, times)
right.subB.probs <- GetProbs(yini, times)
nodeR <- c(right.subA.probs[3:4] *  right.subB.probs[3:4] * c(0.1, 0.2))
#Again, arbritarily using the extinction probabilities from side A:
phi_A <- right.subA.probs[1:2]

times=c(10, 20)
yini <-c(E0=phi_A[1], E1=phi_A[2], D0=nodeR[1], D1=nodeR[2])
right.branch.probs <- GetProbs(yini, times)
right.branch.probs[1:2]
left.branch.probs[1:2]
round(left.branch.probs[1:2],4) == round(right.branch.probs[1:2],4)

