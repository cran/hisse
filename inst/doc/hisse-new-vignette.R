## ---- eval=TRUE, echo=FALSE----------------------------------------------
suppressWarnings(library(hisse))
suppressWarnings(library(diversitree))

## ---- eval=FALSE---------------------------------------------------------
#  suppressWarnings(library(diversitree))
#  set.seed(4)
#  # Essentially we are setting up a model that models the evolution of two binary characters
#  # Thus, we are assuming the following state combinations 1=00, 2=10, 3=01, 4=11:
#  pars <- c(0.1,0.1,0.1,0.2, rep(0.03, 4), 0.01,0.01,0,0.01,0,0.01,0.01,0,0.01,0,0.01,0.01)
#  phy <- tree.musse(pars, max.taxa=50, x0=1, include.extinct=FALSE)
#  sim.dat <- data.frame(names(phy$tip.state), phy$tip.state)
#  # Now we want to make the states associated with the second character hidden from us. So,
#  # we remove states 3 and 4 and make them 1 and 2
#  sim.dat[sim.dat[,2]==3,2] = 1
#  sim.dat[sim.dat[,2]==4,2] = 2
#  # This next step simply forces the character to be binary:
#  sim.dat[,2] = sim.dat[,2] - 1

## ---- eval=FALSE---------------------------------------------------------
#  turnover <- c(1,1)
#  extinction.fraction <- c(1,1)
#  f <- c(1,1,1,1)

## ---- eval=TRUE----------------------------------------------------------
trans.rates.bisse <-  TransMatMakerHiSSE(hidden.traits=0)
print(trans.rates.bisse)

## ---- eval=FALSE---------------------------------------------------------
#  dull.null <- hisse(phy=phy, data=sim.dat, f=f, turnover=turnover,
#                       eps=extinction.fraction, hidden.states=FALSE,
#                       trans.rate=trans.rates.bisse)

## ---- eval=FALSE---------------------------------------------------------
#  turnover <- c(1,2)
#  extinction.fraction <- c(1,1)
#  BiSSE <- hiSSE(phy=phy, data=sim.dat, f=f, turnover=turnover,
#                       eps=extinction.fraction, hidden.states=FALSE,
#                       trans.rate=trans.rates.bisse)

## ---- eval=FALSE---------------------------------------------------------
#  turnover <- c(1,2,3,4)
#  extinction.fraction <- rep(1, 4)
#  f = c(1,1)

## ---- eval=TRUE----------------------------------------------------------
trans.rate.hisse <- TransMatMakerHiSSE(hidden.traits=1)
print(trans.rate.hisse)

## ---- eval=FALSE---------------------------------------------------------
#  HiSSE <- hisse(phy=phy, data=states.trans, f=f, turnover=turnover,
#                       eps=extinction.fraction, hidden.states=TRUE,
#                       trans.rate=trans.rate.hisse)

## ---- eval=FALSE---------------------------------------------------------
#  turnover <- c(1, 1, 2, 2)
#  extinction.fraction <- rep(1, 4)
#  f = c(1,1)
#  trans.rate <- TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
#  extinction.fraction <- rep(1, 8)
#  trans.rate <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)

