## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache=FALSE)

## ---- eval=TRUE---------------------------------------------------------------
suppressPackageStartupMessages(library(hisse))

## ---- eval=TRUE---------------------------------------------------------------
set.seed(42)
phy <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 0.3, mu = 0.2)[[1]]

## ---- eval=TRUE---------------------------------------------------------------
f <- GetFossils(phy, psi=0.05)
fbd.tree <- ProcessSimSample(phy, f)

## ---- eval=TRUE---------------------------------------------------------------
names(fbd.tree)
head(fbd.tree$k.samples)

## ---- eval=FALSE--------------------------------------------------------------
#  k.samples <- data.frame(taxon1="sp12", taxon2="sp12", timefrompresent=3.164384,
#            state=1, stringsAsFactors=FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  k.samples <- data.frame(taxon1="sp12", taxon2="sp12", timefrompresent=3.164384,
#            state1=0, state2=1, stringsAsFactors=FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  turnover <- c(1)
#  eps <- c(1)
#  one.rate <- MiSSE(fbd.tree$phy, f=1, turnover=turnover, eps=eps,
#            includes.fossils=TRUE, k.samples=fdb.tree$k.samples, sann=TRUE,
#            sann.its=1000)

## ---- eval=FALSE--------------------------------------------------------------
#  one.rate <- MiSSE(fbd.tree$phy, f=1, turnover=turnover, eps=eps,
#            includes.fossils=TRUE, k.samples=NULL, sann=TRUE, sann.its=1000)

## ---- eval=FALSE--------------------------------------------------------------
#  margin.test <- MarginReconMiSSE(phy=fbd.tree$phy, f=1, pars=one.rate$solution,
#            hidden.states=1, includes.fossils=TRUE, k.samples=fbd.tree$k.samples,
#            aic=one.rate$AIC)

## ---- eval=FALSE--------------------------------------------------------------
#  trans.rate <- TransMatMakerHiSSE()
#  pp <- hisse(phy=fbd.tree$phy, data, f=c(1,1), turnover=c(1,1), eps=c(1,1),
#              trans.rate=trans.rate,  k.samples=fbd.tree$k.samples,
#              includes.fossils=TRUE)
#  margin.test <- MarginReconHiSSE(phy=fbd.tree$phy, data=data, f=c(1,1),
#              pars=pp$solution, hidden.states=1, includes.fossils=TRUE,
#              k.samples=fbd.tree$k.samples)

## ---- eval=FALSE--------------------------------------------------------------
#  trans.rate <- TransMatMakerMuHiSSE()
#  pp <- MuHiSSE(phy=fbd.tree$phy, data, f=c(1,1,1,1), turnover=c(1,1,1,1),
#              eps=c(1,1,1,1), trans.rate=trans.rate, k.samples=k.samples,
#              includes.fossils=TRUE)
#  margin.test <- MarginReconMuHiSSE(phy=fbd.tree$phy, data=data, f=c(1,1,1,1),
#              pars=pp$solution, hidden.states=1, includes.fossils=TRUE,
#              k.samples=fbd.tree$k.samples)

## ---- eval=TRUE---------------------------------------------------------------
strat.tree <- ProcessSimStrat(phy, f)

## ---- eval=TRUE---------------------------------------------------------------
names(strat.tree)
head(strat.tree$strat.intervals)

## ---- eval=FALSE--------------------------------------------------------------
#  turnover <- c(1)
#  eps <- c(1)
#  one.rate <- MiSSE(strat.tree$phy, f=1, turnover=turnover, eps=eps,
#              includes.fossils=TRUE, k.samples=NULL,
#              strat.intervals=strat.tree$strat.intervals, sann=TRUE,
#              sann.its=5000)

## ---- eval=FALSE--------------------------------------------------------------
#  margin.test <- MarginReconMiSSE(phy=strat.tree$phy, f=1, pars=one.rate$solution,
#            hidden.states=1, includes.fossils=TRUE, k.samples=NULL,
#            strat.intervals=strat.tree$strat.intervals, aic=one.rate$AIC)

## ---- eval=TRUE---------------------------------------------------------------
plot(ladderize(phy), show.tip.label=FALSE, edge.width=0.75)
### Split the table into m and k for the points ###
extinct.samples <- f[which(f$fossiltype_long=="extinct_terminal" | 
            f$fossiltype_long=="extinct_internal"),]
k.samples.tmp <- extinct.samples[which(extinct.samples$has_sampled_descendant == TRUE),]
extinct.samples <- extinct.samples[which(extinct.samples$has_sampled_descendant == FALSE),]
k.samples <- f[which(f$fossiltype_long == "surviving_terminal" | 
            f$fossiltype_long == "surviving_internal"),]
k.samples <- rbind(k.samples, k.samples.tmp)
AddFossilPoints(ladderize(phy), f=extinct.samples, pch=19, cex=0.8, 
            col="#0D79F2")
AddFossilPoints(ladderize(phy), f=k.samples, pch=19, cex=0.8, col="#F25E0D")

## ---- eval=TRUE---------------------------------------------------------------
plot(ladderize(phy), show.tip.label=FALSE, edge.width=0.75)
AddFossilPoints(ladderize(phy), f=extinct.samples, pch=19, cex=0.8, 
            col="#0D79F2")
AddFossilPoints(ladderize(phy), f=k.samples, pch=19, cex=0.8, col="#F25E0D")
AddStratIntervals(ladderize(phy), f=f, pch=19, cex=0.8, col="purple", lwd=2)

## ---- eval=TRUE---------------------------------------------------------------
library(diversitree)
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(4)
phy <- NULL
while( is.null( phy ) ){
  phy <- tree.bisse(pars, max.t=30, x0=0, include.extinct=TRUE)
}
k.samples <- data.frame(taxon1="sp12", taxon2="sp12", timefrompresent=3.164384, 
            state=1, stringsAsFactors=FALSE)
hidden.states=FALSE
phy.k <- hisse:::AddKNodes(phy, k.samples)
fix.type <- hisse:::GetKSampleMRCA(phy.k, k.samples)
nb.tip <- Ntip(phy.k)
nb.node <- phy.k$Nnode
gen <- hisse:::FindGenerations(phy.k)
    
data <- data.frame(taxon=names(phy$tip.state), phy$tip.state, 
           stringsAsFactors=FALSE)
data <- hisse:::AddKData(data, k.samples)
data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
data.new <- data.new[phy.k$tip.label,]
    
dat.tab <- hisse:::OrganizeDataHiSSE(data.new, phy=phy.k, f=c(1,1), 
             hidden.states=FALSE, includes.fossils=TRUE)
edge_details <- hisse:::GetEdgeDetails(phy=phy.k, 
            intervening.intervals=strat.cache$intervening.intervals)
fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
pars.bisse <- c(0.1+0.03, 0.1+0.03, 0.03/0.1, 0.03/0.1, 0.01, 0.01)
 
model.vec <- numeric(48)
model.vec[1:6] = pars.bisse
phy$node.label = NULL
cache <- hisse:::ParametersToPassfHiSSE(model.vec, hidden.states=hidden.states, 
            nb.tip=Ntip(phy.k), nb.node=Nnode(phy.k), bad.likelihood=-300, 
            f=c(1,1), ode.eps=0)
cache$psi <- 0.01
hisse.full <- hisse:::DownPassHiSSE(dat.tab, gen, cache, root.type="madfitz", 
            condition.on.survival=TRUE, root.p=NULL, node=fix.type$node, 
            state=fix.type$state, fossil.taxa=fossil.taxa, 
            fix.type=fix.type$type)

## ---- eval=TRUE---------------------------------------------------------------
dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy.k, f=1, hidden.states=1, includes.fossils=TRUE)
model.vec <- c(0.1+0.03, 0.03/0.1, rep(0,51))
cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, 
            fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, 
            bad.likelihood=exp(-500), ode.eps=0)#
cache$psi <- 0.01
gen <- hisse:::FindGenerations(phy.k)
MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, 
            condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, 
            fossil.taxa=fossil.taxa, node=fix.type$node, fix.type=fix.type$type)

## ---- eval=TRUE---------------------------------------------------------------
library(corHMM)
char.logL <- corHMM(phy.k, data, rate.cat=1, model = "ER", node.states = "none", 
           fixed.nodes=FALSE, p=0.01, root.p="maddfitz")

## ---- eval=TRUE---------------------------------------------------------------
tot.logL <- char.logL$loglik + MiSSE.logL
comparison <- identical(round(hisse.full,3), round(tot.logL,3))
comparison 

## ---- eval=TRUE---------------------------------------------------------------
library(diversitree)
pars <- c(.1,  .15,  .2, .1,
         .03, .045, .06, 0.03,
         .05, .05, .00,
         .05, .00, .05,
         .05, .00, .05,
         .00, .05, .05)
set.seed(2)
phy <- NULL
while( is.null( phy ) ){
    phy <- tree.musse(pars, 30, x0=1, include.extinct=TRUE)
}
k.samples <- data.frame(taxon1="sp20", taxon2="sp37", timefrompresent=8.54554, 
                        state1=0, state2=1, stringsAsFactors=FALSE)
    
phy.k <- hisse:::AddKNodes(phy, k.samples)
fix.type <- hisse:::GetKSampleMRCA(phy.k, k.samples)
nb.tip <- Ntip(phy.k)
nb.node <- phy.k$Nnode
gen <- hisse:::FindGenerations(phy.k)

states <- phy$tip.state
states <- data.frame(phy$tip.state, phy$tip.state,
row.names=names(phy$tip.state))
states <- states[phy$tip.label,]
states.trans <- states
for(i in 1:Ntip(phy)){
    if(states[i,1] == 1){
        states.trans[i,1] = 0
        states.trans[i,2] = 0
    }
    if(states[i,1] == 2){
        states.trans[i,1] = 0
        states.trans[i,2] = 1
    }
    if(states[i,1] == 3){
        states.trans[i,1] = 1
        states.trans[i,2] = 0
    }
    if(states[i,1] == 4){
        states.trans[i,1] = 1
        states.trans[i,2] = 1
    }
}
    
data <- data.frame(taxon=names(phy$tip.state), 
                   states.trans[,1], states.trans[,2], stringsAsFactors=FALSE)
data <- hisse:::AddKData(data, k.samples, muhisse=TRUE)
data.new <- data.frame(data[,2], data[,3], row.names=data[,1])
data.new <- data.new[phy.k$tip.label,]
    
pars.muhisse <- c(rep(0.1+0.03,4), rep(0.03/.1, 4), 0.05,0.05,0, 0.05,0,0.05, 
                  0.05,0,.05, 0,0.05,.05)
model.vec = rep(0,384)
model.vec[1:20] = pars.muhisse
cache <- hisse:::ParametersToPassMuHiSSE(model.vec=model.vec, hidden.states=FALSE, 
            nb.tip=Ntip(phy.k), nb.node=Nnode(phy.k), bad.likelihood=exp(-500), 
            f=c(1,1,1,1), ode.eps=0)
cache$psi <- 0.01
gen <- hisse:::FindGenerations(phy.k)
dat.tab <- hisse:::OrganizeData(data.new, phy.k, f=c(1,1,1,1), 
            hidden.states=FALSE, includes.fossils=TRUE)
fossil.taxa <- which(dat.tab$branch.type == 1)
    
muhisse.full <- hisse:::DownPassMuHisse(dat.tab, gen=gen, cache=cache, 
            root.type="madfitz", condition.on.survival=TRUE, root.p=NULL, 
            node=fix.type$node, state=fix.type$state, fossil.taxa=fossil.taxa, 
            fix.type=fix.type$type)

## Trait independent model should be loglik_tree + loglik_character ##
#Part 1: MiSSE loglik:
dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy.k, f=1, hidden.states=1,   includes.fossils=TRUE)
model.vec <- c(0.1+0.03, 0.03/0.1, rep(0,51))
cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, 
            fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, 
            bad.likelihood=exp(-500), ode.eps=0)#
cache$psi <- 0.01
edge_details <- hisse:::GetEdgeDetails(phy=phy.k, 
            intervening.intervals=strat.cache$intervening.intervals)
fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
gen <- hisse:::FindGenerations(phy.k)
MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, 
            condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, 
            fossil.taxa=fossil.taxa,node=fix.type$node, fix.type=fix.type$type)
    
#Part 2: corHMM loglik:
char.logL <- corHMM(phy.k, data, rate.cat=1, model = "ER", node.states = "none", 
            fixed.nodes=FALSE, p=0.05, root.p="maddfitz")
tot.logL <- char.logL$loglik + MiSSE.logL
    
comparison <- identical(round(muhisse.full,3), round(tot.logL,3))
comparison

