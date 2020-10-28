## ---- eval=FALSE---------------------------------------------------------
#  library( devtools )
#  install_github(repo = "thej022214/hisse", ref = "master")

## ---- eval=TRUE----------------------------------------------------------
suppressWarnings(library(hisse))
suppressWarnings(library(diversitree))

## ---- eval=TRUE----------------------------------------------------------
## Generate a list with the parameters of the model:
pars <- SimulateGeoHiSSE(hidden.traits = 1, return.GeoHiSSE_pars = TRUE)
pars

## ---- eval=TRUE----------------------------------------------------------
pars$model.pars[,1] <- c(0.1, 0.1, 0.1, 0.03, 0.03, 0.05, 0.05)
pars$model.pars[,2] <- c(0.2, 0.2, 0.2, 0.03, 0.03, 0.05, 0.05)
pars$q.01[1,2] <- pars$q.01[2,1] <- 0.005
pars$q.0[1,2] <- pars$q.0[2,1] <- 0.005
pars$q.1[1,2] <- pars$q.1[2,1] <- 0.005
pars

## ---- eval=TRUE----------------------------------------------------------
set.seed(42)
sim.geohisse <- SimulateGeoHiSSE(pars=pars, hidden.traits = 1, x0 = "01A", max.taxa = 500)
phy <- sim.geohisse$phy
phy$node.labels <- NULL
sim.dat <- data.frame(taxon=sim.geohisse$data[,1], ranges=as.numeric(sim.geohisse$data[,2]))

## ---- eval=FALSE---------------------------------------------------------
#  ## Model 1 - Dispersal parameters vary only, no range-dependent diversification.
#  turnover <- c(1,1,0)
#  eps <- c(1,1)
#  trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
#  trans.rate.mod <- ParEqual(trans.rate, c(1,2))
#  mod1 <- GeoHiSSE(phy = phy, data = sim.dat, f=c(1,1,1),
#                    turnover=turnover, eps=eps,
#                    hidden.states=FALSE, trans.rate=trans.rate.mod,
#                    turnover.upper=100, trans.upper=10)

## ---- eval=FALSE---------------------------------------------------------
#  ## Model 2. Canonical GeoSSE model, range effect on diversification
#  turnover <- c(1,2,3)
#  eps <- c(1,1)
#  trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
#  trans.rate.mod <- ParEqual(trans.rate, c(1,2))
#  mod2 <- GeoHiSSE(phy = phy, data = sim.dat, f=c(1,1,1),
#                    turnover=turnover, eps=eps,
#                    hidden.states=FALSE, trans.rate=trans.rate.mod,
#                    turnover.upper=100, trans.upper=10)

## ---- eval=FALSE---------------------------------------------------------
#  ## Model 3. GeoHiSSE model with 1 hidden trait, no range-dependent diversification.
#  ## Note below how parameters vary among hidden classes but are the same within each
#  ##      hidden class.
#  turnover <- c(1,1,0,2,2,0)
#  eps <- c(1,1,1,1)
#  trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, make.null=TRUE)
#  trans.rate.mod <- ParEqual(trans.rate, c(1,2))
#  mod3 <- GeoHiSSE(phy = phy, data = sim.dat, f=c(1,1,1),
#                    turnover=turnover, eps=eps,
#                    hidden.states=TRUE, trans.rate=trans.rate.mod,
#                    turnover.upper=100, trans.upper=10)

## ---- eval=FALSE---------------------------------------------------------
#  ## Model 4. GeoHiSSE model with 1 hidden trait, no range-dependent diversification.
#  turnover <- c(1,2,3,4,5,6)
#  eps <- c(1,1,1,1)
#  trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1)
#  trans.rate.mod <- ParEqual(trans.rate, c(1,2))
#  mod4 <- GeoHiSSE(phy = phy, data = sim.dat, f=c(1,1,1),
#                    turnover=turnover, eps=eps,
#                    hidden.states=TRUE, trans.rate=trans.rate.mod,
#                    turnover.upper=100, trans.upper=10)

## ---- eval=FALSE---------------------------------------------------------
#  ## Model 5. MuSSE-like model with no hidden trait, no cladogenetic effects.
#  turnover <- c(1,2,0)
#  eps <- c(1,1)
#  trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, make.null=FALSE,
#                                      separate.extirpation = TRUE)
#  trans.rate.mod <- ParEqual(trans.rate, c(1,2))
#  trans.rate.mod <- ParEqual(trans.rate.mod, c(2,3))
#  mod5 <- GeoHiSSE(phy = phy, data = sim.dat, f=c(1,1,1),
#                    turnover=turnover, eps=eps,
#                    hidden.states=FALSE, trans.rate=trans.rate.mod,
#                    turnover.upper=100, trans.upper=10,
#                    assume.cladogenetic = FALSE)

## ---- eval=TRUE----------------------------------------------------------
load( "geohisse_new_vignette.Rsave" )
GetAICWeights(list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4), criterion="AIC")
## As the number of models in the set grows, naming each model in the set can become hard.
## So one can use a list (created by some automated code) as an input also:
list.geohisse <- list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4)
GetAICWeights(list.geohisse, criterion="AIC")

## ---- eval=FALSE---------------------------------------------------------
#  recon.mod1 <- MarginReconGeoSSE(phy = mod1$phy, data = mod1$data, f = mod1$f,
#                                   pars = mod1$solution, hidden.states = 1,
#                                   root.type = mod1$root.type, root.p = mod1$root.p,
#                                   aic = mod1$AIC, n.cores = 4)
#  recon.mod2 <- MarginReconGeoSSE(phy = mod2$phy, data = mod2$data, f = mod2$f,
#                                   pars = mod2$solution, hidden.states = 1,
#                                   root.type = mod2$root.type, root.p = mod2$root.p,
#                                   aic = mod2$AIC, n.cores = 4)
#  recon.mod3 <- MarginReconGeoSSE(phy = mod3$phy, data = mod3$data, f = mod3$f,
#                                   pars = mod3$solution, hidden.states = 2,
#                                   root.type = mod3$root.type, root.p = mod3$root.p,
#                                   aic = mod3$AIC, n.cores = 4)
#  recon.mod4 <- MarginReconGeoSSE(phy = mod4$phy, data = mod4$data, f = mod4$f,
#                                   pars = mod4$solution, hidden.states = 2,
#                                   root.type = mod4$root.type, root.p = mod4$root.p,
#                                   aic = mod4$AIC, n.cores = 4)

## ---- eval=TRUE----------------------------------------------------------
## Load previous results:
load( "geohisse_recons_new_vignette.Rsave" )

## ---- eval=TRUE----------------------------------------------------------
recon.models <- list(recon.mod1, recon.mod2, recon.mod3, recon.mod4)
model.ave.rates <- GetModelAveRates(x = recon.models, type = "tips")

## ---- eval=TRUE----------------------------------------------------------
head( model.ave.rates )

## ----fig1, fig.height = 15, fig.width = 5--------------------------------
plot.geohisse.states(x = recon.models, rate.param = "net.div", type = "fan",
                     show.tip.label = FALSE, legend = FALSE)

