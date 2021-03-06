\name{GeoHiSSE}
\alias{GeoHiSSE}
\title{Hidden Geographic State Speciation and Extinction}
\description{Sets up and executes a GeoHiSSE model (Hidden Geographic State Speciation
 and Extinction) on a very large phylogeny and character distribution.}
\usage{
GeoHiSSE(phy, data, f=c(1,1,1), turnover=c(1,2,3), eps=c(1,2), 
hidden.states=FALSE, trans.rate=NULL, assume.cladogenetic=TRUE, 
condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=TRUE,
sann.its=1000, bounded.search=TRUE,  max.tol=.Machine$double.eps^.50,
mag.san.start=0.5, starting.vals=NULL, turnover.upper=1000, 
eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0, dt.threads=1)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and
    with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a matrix (or dataframe) with two columns containing species
  information. First column has the species names and second column has
  area codes. Values for the areas need to be 0, 1, or 2, where 0 is the 
  widespread area '01', 1 is endemic area '00' and 2 is endemic area '11'. See 'Details'.}
\item{f}{vector of length 3 with the estimated proportion of extant
  species in state 1 (area '00'), state 2 (area '1'), and state 0 (widespread area '01') 
  that are included in the phylogeny. A value of c(0.25, 0.25, 0.5) means that 25 percent
  of species in areas '00' and '11' and 50 percent of species in area '01' are
  included in the phylogeny. By default all species are assumed to be
  sampled.}
\item{turnover}{a numeric vector of length equal to 3+(number of
  \code{hidden.states} * 3). A GeoSSE model has 3 turnover parameters:
  tau00, tau11 and tau01. A GeoHiSSE model with one hidden area has 6 speciation
parameters: tau00A, tau11A, tau01A, tau00B, tau11B, and tau01B, and so on. The length of the numeric vector needs to match the number of speciation parameters in the model.}
\item{eps}{a numeric vector of length equal to 2+(number of
  \code{hidden.states} * 2). A GeoSSE model has 2 extinct fraction parameters:
  ef00 and ef11. A GeoHiSSE model with one hidden area has 4 extinct.frac
parameters: ef00A, ef11A, ef00B, and ef11B, and so on. The length of the numeric
vector needs to match the number of extinct.frac parameters in the
model.}
\item{hidden.states}{a logical indicating whether the model includes hidden.states.
 The default is \code{FALSE}.} 
\item{trans.rate}{provides the transition rate model. See function
  \code{TransMatMakerGeoHiSSE}.}
\item{assume.cladogenetic}{assumes that cladogenetic events occur at nodes. The 
	default is \code{TRUE}.}
\item{condition.on.survival}{a logical indicating whether the likelihood
  should be conditioned on the survival of two lineages and the
  speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root summarization follow the procedure described by FitzJohn et al. 2009, \dQuote{madfitz} or Herrera-Alsina et al. 2018, \dQuote{herr_als}.}
\item{root.p}{a vector indicating fixed root state probabilities. The
  default is \code{NULL}. Order of the areas in the vector need to follow:
  root.p[1] = 1 (endemic area '0'); root.p[2] = 2 (endemic
  area '1'); root.p[3] = 0 (widespread area '01').}
\item{sann}{a logical indicating whether a two-step optimization
  procedure is to be used. The first includes a simulate annealing
  approach, with the second involving a refinement using
  \code{subplex}. The default is \code{TRUE}.}
\item{sann.its}{a numeric indicating the number of times the simulated
  annealing algorithm should call the objective function.}
\item{bounded.search}{a logical indicating whether or not bounds should
  be enforced during optimization. The default is is \code{TRUE}.}
\item{max.tol}{supplies the relative optimization tolerance to
  \code{subplex}.}
\item{mag.san.start}{Sets the extinction fraction to estimate the starting values 
  for the diversification parameters. The equation used is based on Magallon and 
  Sanderson (2001), and follows the procedure used in the original GeoSSE 
  implementation.}
\item{starting.vals}{a numeric vector of length 3 with starting values for the
  model for all areas and hidden states. Position [1] sets turnover, [2]
  sets extinction fraction, and [3] dispersal rates.}
\item{turnover.upper}{sets the upper bound for the speciation parameters.}
\item{eps.upper}{sets the upper bound for the extirpation parameters.}
\item{trans.upper}{sets the upper bound for the transition rate parameters.}
\item{restart.obj}{an object that contains everything to restart an optimization.}
\item{ode.eps}{sets the tolerance for the integration at the end of a
  branch. Essentially if the sum of compD is less than this tolerance,
  then it assumes the results are unstable and discards them. The
  default is set to zero, but in testing a value of 1e-8 can sometimes
  produce stable solutions for both easy and very difficult optimization problems.}
\item{dt.threads}{sets the number of threads available to data.table. In practice this need not change from the default of 1 thread, as we have not seen any speedup from allowing more threads.}
}
\details{
This function sets up and executes a more complex and faster version of the GeoHiSSE model (for the original function see \code{GeoHisse.old}). One of the main differences here is that the model allows up to 10 hidden categories, and implements a more efficient means of carrying out the branch calculation. Specifically, we break up the tree into carry out all descendent branch calculations simultaneously, combine the probabilities based on their shared ancestry, then repeat for the next set of descendent . In testing, we've found that as the number of taxa increases, the calculation becomes much more efficient. In future versions, we will likely allow for multicore processing of these calculations to further improve speed. Also, note this function has replaced the version of \code{GeoSSE} that is currently available (see \code{GeoHisse.old}).

The other main difference is that, like \code{HiSSE}, we employ a modified optimization procedure. In other words, rather than optimizing birth and death separately, \code{GeoHisse} optimizes orthogonal transformations of these variables: we let tau = birth+death define "net turnover", and we let eps = death/birth define the \dQuote{extinction fraction}. This reparameterization alleviates problems associated with overfitting when birth and death are highly correlated, but both matter in explaining the diversity pattern. 

To setup a model, users input vectors containing values to indicate how
many free parameters are to be estimated for each of the variables in
the model. This is done using the \code{turnover} and
\code{extinct.frac} parameters. One needs to specify a value for each of
the parameters of the model, when two parameters show the same value,
then the parameters are set to be linked during the estimation of the
model. For example, a GeoHiSSE model with 1 hidden area and all free
parameters has \code{turnover = 1:6}. The same model with
speciation rates constrained to be the same for all hidden states has
\code{turnover = c(1,2,3,1,2,3)}. This same format applies to
\code{extinct.frac}.

Once the model is specified, the parameters can be estimated using the
subplex routine (default), or use a two-step process (i.e., sann=TRUE)
that first employs a stochastic simulated annealing procedure, which is
later refined using the subplex routine.

The \dQuote{trans.rate} input is the transition model and has an
entirely different setup than speciation and extirpation rates.
See \code{TransMatMakerGeoHiSSE} function for more details. 

For user-specified \dQuote{root.p}, you should specify the probability
for each area. If you are doing a hidden model, there will be six areas:
0A, 1A, 2A, 0B, 1B, 2B. So if you wanted to say the root had to be in
area 0 (widespread distribution), you would specify \dQuote{root.p =
  c(0.5, 0, 0, 0.5, 0, 0)}. In other words, the root has a 50\% chance
to be in one of the areas 0A or 0B.

For the \dQuote{root.type} option, we are currently maintaining the previous default of 
\dQuote{madfitz}. However, it was recently pointed out by Herrera-Alsina et al. (2018)
that at the root, the individual likelihoods for each possible state should be conditioned
prior to averaging the individual likelihoods across states. This can be set doing 
\dQuote{herr_als}. It is unclear to us which is exactly correct, but it does seem that 
both \dQuote{madfitz} and \dQuote{herr_als} behave exactly as they should in the case of 
character-independent diversification (i.e., reduces to likelihood of tree + likelihood of
trait model). We've also tested the behavior and the likelihood differences are very 
subtle and the parameter estimates in simulation are nearly indistinguishable from the
\dQuote{madfitz} conditioning scheme. We provide both options and encourage users to try 
both and let us know conditions in which the result vary dramatically under the two root
implementations. We suspect they do not. 

Also, note, that in the case of \dQuote{root.type=user} and \dQuote{root.type=equal} are
no longer explicit \dQuote{root.type} options. Instead, either \dQuote{madfitz} or 
\dQuote{herr_als} are specified and the \dQuote{root.p} can be set to allow for custom
root options.  
}
\value{
\code{GeoHiSSE} returns an object of class \code{geohisse.fit}. This is a list with 
elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the model 
parameters.}
\item{$index.par}{an index matrix of the parameters being estimated.}
\item{$f}{user-supplied sampling frequencies.}
\item{$hidden.states}{a logical indicating whether hidden states were included in the 
model.}
\item{$assume.cladogenetic}{a logical indicating whether cladogenetic events were allowed
at nodes.}
\item{$condition.on.surivival}{a logical indicating whether the
  likelihood was conditioned on the survival of two lineages and the speciation event 
subtending them.}
\item{$root.type}{indicates the user-specified root prior assumption.}
\item{$root.p}{indicates whether the user-specified fixed root probabilities.}
\item{$phy}{user-supplied tree}
\item{$data}{user-supplied dataset}
\item{$trans.matrix}{the user-supplied transition matrix}
\item{$max.tol}{relative optimization tolerance.}
\item{$starting.vals}{The starting values for the optimization.}
\item{$upper.bounds}{the vector of upper limits to the optimization search.}
\item{$lower.bounds}{the vector of lower limits to the optimization
  search.}
\item{$ode.eps}{The ode.eps value used for the estimation.}
}
\references{
Caetano, D.S., B.C. O'Meara, and J.M. Beaulieu. 2018. Hidden state models improve state-dependent diversification approaches, including biogeographic models. Evolution, 72:2308-2324.

Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., W.P. Maddison, and S.P. Otto. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Maddison W.P., P.E. Midford, and S.P. Otto. 2007. Estimating a binary characters effect on speciation and extinction. Syst. Biol. 56:701-710.

Nee S., R.M. May, and P.H. Harvey. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}

\author{Jeremy M. Beaulieu}
\keyword{models}
