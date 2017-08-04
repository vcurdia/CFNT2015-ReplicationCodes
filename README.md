# CFNT2015-USTrackRe

[![License](https://img.shields.io/badge/license-BSD%203--clause-green.svg)](https://github.com/vcurdia/CFNT2015-USTrackRe/blob/master/LICENSE)

These codes reproduce the results in:

**Cúrdia, V., Ferrero, A., Ng, G. C., and Tambalotti, A. (2015)**  
[Has U.S. monetary policy tracked the efficient interest rate?](http://www.sciencedirect.com/science/article/pii/S0304393214001457)  
*Journal of Monetary Economics*, 70, pp. 72-83.

[Technical Appendix](https://github.com/vcurdia/CFNT2015-USTrackRe/blob/master/CFNT2015_Appendix.pdf)

These replication codes are available online at:  
https://github.com/vcurdia/CFNT2015-USTrackRe

# Requirements

## Matlab (R)
The codes were tested using Matlab (R) R2014a with the following toolboxes
- Symbolic Toolbox
- Statistical Toolbox
- Optimization Toolbox

## LaTeX
LaTeX is used by some tools to compile certain documents.

`epstopdf`, included in most LaTeX releases, is used by some tools.

## Additional codes and packages

Codes from 
[Vasco Cúrdia](http://www.frbsf.org/economic-research/economists/vasco-curdia/):
- [VC-Tools](https://github.com/vcurdia/VC-Tools), 
  version 
  [v1.8.0](https://github.com/vcurdia/VC-Tools/releases/tag/v1.8.0)
- [VC-BayesianEstimation](https://github.com/vcurdia/VC-BayesianEstimation),
  version
  [v1.6.1](https://github.com/vcurdia/VC-BayesianEstimation/releases/tag/v1.6.1)
  
Codes from [Chris Sims](http://www.princeton.edu/~sims/):
- [gensys](http://sims.princeton.edu/yftp/gensys/)
- [optimize](http://dge.repec.org/codes/sims/optimize/)
- [KF](http://sims.princeton.edu/yftp/Times09/KFmatlab/)

All auxiliary codes included in this repository in subfolders.


# Description of Replication Codes

## Estimation

**Note:** it may take a long time to run with the current settings without
using a parallel computing cluster because with default the following
estimation stages are implemented:

1. Run 20 posterior numerical mode searches starting at different guess
vectors. For each of those it will restart the search up to 30 more times to
confirm that it is not a local mode.

2. Will run three stages of Markov-Chain-Monte-Carlo (MCMC) each for four
   separate chains:  
	2.1. First stage uses the negative of the inverse hessian at the posterior
peak to form the covariance matrix for MCMC draws, scaled to have a rejection
rate between 70 and 80 percent (numerical search). It generates
100,000 draws using a Metropolis algorithm.  
	2.2. After the previous stage, the code discards the first quarter of the
draws for each chain, combines the remain from each chain and computes the
covariance matrix, which is then used as the new covariance for the next stage
of MCMC. It is again numerically rescaled to yield a rejection rate between 70
and 80 percent. It generates 200,000 draws.  
	2.3. We repeat step 2.2. one more time.

In each of the MCMC stages a full report with diagnostics and some inference is
generated. 


## Data

`Baseline\data_greenspan_bernanke_20091204.mat`  
Data for Baseline model

`JPT\JPTData_1954q3_2009q3_3obs.mat`  
Data for JPT model

`Recessions_1987q3_2009q3.mat`  
Data for recession dates: 


## Scripts to estimate specifications reported in paper and appendix

Baseline
- W: `Baseline\NoGapReSetDSGE.m`
- T: `Baseline\BaselineSetDSGE.m`
- W&T: `Baseline\ReSetDSGE.m`
- T with HP Gap: `Baseline\HPSetDSGE.m`
- T with Growth: `Baseline\GrowthSetDSGE.m`
- T with 4Q Growth: `Baseline\Growth4QSetDSGE.m`
- T with HP Gap, estimated lambda: `Baseline\HPEstLambdaSetDSGE.m`
- T with HP Gap, high lambda: `Baseline\HPHighLambdaSetDSGE.m`
- T with Exp Gap: `Baseline\ExpSetDSGE.m`
- T with Exp Gap, estimated lambda: `Baseline\ExpEstLambdaSetDSGE.m`

TVIT
- W: `Baseline\NoGapRePistarSetDSGE.m`
- T: `Baseline\PistarSetDSGE.m`
- W&T: `Baseline\RePistarSetDSGE.m`
- T with HP Gap: `Baseline\HPPistarSetDSGE.m`
- T with Growth: `Baseline\GrowthPistarSetDSGE.m`
- T with 4Q Growth: `Baseline\Growth4QPistarSetDSGE.m`
- T with HP Gap, estimated lambda: `Baseline\HPEstLambdaPistarSetDSGE.m`
- T with HP Gap, high lambda: `Baseline\HPHighLambdaPistarSetDSGE.m`
- T with Exp Gap: `Baseline\ExpPistarSetDSGE.m`
- T with Exp Gap, estimated lambda: `Baseline\ExpEstLambdaPistarSetDSGE.m`

Forward-Looking
- W: `Baseline\EPiNoGapReSetDSGE.m`
- T: `Baseline\EPiExSetDSGE.m`
- W&T: `Baseline\EPiExReSetDSGE.m`
- T with HP Gap: `Baseline\EPiExHPSetDSGE.m`
- T with Growth: `Baseline\EPiEGrowthSetDSGE.m`

Forward-looking with TVIT
- W: `Baseline\Baseline\EPiNoGapRePistarSetDSGE.m`
- T: `Baseline\EPiExPistarSetDSGE.m`
- W&T: `Baseline\EPiExRePistarSetDSGE.m`
- T with HP Gap: `Baseline\EPiExHPPistarSetDSGE.m`
- T with Growth: `Baseline\EPiEGrowthPistarSetDSGE.m`

Four-Quarter Inflation
- W: `Baseline\Pi4qNoGapSetDSGE.m`
- T: `Baseline\Pi4qSetDSGE.m`
- W&T: `Baseline\Pi4qWTSetDSGE.m`
- T with HP Gap: `Baseline\Pi4qHPSetDSGE.m`
- T with Growth: `Baseline\Pi4qGrowthSetDSGE.m`

Four-Quarter Inflation with TVIT
- W: `Baseline\Pi4qNoGapRePistarSetDSGE.m`
- T: `Baseline\Pi4qPistarSetDSGE.m`
- W&T: `Baseline\Pi4qRePistarSetDSGE.m`
- T with HP Gap: `Baseline\Pi4qHPPistarSetDSGE.m`
- T with Growth: `Baseline\Pi4qGrowthPistarSetDSGE.m`

JPT
- W: `JPT\NoGapReSetDSGE.m`
- T: `JPT\BaselineSetDSGE.m`
- W&T: `JPT\ReSetDSGE.m`

JPT with TVIT
- W: `JPT\NoGapRePistarSetDSGE.m`
- T: `JPT\PistarSetDSGE.m`
- W&T: `JPT\RePistarSetDSGE.m`


## Tables and Figures

Tables 1-3 in the paper, D.1 in online appendix  
These tables use elements from each specification's report produced after the
MCMC draws are generated.

Figures B.1 and B.2 in online appendix  
These figures are subsets of the prior-posterior marginal distributions shown
for the individual estimation reports, for the respective specifications.

Figures 1 and 2 in paper  
These two figures are generated by `Make_Figure_1_2.m`

Figure 3 in paper  
This figure is generated by `Make_Figure_3.m`

Figure E1 in appendix  
This figure is generated by `Make_Figure_E1.m`

