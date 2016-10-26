library(nimble)
source('benchmark_suite.R')

#'  By default, benchmarkSuite() runs three different models:
#'  'seeds' is the 'seeds' BUGS example model
#'  'birats' is the bivariate normal 'birats' BUGS example model
#'  'linearGaussianSSM' is a univariate linear gaussian state space model with 100 time points
#'  
#'  To run any subset of these models, use the arguments MCMCtests and PFtests.
#'  For example, benchmarkSuite(MCMCtests = "seeds", PFtests = NULL) 
#'  will run only the 'seeds' MCMC test

## run with basic settings 
tst <- benchmarkSuite(niter = c(100000, 10000), nparticles = 100000, nreps = 10, summarize = T, debug = F)
  
## only run MCMC benchmarks, using different numbers of iterations for each benchmark
tst <- benchmarkSuite(PFtests = NULL, niter = c(100000, 10000), nreps = 10, summarize = T, debug = F)

## specify PF and MCMC options
tst <- benchmarkSuite(MCMCs = 'nimble_RW', PFs = 'auxiliary', niter = 100000, nparticles = 100000, nreps = 10, summarize = T, debug = F)
