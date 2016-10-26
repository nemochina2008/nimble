
benchmarkSuite <- function(
            MCMCtests           = "all",
            MCMCs               = c('nimble', 'nimble'),
            niter               = 10000,
            PFtests             = "all",
            PFs                 = 'bootstrap',
            pfControlList       = list(),
            nparticles          = 100000,
            nreps               = 1,
            setSeed             = TRUE,
            check               = getNimbleOption('checkModel'),
            summarize           = TRUE,
            debug               = FALSE
) {
    ## aliased in MCMCsuiteClass
    suite <- nimBenchmarkClass(MCMCtests, MCMCs, niter, PFtests, PFs,
                               pfControlList, nparticles, nreps,
                               setSeed, check, summarize, debug)
    return(suite$output)
}

nimBenchmarkClass <- setRefClass(

    Class = 'nimBenchmarkClass',
    
    fields = list(
        ## set in initialize()
        MCMCtests = 'character',   ## parsed expression for the model code; must be contained in { ... }    --- ORIGINAL ARGUMENT
        PFtests = 'character',
        RMCMCmodels = 'list',   ## Rmodel object
        RPFmodels = 'list',
        monitors = 'list',    ## the original character vector argument to initialize()    --- ORIGINAL ARGUMENT --- SLIGHTLY MODIFIED
        
        ## set in initialize()
        niter = 'numeric',    ## number of MCMC iterations to run    --- ORIGINAL ARGUMENT
        burnin = 'numeric',   ## burn-in period, the number of initial samples to discard, prior to thinning    --- ORIGINAL ARGUMENT
        thin = 'numeric',   ## thinning interval    --- ORIGINAL ARGUMENT
        nkeep = 'numeric',   ## number of samples we'll keep. equal to (niter/thin - burnin)
        burninFraction = 'numeric',  ## fraction of total sampling effort spent on burnin (burnin / (nkeep + burnin))
        nparticles = 'numeric',
        pfControlList = 'list',
        MCMCs = 'character',   ## character vector of the MCMC analyses.  'winbugs', 'openbugs', 'jags', 'stan', or anything else is nimble    --- ORIGINAL ARGUMENT
        PFs   = 'character',
        latents = 'character',
        
        ## setMCMCdefs()
        MCMCdefs = 'list',   ## named list of {} expression code blocks, corresponding the setup for nimble MCMCs    --- ORIGINAL ARGUMENT --- SLIGHTLY MODIFIED
        MCMCdefNames = 'list',   ## names of the MCMCdefs list
        
        ## setPFdefs()
        PFdefs = 'list',
        PFdefNames = 'list',
        
        setSeed = 'logical',   ## whether to setSeed(0) prior to running each algorithm    --- ORIGINAL ARGUMENT
        debug = 'logical',   ## whether to enter browser() before running each algorithm    --- ORIGINAL ARGUMENT
        nreps = 'numeric',
        runMCMCs = 'logical',
        runPFs   = 'logical',
        output = 'list'   ## list of numeric outputs: samples, summary, timing
    ),
    
    methods = list(
        initialize = function(
          MCMCtestArgs        = 'all',
          MCMCs               = 'nimble',
          niter               = 10000,
          PFtestArgs          = 'all',
          PFs                 = "bootstrap",
          pfControlList       = list(),
          nparticles          = 100000,
          nreps               = 1,
          setSeed             = TRUE,
          check               = getNimbleOption('checkModel'),
          summarize           = TRUE,
          debug               = FALSE) {
          
          runMCMCs <<- if(is.null(MCMCtestArgs)) FALSE else TRUE
          runPFs <<-  if(is.null(PFtestArgs)) FALSE else TRUE
          
          if(debug) browser()
          if(runMCMCs){
            MCMCtests <<- MCMCtestArgs
            if(MCMCtestArgs == 'all')
              MCMCtests <<- c("seeds", "birats")  ## need better test names
            if(length(MCMCs == 1)) MCMCs <<- rep(MCMCs, length(MCMCtests))
            else MCMCs <<- MCMCs
            if(length(niter) == 1) niter <<- rep(niter, length(MCMCtests))
            else niter <<- niter
            RMCMCmodels <<- list()
            monitors <<- list() ## in future, make customizable monitors
            for(i in seq_along(MCMCtests)){
              RMCMCmodels[[i]] <<- nimbleBenchmarkModel(MCMCtests[i])
              monitors[[i]] <<- RMCMCmodels[[i]]$getNodeNames(topOnly = TRUE,  ## in future, make monitors only for MCMC benchmarks
                                                              stochOnly = TRUE)
            }
            setMCMCdefs()
          }
            
          if(runPFs){
            PFtests <<- PFtestArgs
            if(PFtestArgs == 'all')
              PFtests <<- c("linearGaussianSSM")  ## need better test names
            if(length(PFs == 1)) PFs <<- rep(PFs, length(PFs))
            else PFs <<- PFs
            if(length(pfControlList) <= 1) {
              tmpList <- list()
              for(i in 1:length(PFtests))
                tmpList[[i]] <- pfControlList
              
              pfControlList <<- tmpList
            }
            else pfControlList <<- pfControlList
            if(length(nparticles) == 1) nparticles <<- rep(nparticles, length(PFtests))
            else nparticles <<- nparticles
            
            RPFmodels <<- list()
            for(i in seq_along(PFtests)){
              RPFmodels[[i]] <<- nimbleBenchmarkModel(PFtests[i])
              latents[[i]] <<- nimbleLatentName(PFtests[i])
            }
            
            setPFdefs()
          }
          
            nreps <<- nreps
            
            setSeed <<- setSeed
            debug <<- debug

            ## run
            init_output()
            if(debug)              browser()
            run_nimble()
            if(summarize == T) summarizeTimes()
        },
        
        nimbleBenchmarkModel = function(modelName){
          if(modelName == "seeds"){
            dir = nimble:::getBUGSexampleDir('seeds')
            Rmodel <- readBUGSmodel("seeds", dir = dir, data = NULL, inits = NULL, useInits = TRUE,
                                    check = FALSE)
          }
          else if(modelName == "birats"){
            dir = nimble:::getBUGSexampleDir('birats')
            Rmodel <- readBUGSmodel("birats2", dir = dir, data = "birats-data.R", useInits = FALSE,
                                    check = FALSE)
          }
          else if(modelName == "linearGaussianSSM"){
            code <- nimbleCode({
              x[1] ~ dnorm(mean = mu0, sd = sigma_x);
              y[1] ~ dnorm(x[1], sd=sigma_y);
              for(i in 2:N){
                x[i] ~ dnorm(mean = x[i-1], sd = sigma_x);
                y[i] ~ dnorm(mean = x[i], sd=sigma_y);
              }
              sigma_x ~ T(dnorm(1, sd = .5), 0,);
              sigma_y ~ T(dnorm(.1, sd = .5), 0,);
              mu0 <- 0
            })
            
            set.seed(0)
            
            N <- 100
            sigma_x <- 1
            sigma_y <- .1
            x <- rep(NA, N)
            y <- x
            x[1] <- rnorm(1,0,sigma_x)
            y[1] <- rnorm(1,x[1], sigma_y)
            for(i in 2:N){
              x[i] <- rnorm(1,x[i-1], sigma_x)
              y[i] <- rnorm(1,x[i], sigma_y)
            }
            consts <- list(N=N)
            testdata <- list(y=y)
            inits <- list(sigma_x=1, sigma_y=.1, x = x)
            Rmodel <- nimbleModel(code, constants = consts, data = testdata, inits = inits, check = FALSE)
          }
          else stop("not a valid benchmark test name")
          return(Rmodel)
        },
        
        nimbleLatentName = function(modelName){
          if(modelName == 'linearGaussianSSM')
            return('x')
          else
            return(NULL)
        },
        
        setMCMCdefs = function() {
          for(i in seq_along(MCMCtests)){
            MCMCdefs[[i]] <<- list(nimble        = quote(configureMCMC(RMCMCmodels[[MCMCtest]])),
                                   nimble_noConj = quote(configureMCMC(RMCMCmodels[[MCMCtest]], useConjugacy = FALSE)),
                                   nimble_RW     = quote(configureMCMC(RMCMCmodels[[MCMCtest]], onlyRW       = TRUE)),
                                   nimble_slice  = quote(configureMCMC(RMCMCmodels[[MCMCtest]], onlySlice    = TRUE)),
                                   autoBlock     = quote(configureMCMC(RMCMCmodels[[MCMCtest]], autoBlock    = TRUE)))
            MCMCdefNames[[i]] <<- names(MCMCdefs[[i]])
          }
        },
        
        setPFdefs = function(newPFdefs) {
          for(i in seq_along(PFtests)){
            PFdefs[[i]] <<-   list(bootstrap     = quote(buildBootstrapFilter(RPFmodels[[PFtest]], nodes = latents[PFtest], 
                                                                              control = pfControlList[[PFtest]])),
                                   auxiliary     = quote(buildAuxiliaryFilter(RPFmodels[[PFtest]], nodes = latents[PFtest], 
                                                                              control = pfControlList[[PFtest]])),                                   nimble_RW     = quote(configureMCMC(Rmodel, onlyRW       = TRUE)),
                                   enkf          = quote(buildEnsembleKF(RPFmodels[[PFtest]], nodes = latents[PFtest], 
                                                                              control = pfControlList[[PFtest]]))
                                   )
            PFdefNames[[i]] <<- names(PFdefs[[i]])
          }
        },
        
        init_output = function() {
          nMCMCtests <- length(MCMCtests) 
          if(runMCMCs){
            for(i in seq_along(MCMCtests)){
              timing <- rep(NA, nreps+1)
              names(timing) <- c(paste0('MCMC_run_', 1:nreps), 'nimble_compile')
              runParams <- c(niter = niter[i]) 
              initialOutput <- list(timing=timing, runParams = runParams)
              output[[MCMCtests[i]]] <<- initialOutput
            }
          }
          if(runPFs){
            for(i in seq_along(PFtests)){
              timing <- rep(NA, nreps+1)
              names(timing) <- c(paste0('PF_run_', 1:nreps), 'nimble_compile')
              runParams <- c(pfControlList[[i]], pfType = PFs[i])
              initialOutput <- list(timing=timing, runParams = runParams)
              output[[PFtests[i]]] <<- initialOutput
            }
          }
        },
        

        run_nimble = function() {
          if(runMCMCs){
            for(MCMCtest in seq_along(MCMCtests)){
                  mcmcTag <- MCMCs[MCMCtest]
                  mcmcDef <- MCMCdefs[[MCMCtest]][[mcmcTag]]
                  mcmcConf <- eval(mcmcDef)
                  mcmcConf$addMonitors(monitors[[MCMCtest]], print = FALSE)
                  RmcmcFunction <- buildMCMC(mcmcConf)
                  timeResult <- system.time({
                      Cmodel <- compileNimble(RMCMCmodels[[MCMCtest]])
                      CmcmcFunction <- compileNimble(RmcmcFunction, project = RMCMCmodels[[MCMCtest]])
                  })
                  addTimeResult(MCMCtests[MCMCtest], 'nimble_compile', timeResult)
                  browser()
                  if(setSeed) set.seed(0)
                  for(runNum in 1:nreps){
                    timeResult <- system.time({ CmcmcFunction$run(niter[MCMCtest])})
                    addTimeResult(MCMCtests[MCMCtest], paste0('MCMC_run_', runNum), timeResult)
                  }
            }
          }
          if(runPFs){
            for(PFtest in seq_along(PFtests)){
                pfTag <- PFs[PFtest]
                pfDef <- PFdefs[[PFtest]][[pfTag]]
                RpfFunction <- eval(pfDef)
                timeResult <- system.time({
                  Cmodel <- compileNimble(RPFmodels[[PFtest]])
                  CpfFunction <- compileNimble(RpfFunction, project = RPFmodels[[PFtest]])
                })
                addTimeResult(PFtests[PFtest], 'nimble_compile', timeResult)
                if(setSeed) set.seed(0)
                for(runNum in 1:nreps){
                  timeResult <- system.time({ CpfFunction$run(nparticles[PFtest])})
                  addTimeResult(PFtests[PFtest], paste0('PF_run_', runNum), timeResult) 
                }
            }
          }
        },
        
        addTimeResult = function(testNum, timingType, timeResult) {
            output[[testNum]]$timing[timingType] <<- timeResult[[3]]
        },
        
        summarizeTimes = function(){
          for(i in seq_along(output)){
            output[[i]]$timing <<- c(round(summary(output[[i]]$timing[1:nreps]), 2),
                                     St.Dev. = round(sd(output[[i]]$timing[1:nreps]), 2),
                                     round(output[[i]]$timing['nimble_compile'], 2))
          }
        }
    )
)