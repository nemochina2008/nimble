
#' Create an MCMC function, from an MCMCspec object
#' 
#' Accepts a single required argument, which must be an object of class MCMCspec.  Returns an MCMC function; see details section.
#' 
#' @param mcmcspec An object of class MCMCspec, which specifys the model, samplers, monitors, and thinning intervals for the resulting MCMC function.  See \code{MCMCspec} for details of this argument.
#' @author Daniel Turek
#' @export
#' @details
#' Calling buildMCMC(mcmcspec) will produce an R mcmc function object, say 'Rmcmc'.
#'
#' The Rmcmc function will have arguments:
#'
#' niter: The number of iterations to run the MCMC.
#'
#' reset: Boolean specifying whether to reset the model and stored samples.  This will simulate into any stochastic nodes with value NA,
#' propagate values through any deterministic nodes, and calculate all model probabilities.
#' This will also reset the internal stored MCMC samples.
#' Specifying reset=FALSE allows the MCMC algorithm to continue running from where it left off.
#' Generally, reset=FALSE should only be used when the MCMC has already been run.  See examples.
#'
#' simulateAll: Boolean specifying whether to simulate into all stochastic nodes.  This will overwrite the current values in all stochastic nodes.
#' 
#' Samples corresponding to the 'monitors' and 'monitors2' from the MCMCspec are stored into the interval variables 'mvSamples' and 'mvSamples2', respectively.
#' These may be accessed using:
#' nfVar(Rmcmc, 'mvSamples')
#' 
#' The Rmcmc function may be compiled to a C MCMC object, taking care to compile in the same project as the R model object, using:
#' Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
#' 
#' The Cmcmc function will function identically as the Rmcmc object, except acting on the C model object.
#' @examples
#' mCode <- modelCode({
#'  mu ~ dnorm(0, 1)
#'  x ~ dnorm(mu, 1)
#' })
#' Rmodel <- nimbleModel(mCode)
#' mcmcspec <- MCMCspec(Rmodel)
#' Rmcmc <- buildMCMC(mcmcspec)
#' Rmcmc(10)
#' samples <- nfVar(Rmcmc, 'mvSamples')
#' samples[['x']]
#' Rmcmc(100, reset = FALSE)
buildMCMC <- nimbleFunction(
    
    setup = function(mcmcspec) {
        model <- mcmcspec$model
        
        RHSonlyNodes <- model$getMaps()$nodeNamesRHSonly
        RHSonlyInitFunctions <- nimbleFunctionList(RHSonlyInit_virtual)
        for(i in seq_along(RHSonlyNodes))     { RHSonlyInitFunctions[[i]] <- RHSonlyInit(model, RHSonlyNodes[i]) }
        
        modelNodes <- model$getNodeNames()
        nodeInitFunctions <- nimbleFunctionList(mcmcNodeInit_virtual)
        for(i in seq_along(modelNodes))     { nodeInitFunctions[[i]] <- mcmcNodeInit(model, modelNodes[i]) }
        
        mvSaved <- modelValues(model)
        samplerFunctions <- nimbleFunctionList(sampler_BASE)
        for(i in seq_along(mcmcspec$samplerSpecs))     { samplerFunctions[[i]] <- mcmcspec$samplerSpecs[[i]]$buildSampler(model=model, mvSaved=mvSaved) }
        
        monitors  <- mcmcspec$monitors
        monitors2 <- mcmcspec$monitors2
        thin  <- mcmcspec$thin
        thin2 <- mcmcspec$thin2
        mvSamples  <- mcmcspec$newMvSamples()
        mvSamples2 <- mcmcspec$newMvSamples2()
    },
    
    run = function(niter = integer(), reset = logical(default=TRUE), simulateAll = logical(default=FALSE)) {
        if(simulateAll)     simulate(model)    ## default behavior excludes data nodes
        if(reset) {
            for(i in seq_along(RHSonlyInitFunctions))  {   RHSonlyInitFunctions[[i]]()   }
            for(i in seq_along(nodeInitFunctions))     {   nodeInitFunctions[[i]]()      }
            nimCopy(from = model, to = mvSaved, row = 1, logProb = TRUE)
            for(i in seq_along(samplerFunctions))      {   nfMethod(samplerFunctions[[i]], 'reset')()   }
            mvSamples_offset  <- 0
            mvSamples2_offset <- 0
            resize(mvSamples,  niter/thin)
            resize(mvSamples2, niter/thin2)
        } else {
            mvSamples_offset  <- getsize(mvSamples)
            mvSamples2_offset <- getsize(mvSamples2)
            resize(mvSamples,  mvSamples_offset  + niter/thin)
            resize(mvSamples2, mvSamples2_offset + niter/thin2)
        }
        
        for(iter in 1:niter) {
            for(i in seq_along(samplerFunctions))      {   samplerFunctions[[i]]()   }
            if(iter %% thin  == 0) { nimCopy(from = model, to = mvSamples,  row = mvSamples_offset  + iter/thin,  nodes = monitors) }
            if(iter %% thin2 == 0) { nimCopy(from = model, to = mvSamples2, row = mvSamples2_offset + iter/thin2, nodes = monitors2) }
        }
    },  where = getLoadingNamespace()
)



##### OLD (v0.1) machinery for handling samplerOrdering (e.g., samplerCalls: c(1,2,1,3,1,4,1))
##### this code below was in the setup():
## create a list of the unique sampler functions
# samplerOrder       <- mcmcspec$samplerOrder
# uniqueSamplerOrder <- unique(samplerOrder)
# numSamplers        <- length(uniqueSamplerOrder)
# samplerFunctions   <- vector('list', numSamplers)
# for(i in seq_along(uniqueSamplerOrder)) {
#     uspec <- mcmcspec$samplerSpecs[[uniqueSamplerOrder[i]]]
#     samplerFunctions[[i]] <- uspec$buildSampler(model=model, mvSaved=mvSaved)
# }
# ## create a list of sampler function calls
# numSamplerCalls <- length(samplerOrder)
# samplerMap <- array(0, c(numSamplerCalls,1))
# for(i in 1:numSamplerCalls) samplerMap[i,1] <- which(uniqueSamplerOrder==samplerOrder[i])
#
#
#### this code in runtime, calling of the samplers:
# for(fc in seq_along(samplerOrder)) {
#     samplerFunctions[[samplerMap[fc,1]]](scale = scales[samplerMap[fc,1]]) }





