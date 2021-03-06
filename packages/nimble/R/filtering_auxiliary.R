##  Contains code to run auxiliary particle filters.
##  We have a build function (buildAuxiliaryFilter),
##  and step function (auxFStep)
##
##  This version of the APF is based on 
##  Pitt et al., 2012

auxStepVirtual <- nimbleFunctionVirtual(
  run = function(m = integer()) {
    returnType(double())
  },
  methods = list(
    returnESS = function() {
      returnType(double())
    }
  )
)

auxFuncVirtual <- nimbleFunctionVirtual(
  methods = list(
    lookahead = function(){}
  )
)

auxLookFunc = nimbleFunction(
    name = 'auxLookFunc',
  contains = auxFuncVirtual,
  setup = function(model, node){
  },
  methods = list(
    lookahead = function(){
    model[[node]] <<- model$getParam(node, 'mean')
    }),
  where = getLoadingNamespace()
)


auxSimFunc = nimbleFunction(
    name = 'auxSimFunc',
  contains = auxFuncVirtual,
  setup = function(model, node){},
  methods = list(
    lookahead = function(){
      simulate(model, node)
    }), where = getLoadingNamespace()
)

auxFStep <- nimbleFunction(
    name = 'auxFStep',
  contains = auxStepVirtual,
  setup = function(model, mvEWSamples, mvWSamples, nodes, iNode, names,
                   saveAll, smoothing, lookahead, silent = TRUE) {
    notFirst <- iNode != 1
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisNode <- nodes[iNode]
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE) 
    #current time point
    t <- iNode
    # Get names of x and xs node for current and previous time point,
    # will be different depending on whether we are saving all time points
    # or only the most recent
    if(saveAll == 1){
      allPrevNodes <- model$expandNodeNames(nodes[1:(iNode-1)])
      prevXName <- prevNode    
      thisXName <- thisNode
      currInd <- t
      prevInd <- t-1
      if(smoothing == TRUE){
        currInd <- 1
        prevInd <- 1
      }
    }
    else{
      allPrevNodes <- names
      prevXName <- names    
      thisXName <- names
      currInd <- 1
      prevInd <- 1 
    }
    isLast <- (iNode == length(nodes))

    auxFuncList <- nimbleFunctionList(auxFuncVirtual) 
    allLatentNodes <- model$expandNodeNames(thisNode)
    numLatentNodes <- length(allLatentNodes)
    if(lookahead == "mean"){
       for(i in 1:numLatentNodes)
         auxFuncList[[i]] <- auxLookFunc(model, allLatentNodes[i])
    }
    else{
      for(i in 1:numLatentNodes)
        auxFuncList[[i]] <- auxSimFunc(model,  allLatentNodes)
    }
    ess <- 0
  },
  run = function(m = integer()) {
    returnType(double())
    ##declare(auxll, double(1,m))
    auxll <- numeric(m, init=FALSE)
    ##declare(auxWts, double(1,m))
    auxWts <- numeric(m, init=FALSE)
    ##declare(wts, double(1,m))
    wts <- numeric(m, init=FALSE)
    ids <- integer(m, 0)
    ##declare(ll, double(1,m))
    ll <- numeric(m, init=FALSE)
    
    ## This is the look-ahead step, not conducted for first time-point
    if(notFirst){ 
      for(i in 1:m) {
        if(smoothing == 1){
          copy(mvEWSamples, mvWSamples, nodes = allPrevNodes, nodesTo = allPrevNodes, row = i, rowTo=i)
        }
        copy(mvWSamples, model, prevXName, prevNode, row=i)        
        calculate(model, prevDeterm)
        if(lookahead == "mean"){
          for(j in 1:numLatentNodes)
              auxFuncList[[j]]$lookahead()
        }
        else
          auxFuncList[[1]]$lookahead()
        
        calculate(model, thisDeterm)
        auxll[i] <- calculate(model, thisData)  # get p(y_t+1 | x_t+1)
        if(is.nan(auxll[i])){
          return(-Inf)
        }
        auxll[i] <- auxll[i]+calculate(model, thisNode) # multiply by p(x_t+1 | x_t)
        auxWts[i] <- auxll[i] + mvWSamples['wts',i][prevInd] # multiply by weight from time t
      }
      normAuxWts <- exp(auxWts)/sum(exp(auxWts))  # normalize weights and resample
      rankSample(normAuxWts, m, ids, silent)
    }   
    for(i in 1:m) {
      if(notFirst) {
        copy(mvWSamples, model, nodes = prevXName, nodesTo = prevNode, row = ids[i])
        calculate(model, prevDeterm) 
      }
      simulate(model, thisNode)  # simulate from x_t+1 | x_t
      copy(model, mvWSamples, nodes = thisNode, nodesTo = thisXName, row=i)
      calculate(model, thisDeterm)
      ll[i]  <- calculate(model, thisData)  # get p(y_t+1 | x_t+1)
      if(is.nan(ll[i])){
        return(-Inf)
      }
      if(notFirst){
        wts[i] <- ll[i]-auxll[ids[i]]  # construct weight following step 4 of paper    
      }
      else{
        wts[i] <- ll[i]  # First step has no auxiliary weights
      }
    }
    
    normWts <- exp(wts)/sum(exp(wts))
    ess <<- 1/sum(normWts^2) 
    
    for(i in 1:m){
      ##  save weights for use in next timepoint's look-ahead step
      mvWSamples['wts', i][currInd] <<- log(normWts[i])   
    }
    rankSample(normWts, m, ids, silent)
    for(i in 1:m){
      if(smoothing == 1){
        copy(mvWSamples, mvEWSamples, nodes = allPrevNodes, nodesTo = allPrevNodes, row = ids[i], rowTo=i)
      }
      copy(mvWSamples, mvEWSamples, thisXName, thisXName, ids[i], i)
    }
    
    ##  Calculate likelihood p(y_t+1 | y_1:t) as in equation (3) of paper
    if(notFirst){
      outLL <- sum(exp(wts))/m
      outLL <- outLL*sum(exp(auxWts))
    }
    else{
      outLL <- sum(exp(wts))/m
    }
    return(log(outLL))
return(0)
  }, 
  methods = list(
    returnESS = function() {
      returnType(double(0))
      return(ess)
    }
  ),
  where = getLoadingNamespace()
)


#' Create an auxiliary particle filter algorithm to estimate log-likelihood.
#' 
#' @description Create an auxiliary particle filter algorithm for a given NIMBLE state space model.  
#'
#' @param model A NIMBLE model object, typically representing a state space model or a hidden Markov model
#' @param nodes A character vector specifying the latent model nodes over which the particle filter will stochastically integrate to estimate the log-likelihood function
#' @param control  A list specifying different control options for the particle filter.  Options are described in the details section below.
#' @author  Nicholas Michaud
#' @details 
#' The \code{control} argument can be specified as follows.
#' \describe{
#' \item{lookahead}{The lookahead function used to calculate auxiliary weights.  Can choose between \code{'mean'} and \code{'simulate'}.
#'  Defaults to \code{'simulate'}.}
#'  \item{saveAll}{Indicates whether to save state samples for all time points (\code{TRUE}), or only for the most recent time point (\code{FALSE})}
#'  \item{smoothing}{Decides whether to save smoothed estimates of latent states, i.e., samples from f(x[1:t]|y[1:t]) if \code{smoothing = TRUE}, or instead to save filtered samples from f(x[t]|y[1:t]) if \code{smoothing = FALSE}. \code{smoothing = TRUE} only works if \code{saveAll = TRUE}.}
#'  \item{timeIndex}{An integer used to manually specify which dimension of the latent state variable indexes time. This need only be set if the number of time points is less than or equal to the size of the latent state at each time point.}
#'  \item{initModel}{A logical value indicating whether to initialize the model before running the filtering algorithm.  Defaults to TRUE.}
#'  }
#' 
#' The auxiliary particle filter modifies the bootstrap filter (\code{\link{buildBootstrapFilter}})
#' by adding a lookahead step to the algorithm: before propagating particles from one time
#' point to the next via the transition equation, the auxiliary filter calculates a weight
#' for each pre-propogated particle by predicting how well the particle will agree with the
#' next data point.  These pre-weights are used to conduct an initial resampling step before
#' propagation. 
#' 
#'  The resulting specialized particle filter algorthm will accept a
#'  single integer argument (\code{m}, default 10,000), which specifies the number
#'  of random \'particles\' to use for estimating the log-likelihood.  The algorithm 
#'  returns the estimated log-likelihood value, and saves
#'  unequally weighted samples from the posterior distribution of the latent
#'  states in the \code{mvWSamples} modelValues object, with corresponding logged weights in \code{mvWSamples['wts',]}.
#'  An equally weighted sample from the posterior can be found in the \code{mvEWsamp} modelValues object.  
#'  
#'   The auxiliary particle filter uses a lookahead function to select promising particles before propagation.  This function can eithre be the expected
#'   value of the latent state at the next time point (\code{lookahead = 'mean'}) or a simulation from the distribution of the latent state at the next time point (\code{lookahead = 'simulate'}), conditioned on the current particle.
#'   
#'  @section \code{returnESS()} Method:
#'  Calling the \code{returnESS()} method of an auxiliary particle filter after that filter has been \code{run()} for a given model will return a vector of ESS (effective
#'  sample size) values, one value for each time point.

#' @export
#' 
#' @family particle filtering methods
#' @references Pitt, M.K., and Shephard, N. (1999). Filtering via simulation: Auxiliary particle filters. \emph{Journal of the American Statistical Association} 94(446): 590-599.
#'   
#' @examples
#' \dontrun{
#' model <- nimbleModel(code = ...)
#' my_AuxF <- buildAuxiliaryFilter(model, 'x[1:100]',
#'    control = list(saveAll = TRUE, lookahead = 'mean'))
#' Cmodel <- compileNimble(model)
#' Cmy_AuxF <- compileNimble(my_AuxF, project = model)
#' logLike <- Cmy_AuxF$run(m = 100000)
#' ESS <- Cmy_AuxF$returnESS(m = 100000)
#' hist(as.matrix(Cmy_Auxf$mvEWSamples, 'x'))
#' }
buildAuxiliaryFilter <- nimbleFunction(
    name = 'buildAuxiliaryFilter',
  setup = function(model, nodes, control = list()) {
    
   
    
    #control list extraction
    saveAll <- control[['saveAll']]
    smoothing <- control[['smoothing']]
    silent <- control[['silent']]
    timeIndex <- control[['timeIndex']]
    lookahead <- control[['lookahead']]
    initModel <- control[['initModel']]
    if(is.null(silent)) silent <- TRUE
    if(is.null(saveAll)) saveAll <- FALSE
    if(is.null(smoothing)) smoothing <- FALSE
    if(is.null(lookahead)) lookahead <- 'simulate'
    if(is.null(initModel)) initModel <- TRUE
    if(!saveAll & smoothing) stop("must have saveAll = TRUE for smoothing to work")
    if(lookahead == "mean"){
      errors <- sapply(model$expandNodeNames(nodes), function(node){tryCatch(getParam(model, node, 'mean'), error=function(a){return("error")})})
      if(any(errors == "error", na.rm=TRUE)) stop("cannot use 'mean' lookahead for this model, try 'simulate'")
    } 
    else if(lookahead != "simulate"){
      stop("lookahead argument must be either 'simulate' or 'mean'")
    }
    
    #latent state info
    varName <- sapply(nodes, function(x){return(model$getVarNames(nodes = x))})
    if(length(unique(varName))>1){
      stop("all latent nodes must come from same variable")
    }
    varName <- varName[1]
    info <- model$getVarInfo(varName)
    latentDims <- info$nDim
    if(is.null(timeIndex)){
      timeIndex <- which.max(info$maxs)
      timeLength <- max(info$maxs)
      if(sum(info$maxs==timeLength)>1) # check if multiple dimensions share the max index size
        stop("unable to determine which dimension indexes time. 
             Specify manually using the 'timeIndex' control list argument")
    } else{
      timeLength <- info$maxs[timeIndex]
    }
    nodes <- paste(info$varName,"[",rep(",", timeIndex-1), 1:timeLength,
                   rep(",", info$nDim - timeIndex),"]", sep="")
    dims <- lapply(nodes, function(n) nimDim(model[[n]]))
    if(length(unique(dims)) > 1) stop('sizes or dimensions of latent states varies')
    vars <- model$getVarNames(nodes =  nodes)  # need var names too
    
    my_initializeModel <- initializeModel(model, silent = silent)
    
    
    # Create mv variables for x state and sampled x states.  If saveAll=TRUE, 
    # the sampled x states will be recorded at each time point. 
    modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[vars]
    if(saveAll){
      
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x)return(x$size))
      mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                              types = type,
                                              sizes = size))
      
      names <- c(names, "wts")
      type <- c(type, "double")
      size$wts <- length(dims)
      if(smoothing == T){
        size$wts <- 1 ##  only need one weight per particle (at time T) if smoothing == TRUE
      }
      mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                              types = type,
                                              sizes = size))
      
    }
    else{
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x)return(x$size))
      size[[1]] <- as.numeric(dims[[1]])
      
      mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                              types = type,
                                              sizes = size))
      
      names <- c(names, "wts")
      type <- c(type, "double")
      size$wts <- 1
      mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                              types = type,
                                              sizes = size))
    }
    
    names <- names[1]
    auxStepFunctions <- nimbleFunctionList(auxStepVirtual)
    for(iNode in seq_along(nodes))
      auxStepFunctions[[iNode]] <- auxFStep(model, mvEWSamples, mvWSamples, nodes,
                                            iNode, names, saveAll, smoothing, lookahead, silent)
    
    essVals <- rep(0, length(nodes))
  },
  run = function(m = integer(default = 10000)) {
    returnType(double())
    declare(logL, double())
    if(initModel == TRUE) my_initializeModel$run()
    resize(mvEWSamples, m) 
    resize(mvWSamples, m)  
    logL <- 0
    for(iNode in seq_along(auxStepFunctions)) {
      logL <- logL + auxStepFunctions[[iNode]]$run(m)
      essVals[iNode] <<- auxStepFunctions[[iNode]]$returnESS()
      
      # when all particles have 0 weight, likelihood becomes NAN
      # this happens if top-level params have bad values - possible
      # during pmcmc for example
      if(is.nan(logL)) return(-Inf)
      if(logL == -Inf) return(logL) 
      if(logL == Inf) return(-Inf) 
    }
  
    return(logL)
  },
  methods = list(
    returnESS = function(){
      returnType(double(1))
      return(essVals)
    }
  ),where = getLoadingNamespace()
)


