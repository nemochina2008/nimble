#' Executes multiple MCMC algorithms and organizes results.
#'
#' Creates, runs, and organizes output from a suite of MCMC algorithms, all applied to the same model, data, and initial values.
#' This can include WinBUGS, OpenBUGS, JAGS and Stan MCMCs, as well as NIMBLE MCMC algorithms.
#' Trace plots and density plots for the MCMC samples may also be generated and saved.
#'
#' @details
#' Creates and runs an MCMC Suite.
#' By default, this will execute the specified MCMCs, record all samples, generate summary statistics, and create and save trace plots and posterior density plots.
#' This default behavior can ben altered via a variety of arguments.
#' Following execution of the MCMC algorithms, returns a named list containing \code{samples}, \code{summary}, and \code{timing} elements.
#' See the NIMBLE User Manual for more information about the organization of the return object.
#' 
#' @param code The quoted code expression representing the model, such as the return value from a call to \code{nimbleCode}).
#' No default value, this is a required argument.
#'
#' @param constants A named list giving values of constants for the model.
#' This is the same as the \code{constants} argument which would be passed to \code{nimbleModel}.
#' Default value is list().
#'
#' @param data A named list giving the data values for the model.
#' This is the same as the \code{data} argument which would be passed to \code{nimbleModel} or \code{model$setData}.
#' Default value is \code{list()}.
#' 
#' @param inits A named list giving the initial values for the model.
#' This is the same as the \code{inits} argument which would be passed to \code{nimbleModel} or \code{model$setInits}.
#' Default value is \code{list()}.
#'
#' @param monitors A character vector giving the node names or variable names to monitor.
#' The samples corresponding to these nodes will be stored in the output samples, will have summary statistics calculated, and density and trace plots generated.
#' Default value is all top-level stochastic nodes of the model.
#' 
#' @param niter Number of MCMC iterations to run.
#' This applies to all MCMC algorithms in the suite.
#' Default value is 10,000.
#'
#' @param burnin Number of initial, post-thinning, MCMC iterations to discard.
#' Default value is 2,000.
#' 
#' @param thin Thinning interval for the MCMC samples.
#' This applies to all MCMC algorithms in the suite.  The thinning occurs prior to the burnin samples being discarded.
#' Default value is 1.
#' 
#' @param summaryStats A character vector, specifying the summary statistics to calculate on the MCMC samples.
#' Each element may be the character name of an exisiting R function (possibly user-defined) which acts on a numeric vector and returns a scalar (e.g., \code{mean} or \code{sd},
#' or a character string which when parsed and evaluted will define such a function (e.g., \code{function(x) mean(sqrt(x))}).
#' Default value is \code{c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp')}, where the final two elements are functions which calculate the limits of a 95 percent Bayesian credible interval.
#' 
#' @param calculateEfficiency A logical, specifying whether to calculate the efficiency for each MCMC algorithm.  Efficiency is defined as the effective sample size (ESS) of each model parameter divided by the algorithm runtime (in seconds).  Default is FALSE.
#'
#' @param MCMCs A character vector specifying the MCMC algorithms to run.
#' \code{'winbugs'} specifies WinBUGS;
#' \code{'openbugs'} specifies OpenBUGS;
#' \code{'jags'} specifies JAGS;
#' \code{'stan'} specifies Stan; in this case, must also provide the \code{'stan_model'} argument;
#' \code{'nimble'} specifies NIMBLE's default MCMC algorithm;
#' \code{'nimble_noConj'} specifies NIMBLE's default MCMC algorithm without the use of any conjugate Gibbs sampling;
#' \code{'nimble_RW'} specifies NIMBLE MCMC algorithm using only random walk Metropolis-Hastings (\code{'RW'}) samplers;
#' \code{'nimble_slice'} specifies NIMBLE MCMC algorithm using only slice (\code{'slice'}) samplers;
#' \code{'autoBlock'} specifies NIMBLE MCMC algorithm with block sampling of dynamically determined parameter groups attempting to maximize sampling efficiency;
#' Anything else will be interpreted as NIMBLE MCMC algorithms, and must have associated entries in the MCMCdefs argument.
#' Default value is \code{'nimble'}, which specifies NIMBLE's default MCMC algorithm.
#' 
#' @param MCMCdefs A named list of MCMC definitions.  The names of list elements should corespond to any custom MCMC algorithms specified in the \code{MCMCs} argument.
#' The list elements should be quoted expressions, enclosed in {} braces.  When executed, the internal code must return an MCMC configuration object, 
#' specifying the corresponding MCMC algorithm; in particular, setting the appropriate samplers.  The code may assume existance of the R model object \code{Rmodel},
#' and must *return* the MCMC configuration object.  Therefore, the final line of such a code block would frequently be a standalone \code{MCMCconf}, to return this object.
#' 
#' @param winbugs_directory A character string giving the directory of the executable WinBUGS program for the WinBUGS MCMC.
#' This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
#' Default value is \code{'C:/WinBUGS14'}.
#' 
#' @param winbugs_program A character string giving the name of the WinBUGS program, for the WinBUGS MCMC.
#' This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
#' Default value is \code{'WinBUGS'}.
#'
#' @param openbugs_directory A character string giving the directory of the executable OpenBUGS program for the OpenBUGS MCMC.
#' This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
#' Default value is \code{'C:/OpenBUGS323'}.
#' 
#' @param openbugs_program A character string giving the name of the OpenBUGS program, for the OpenBUGS MCMC.
#' This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
#' Default value is \code{'OpenBUGS'}.
#' 
#' @param stan_model A character string specifying the location and name of the model file (\code{'modelName.stan'}) for use with the Stan MCMC program.
#' This argument must include the \code{'.stan'} extension, and must be provided whenever the \code{MCMCs} argument includes \code{'stan'}.
#' 
#' @param stan_inits A character string specifying the location and name of the inits file (\code{'modelName.init.R'}) for use with the Stan MCMC program.
#' This argument must include the \code{'.init.R'} extension, and must be provided whenever the \code{MCMCs} argument includes \code{'stan'}.
#' If omitted, it will attempt to locate an inits file in the same directory as the Stan model file.
#'
#' @param stan_data A character string specifying the location and name of the data file (in the form \code{'modelName.data.R'}) for use with the Stan MCMC program.
#' This argument must include the \code{'.data.R'} extension, and must be provided whenever the \code{MCMCs} argument includes \code{'stan'}.
#' If omitted, it will attempt to locate a data file in the same directory as the Stan model file.
#'
#' @param stanNameMaps A list specifying name mappings between Stan and WinBUGS/OpenBUGS.
#' The syntax for list elements is list(BUGS_PARAM_NAME = list(StanSourceName = 'STAN_PARAM_NAME', transform = function(x) TRANSFORMATION_FUNCTION(x))).
#' The transformation is optional.
#' 
#' @param makePlot Logical argument, specifying whether to generate the trace plots and posterior density plots, for each monitored node.
#' Default value is \code{TRUE}.
#' 
#' @param savePlot Logical argument, specifying whether to save the trace plots and density plots.
#' Plots will be saved into the current working directory.
#' Only used when \code{makePlot == TRUE}.
#' Default value is \code{TRUE}.
#' 
#' @param plotName Character string, giving the file name for saving the trace plots and density plots.
#' Only used when \code{makePlot == TRUE} and \code{savePlot == TRUE}.
#' Default value is \code{'MCMCsuite'}.
#'
#' @param setSeed Logical argument, specifying whether to set.seed(0) prior to MCMC sampling.
#' Default value is \code{TRUE}.
#' 
#' @param check Logical argument, specifying whether to check the model object for missing or invalid values.  Default is given by the NIMBLE option 'checkModel', see help on \code{nimbleOptions} for details.
#' 
#' @param debug Logical argument, specifying whether to enter a \code{browser()} at the onset of executing each MCMC algrithm.
#' For use in debugging individual MCMC algorithms, if necessary.
#' Default value is FALSE.
#'
#' @param ... For internal use only
#'
#' @return Returns a named list containing elements:
#' samples: A 3-dimensional array containing samples from each MCMC algorithm.
#' summary: A 3-dimensional array containing summary statistics for each variable and algorithm.
#' timing: A numeric vector containing timing information.
#' efficiency: Minimum and mean sampling efficiencies for each algorithm (only provided if option calculateEfficiency = TRUE).
#' See the NIMBLE User Manual for more information about the organization of the return object.
#'  
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#' })
#' output <- MCMCsuite(code,
#'                     data = list(x=3),
#'                     inits = list(mu=0),
#'                     niter = 10000,
#'                     monitors = 'mu',
#'                     MCMCs = c('nimble', 'nimble_RW'),
#'                     summaryStats = c('mean', 'sd', 'max', 'function(x) max(abs(x))'),
#'                     makePlot = FALSE)
#' }
#' 
#' @author Daniel Turek
#' @export
MCMC_CODESS <- function(
  code,
  constants           = list(),
  data                = list(),
  inits               = list(),
  monitors            = character(),
  niter               = 10000,
  burnin              = 2000,
  thin                = 1,
  tuning		  = list(),
  summaryStats        = c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp'),
  calculateEfficiency = FALSE,
  MCMCs               = 'nimble',
  MCMCdefs            = list(),
  targetNames        =list(),
  winbugs_directory   = 'C:/WinBUGS14',
  winbugs_program     = 'WinBUGS',
  openbugs_directory  = 'C:/OpenBUGS323',
  openbugs_program    = 'OpenBUGS',
  stan_model          = '',
  stan_inits          = NULL,
  stan_data           = NULL,
  stanNameMaps        = list(),
  makePlot            = TRUE,
  savePlot            = TRUE,
  plotName            = 'MCMC_CODESS',
  setSeed             = TRUE,
  check               = getNimbleOption('checkModel'),
  debug               = FALSE) {
  ## aliased in MCMC_CODESSClass
  CODESS <- MCMC_CODESSClass(code, constants, data, inits, monitors, niter, burnin, thin, tuning, summaryStats, calculateEfficiency,
                             MCMCs, MCMCdefs, targetNames, winbugs_directory, winbugs_program, openbugs_directory, openbugs_program,
                             stan_model, stan_inits, stan_data, stanNameMaps, makePlot, savePlot, plotName, setSeed,
                             check, debug)
  return(CODESS$output)
}

codess<-function(x, tuning){
  Nchain = round(length(x)/2)
  N <- tuning$gridsize
  gridsize = c(N,N)
  bandwidth=c(tuning$bw, tuning$bw)
  range.x = tuning$range.x
  
  est <- bkde2D(x, bandwidth=bandwidth, gridsize = gridsize, range.x = range.x)
  K <- est$fhat
  K <- t(K) 
  for(j in 1:N) K[,j] <- K[,j] / sum(K[,j])
  
  ### Generating a chain for the discrete kernel
  x <- est$x1
  
  discreteChainIndices <- integer(Nchain)
  discreteChainIndices[1] <- round(length(x)/2)
  for(i in 2:Nchain) {
    discreteChainIndices[i] <- sample(1:N, 1,
                                      prob = K[,discreteChainIndices[i-1]])
  }
  discreteChain <- x[discreteChainIndices]
  #effectiveSize(discreteChain)
  
  
}

sampler_record_wrapper <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control){
      numSamples <- 0
      before <- c(0, 0)
      after <- c(0, 0)
      samplerFunctionList <- nimbleFunctionList(sampler_BASE)
    ###### make sure to provide *named* arguments to this function
    ###### shouldn't require anything in control$control, if you don't want
    controlListForNestedSampler <- mcmc_generateControlListArgument(samplerFunction = control$sampler_function, control = control$control)
    samplerFunctionList[[1]] <- eval(call( control$sampler_function, model = model, mvSaved = mvSaved, target = target, control =  controlListForNestedSampler))}, 
    run = function() {
      ## these lines are new:
      numSamples <<- numSamples + 1
      setSize(before, numSamples)
      setSize(after, numSamples)
      before[numSamples] <<- model[[target]]
      ## back to the original sampler function code:
      
      samplerFunctionList[[1]]$run()
      
      ## this line new:
      after[numSamples] <<- model[[target]]
    },
    methods = list(
      reset = function() {samplerFunctionList[[1]]$reset()}
    ))


sampler_conjugate_wrapper <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control){
      numSamples <- 0
      before <- c(0, 0)
      after <- c(0, 0)
      samplerFunctionList <- nimbleFunctionList(sampler_BASE)
    ###### make sure to provide *named* arguments to this function
    ###### shouldn't require anything in control$control, if you don't want
    
      conjugacyResult <- model$checkConjugacy2(target)[[target]]
      if(is.null(conjugacyResult))     stop('non-conjugate target \'', target, '\' in conjugate sampler')
      prior <- conjugacyResult$prior
      dependentCounts <- sapply(conjugacyResult$control, length)
      names(dependentCounts) <- gsub('^dep_', '', names(dependentCounts))
      conjSamplerName <- createDynamicConjugateSamplerName(prior = prior, dependentCounts = dependentCounts)
      if(!dynamicConjugateSamplerExists(conjSamplerName)) {
        conjSamplerDef <- conjugacyRelationshipsObject$generateDynamicConjugateSamplerDefinition(prior = prior, dependentCounts = dependentCounts)
        dynamicConjugateSamplerAdd(conjSamplerName, conjSamplerDef)
      }
      conjSamplerFunction <- dynamicConjugateSamplerGet(conjSamplerName)
      
      
      controlListForNestedSampler <- mcmc_generateControlListArgument(samplerFunction = conjSamplerFunction, control = conjugacyResult$control)
    samplerFunctionList[[1]] <- eval(call( "conjSamplerFunction", model = model, mvSaved = mvSaved, target = target, control =  controlListForNestedSampler))}, 
    run = function() {
      ## these lines are new:
      numSamples <<- numSamples + 1
      setSize(before, numSamples)
      setSize(after, numSamples)
      before[numSamples] <<- model[[target]]
      ## back to the original sampler function code:
      
      samplerFunctionList[[1]]$run()
      
      ## this line new:
      after[numSamples] <<- model[[target]]
    },
    methods = list(
      reset = function() {samplerFunctionList[[1]]$reset()}
    ))


sampler_conjugate <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control){
      samplerFunctionList <- nimbleFunctionList(sampler_BASE)
    ###### make sure to provide *named* arguments to this function
    ###### shouldn't require anything in control$control, if you don't want
    
      conjugacyResult <- model$checkConjugacy2(target)[[target]]
      if(is.null(conjugacyResult))     stop('non-conjugate target \'', target, '\' in conjugate sampler')
      prior <- conjugacyResult$prior
      dependentCounts <- sapply(conjugacyResult$control, length)
      names(dependentCounts) <- gsub('^dep_', '', names(dependentCounts))
      conjSamplerName <- createDynamicConjugateSamplerName(prior = prior, dependentCounts = dependentCounts)
      if(!dynamicConjugateSamplerExists(conjSamplerName)) {
        conjSamplerDef <- conjugacyRelationshipsObject$generateDynamicConjugateSamplerDefinition(prior = prior, dependentCounts = dependentCounts)
        dynamicConjugateSamplerAdd(conjSamplerName, conjSamplerDef)
      }
      conjSamplerFunction <- dynamicConjugateSamplerGet(conjSamplerName)
      
      
      controlListForNestedSampler <- mcmc_generateControlListArgument(samplerFunction = conjSamplerFunction, control = conjugacyResult$control)
    samplerFunctionList[[1]] <- eval(call( "conjSamplerFunction", model = model, mvSaved = mvSaved, target = target, control =  controlListForNestedSampler))}, 
    run = function() {
      
      samplerFunctionList[[1]]$run()
      
      
    },
    methods = list(
      reset = function() {samplerFunctionList[[1]]$reset()}
    ))




## Construct MCMCconf from SamplerList
BuildMCMCconfOld <- function(CandidateSamplerList, targetNames, monitor){
  n = length(CandidateSamplerList)
  p = length(monitor)
  MCMCdefs<-vector(mode="list", length=(n+1))
  names(MCMCdefs)[1]<-'Combine'
  
  str1 <- "{
  mcmcConf <- configureMCMC(Rmodel)\n"
  
  str1<- paste0(str1,"mcmcConf$removeSamplers('",targetNames[[1]],"')\n")
  for (i in 1:n){
    str1<-paste0(str1,"mcmcConf$addSampler(target = '",targetNames[[1]],"', type = sampler_record_wrapper, control = list(sampler_function = '",CandidateSamplerList[[i]]$type,"', control = list(")
    m = length(CandidateSamplerList[[i]]$control)
    if (m==1){
      str1 <-paste0(str1, names(CandidateSamplerList[[i]]$control[1]),"=",CandidateSamplerList[[i]]$control[1])
    } else if (m>1){
      str1 <-paste0(str1, names(CandidateSamplerList[[i]]$control[1]),"=",CandidateSamplerList[[i]]$control[1])
      for(j in 1:(m-1)){
        str1 <-paste0(str1,",", names(CandidateSamplerList[[i]]$control[j+1]),"=",CandidateSamplerList[[i]]$control[j+1])
        
      }
      
    }
    str1 <-paste0(str1,"), name = '",CandidateSamplerList[[i]]$name,"'))\n")
  }
  str2<-"mcmcConf$addMonitors(c('"
    if (p==1){
      str2 <-paste0(str2, monitor[1])
    } else if (p>1){
      str2 <-paste0(str2, monitor[1]) 
      for(k in 2:p){
        str2 <-paste0(str2,"','",monitor[k]) 
      }
     }
    
    str2<- paste0(str2,"'))\n mcmcConf \n }")
    
  str1<-paste0(str1,str2)
  
  MCMCdefs[[1]] =parse(text=str1) 
  
  for(i in 1:n){
    names(MCMCdefs)[i+1]<-CandidateSamplerList[[i]]$name
    str1 <- "{
    mcmcConf <- configureMCMC(Rmodel)\n"
    
    str1<- paste0(str1,"mcmcConf$removeSamplers('",targetNames[[1]],"')\n")
    str1<-paste0(str1,"mcmcConf$addSampler(target = '",targetNames[[1]],"', type = sampler_record_wrapper, control = list(sampler_function = '",CandidateSamplerList[[i]]$type,"', control = list(")
    m = length(CandidateSamplerList[[i]]$control)
    if (m==1){
      str1 <-paste0(str1, names(CandidateSamplerList[[i]]$control[1]),"=",CandidateSamplerList[[i]]$control[1])
    } else if (m>1){
      str1 <-paste0(str1, names(CandidateSamplerList[[i]]$control[1]),"=",CandidateSamplerList[[i]]$control[1])
      for(j in 1:(m-1)){
        str1 <-paste0(str1,",", names(CandidateSamplerList[[i]]$control[j+1]),"=",CandidateSamplerList[[i]]$control[j+1])
        
      }
      
    }
    str1 <-paste0(str1,"), name = '",CandidateSamplerList[[i]]$name,"'))\n")
    
    str1<-paste0(str1,str2)
    
    MCMCdefs[[i+1]] =parse(text=str1) 
  }  
  return (MCMCdefs)
}

BuildDefaultConf1 <- function(DefaultSamplerList, monitor){
  n = length(monitor)
  str2<-"\nconf$addMonitors(c('"
    if (n==1){
      str2 <-paste0(str2, monitor[1])
    } else if (n>1){
      str2 <-paste0(str2, monitor[1]) 
      for(k in 2:n){
        str2 <-paste0(str2,"','",monitor[k]) 
      }
     }
    
    str2<- paste0(str2,"'))\n")
   
    str1 <- "conf <- configureMCMC(Rmodel)\n"
  
  for (i in 1:n){
    str1<- paste0(str1,"conf$removeSamplers('",monitor[i],"')\n")
  }
  for (i in 1:n){
    if(regexpr('conjugate', DefaultSamplerList[[i]]$type)>0){
      node <- monitor[i]
      conjInfo <- Rmodel$checkConjugacy2(node)[[node]]
      str1<- paste0(str1,"\nconf$addConjugateSampler(",conjInfo,") \n")
      #str1 <-
      #        paste0(
      #          str1, "\n conf$addSampler(target = '",monitor[i],"', type = sampler_conjugate_wrapper, control=list(), name = '",DefaultSamplerList[[i]]$name,"') \n"
 #             )
    }
    else if(DefaultSamplerList[[i]]$type=='sampler_RW_block'){
      if(length(DefaultSamplerList[[i]]$target)>1){
        str1<-paste0(str1,"conf$addSampler(target = c('",DefaultSamplerList[[i]]$target[1])
        for(j in 2: length(DefaultSamplerList[[i]]$target))
          str1<-paste0(str1,"','",DefaultSamplerList[[i]]$target[j])
      } else {
        str1<-paste0(str1,"conf$addSampler(target = c('",monitor[i],"','",monitor[-i][1])
      }
      str1<-paste0(str1,"'), type =", DefaultSamplerList[[i]]$type,", control = list(")
    } else{
      str1<-paste0(str1,"\nconf$addSampler(target = '",monitor[i],"', type ='", DefaultSamplerList[[i]]$type,"', control = list(")
    }
    
    m = length(DefaultSamplerList[[i]]$control)
    if (m==1){
      str1 <-paste0(str1, names(DefaultSamplerList[[i]]$control[1]),"=",DefaultSamplerList[[i]]$control[1])
    } else if (m>1){
      str1 <-paste0(str1, names(DefaultSamplerList[[i]]$control[1]),"=",DefaultSamplerList[[i]]$control[1])
      for(j in 1:(m-1)){
        str1 <-paste0(str1,",", names(DefaultSamplerList[[i]]$control[j+1]),"=",DefaultSamplerList[[i]]$control[j+1])
        
      }
      
    }
    if(!regexpr('conjugate', DefaultSamplerList[[i]]$type)>0)
      str1 <-paste0(str1,"), name = '",names(DefaultSamplerList)[i],"')\n")
  
    
  }
  str1<-paste0(str1,str2)
  eval(str1)
  return (conf)
}


## Construct MCMCconf from SamplerList
BuildMCMCconf <- function(CandidateSamplerList, targetNames, nontargetNames, monitor){
  n = length(CandidateSamplerList$target)
  n1= length(CandidateSamplerList$nontarget)
  
  p = length(monitor)
  str2<-"mcmcConf$addMonitors(c('"
    if (p==1){
      str2 <-paste0(str2, monitor[1])
    } else if (p>1){
      str2 <-paste0(str2, monitor[1]) 
      for(k in 2:p){
        str2 <-paste0(str2,"','",monitor[k]) 
      }
     }
    
    str2<- paste0(str2,"'))\n mcmcConf \n }")
   
  MCMCdefs<-vector(mode="list", length=(n+n1+1))
  
  for(l in 0:n1){
  
    if(l==0){
      names(MCMCdefs)[1]<-'ConjCombine'
        str1 <- "{
  mcmcConf <- configureMCMC(Rmodel)\n"
  
  str1<- paste0(str1,"mcmcConf$removeSamplers('",targetNames[[1]],"')\n")
  for (i in 1:n){
    if(CandidateSamplerList$target[[i]]$type=='sampler_RW_block'){
    str1<-paste0(str1,"mcmcConf$addSampler(target = c('",targetNames[[1]],"','", nontargetNames[[1]],"'), type = sampler_record_wrapper, control = list(sampler_function = '",CandidateSamplerList$target[[i]]$type,"', control = list(")
    } else{
      str1<-paste0(str1,"mcmcConf$addSampler(target = '",targetNames[[1]],"', type = sampler_record_wrapper, control = list(sampler_function = '",CandidateSamplerList$target[[i]]$type,"', control = list(")
    }
    
    m = length(CandidateSamplerList$target[[i]]$control)
    if (m==1){
      str1 <-paste0(str1, names(CandidateSamplerList$target[[i]]$control[1]),"=",CandidateSamplerList$target[[i]]$control[1])
    } else if (m>1){
      str1 <-paste0(str1, names(CandidateSamplerList$target[[i]]$control[1]),"=",CandidateSamplerList$target[[i]]$control[1])
      for(j in 1:(m-1)){
        str1 <-paste0(str1,",", names(CandidateSamplerList$target[[i]]$control[j+1]),"=",CandidateSamplerList$target[[i]]$control[j+1])
        
      }
      
    }
    str1 <-paste0(str1,"), name = '",CandidateSamplerList$target[[i]]$name,"'))\n")
  
    
  }
  str1<-paste0(str1,str2)
  
  MCMCdefs[[1]] =parse(text=str1) 
    }
    else{
      names(MCMCdefs)[l+1]<-paste0(CandidateSamplerList$nontarget[[l]]$name,'Combine')
  str1 <- "{
  mcmcConf <- configureMCMC(Rmodel)\n"
  str1<- paste0(str1,"mcmcConf$removeSamplers('",nontargetNames[[1]],"')\n")

  str1<-paste0(str1,"mcmcConf$addSampler(target = '",nontargetNames[[1]],"', type = sampler_record_wrapper, control = list(sampler_function = '",CandidateSamplerList$nontarget[[l]]$type,"', control = list(")
    m = length(CandidateSamplerList$nontarget[[l]]$control)
    if (m==1){
      str1 <-paste0(str1, names(CandidateSamplerList$nontarget[[l]]$control[1]),"=",CandidateSamplerList$nontarget[[l]]$control[1])
    } else if (m>1){
      str1 <-paste0(str1, names(CandidateSamplerList$nontarget[[l]]$control[1]),"=",CandidateSamplerList$nontarget[[l]]$control[1])
      for(j in 1:(m-1)){
        str1 <-paste0(str1,",", names(CandidateSamplerList$nontarget[[l]]$control[j+1]),"=",CandidateSamplerList$nontarget[[l]]$control[j+1])
        
      }
      
    }
    str1 <-paste0(str1,"), name = '",CandidateSamplerList$target[[l]]$name,"'))\n")
  
    
    
  str1<- paste0(str1,"mcmcConf$removeSamplers('",targetNames[[1]],"')\n")
  for (i in 1:n){
    
    if(CandidateSamplerList$target[[i]]$type=='sampler_RW_block'){
    str1<-paste0(str1,"mcmcConf$addSampler(target = c('",targetNames[[1]],"','", nontargetNames[[1]],"'), type = sampler_record_wrapper, control = list(sampler_function = '",CandidateSamplerList$target[[i]]$type,"', control = list(")
    } else{
      str1<-paste0(str1,"mcmcConf$addSampler(target = '",targetNames[[1]],"', type = sampler_record_wrapper, control = list(sampler_function = '",CandidateSamplerList$target[[i]]$type,"', control = list(")
    }
    m = length(CandidateSamplerList$target[[i]]$control)
    if (m==1){
      str1 <-paste0(str1, names(CandidateSamplerList$target[[i]]$control[1]),"=",CandidateSamplerList$target[[i]]$control[1])
    } else if (m>1){
      str1 <-paste0(str1, names(CandidateSamplerList$target[[i]]$control[1]),"=",CandidateSamplerList$target[[i]]$control[1])
      for(j in 1:(m-1)){
        str1 <-paste0(str1,",", names(CandidateSamplerList$target[[i]]$control[j+1]),"=",CandidateSamplerList$target[[i]]$control[j+1])
        
      }
      
    }
    str1 <-paste0(str1,"), name = '",CandidateSamplerList$target[[i]]$name,"'))\n")
  
  }
  str1<-paste0(str1,str2)
  
  MCMCdefs[[l+1]] =parse(text=str1)
      
      
    
  
     
    }
  }
  for(i in 1:n){
    names(MCMCdefs)[n1+i+1]<-CandidateSamplerList$target[[i]]$name
    str1 <- "{
    mcmcConf <- configureMCMC(Rmodel)\n"
    
    str1<- paste0(str1,"mcmcConf$removeSamplers('",targetNames[[1]],"')\n")
    
    if(CandidateSamplerList$target[[i]]$type=='sampler_RW_block'){
      str1<- paste0(str1,"mcmcConf$removeSamplers('",nontargetNames[[1]],"')\n")
      str1<-paste0(str1,"mcmcConf$addSampler(target = c('",targetNames[[1]],"','",nontargetNames[[1]],"'), type = sampler_record_wrapper, control = list(sampler_function = '",CandidateSamplerList$target[[i]]$type,"', control = list(")
    }
    else{
      str1<-paste0(str1,"mcmcConf$addSampler(target = '",targetNames[[1]],"', type = sampler_record_wrapper, control = list(sampler_function = '",CandidateSamplerList$target[[i]]$type,"', control = list(")
    }
    
    m = length(CandidateSamplerList$target[[i]]$control)
    if (m==1){
      str1 <-paste0(str1, names(CandidateSamplerList$target[[i]]$control[1]),"=",CandidateSamplerList$target[[i]]$control[1])
    } else if (m>1){
      str1 <-paste0(str1, names(CandidateSamplerList$target[[i]]$control[1]),"=",CandidateSamplerList$target[[i]]$control[1])
      for(j in 1:(m-1)){
        str1 <-paste0(str1,",", names(CandidateSamplerList$target[[i]]$control[j+1]),"=",CandidateSamplerList$target[[i]]$control[j+1])
        
      }
      
    }
    str1 <-paste0(str1,"), name = '",CandidateSamplerList$target[[i]]$name,"'))\n")
    
    str1<-paste0(str1,str2)
    
    MCMCdefs[[n1+i+1]] =parse(text=str1) 
  }  
  return (MCMCdefs)
}

## Construct MCMCconf from CandidateSamplerList
BuildCombinedConf <- function(DefaultSamplerList, CandidateSamplerList, targetNames, GroupLM, monitor) {
  MCMCdefs <- vector(mode = "list", length = 1)
  n = length(monitor)
  
  p = length(monitor)
  str2 <- "\n mcmcConf$addMonitors(c('"
  if (p == 1) {
    str2 <- paste0(str2, monitor[1])
  } else if (p > 1) {
    str2 <- paste0(str2, monitor[1])
    for (k in 2:p) {
      str2 <- paste0(str2,"','",monitor[k])
    }
  }
  
  str2 <- paste0(str2,"'))\n mcmcConf \n }")
  
  names(MCMCdefs)[1] <- 'Default'
  str1 <- "{\n
    mcmcConf <- configureMCMC(Rmodel)\n"
  
  for (i in 1:n) {
    str1 <- paste0(str1,"\n mcmcConf$removeSamplers('",monitor[i],"')\n")
  }
  for (i in 1:n) {
    if (monitor[i] == targetNames[[1]]) {
      for (l in 1:length(CandidateSamplerList$target)) {
        if (CandidateSamplerList$target[[l]]$type == 'sampler_RW_block') {
          str1 <- paste0(str1,"\n mcmcConf$addSampler(target = c('", GroupLM[1])
          if (length(GroupLM) > 1) {
            for (j in 2:length(GroupLM))
              str1 <- paste0(str1,"','", GroupLM[j])
          }
          str1 <-
            paste0(
              str1,"'), type = sampler_record_wrapper, control = list(sampler_function = '",CandidateSamplerList$target[[l]]$type,"', control = list("
            )
        } else if (CandidateSamplerList$target[[l]]$type == 'sampler_conjugate') {
          conjInfo <- Rmodel$checkConjugacy2(monitor[i])[[monitor[i]]]
          if (!is.null(conjInfo)) {
            str1 <-
              paste0(
                str1, "\n mcmcConf$addSampler(target = '",monitor[i],"', type = sampler_conjugate_wrapper, control=list(), name = 'conjugate') \n"
              )
          }
        } else{
          str1 <-
            paste0(
              str1,"\n mcmcConf$addSampler(target = '",monitor[i],"', type = sampler_record_wrapper, control = list(sampler_function = '",CandidateSamplerList$target[[l]]$type,"', control = list("
            )
          
        }
        if (CandidateSamplerList$target[[l]]$type != 'sampler_conjugate'){
        m = length(CandidateSamplerList$target[[l]]$control)
        if (m == 1) {
          str1 <-
            paste0(
              str1, names(CandidateSamplerList$target[[l]]$control[1]),"=",CandidateSamplerList$target[[l]]$control[1]
            )
        } else if (m > 1) {
          str1 <-
            paste0(
              str1, names(CandidateSamplerList$target[[l]]$control[1]),"=",CandidateSamplerList$target[[l]]$control[1]
            )
          for (j in 1:(m - 1)) {
            str1 <-
              paste0(
                str1,",", names(CandidateSamplerList$target[[l]]$control[j + 1]),"=",CandidateSamplerList$target[[l]]$control[j +
                                                                                                                                1]
              )
            
          }
          
        }
        
          str1 <-
          paste0(str1,")), name = '",CandidateSamplerList$target[[l]]$name,"')\n")
        }
      }
    }
    else{
      conjInfo <- NULL
      if (regexpr('conjugate', DefaultSamplerList[[i]]$name) > 0) {
        conjInfo <- Rmodel$checkConjugacy2(monitor[i])[[monitor[i]]]
        if (!is.null(conjInfo)) {
          str1 <-
              paste0(
                str1, "\n mcmcConf$addSampler(target = '",monitor[i],"', type = sampler_conjugate_wrapper, control=list(), name = 'conjugate') \n"
              )
        }
        
      } else if (DefaultSamplerList[[i]]$type == 'sampler_RW_block') {
        str1 <-
          paste0(
            str1,"\n mcmcConf$addSampler(target = c('",monitor[i],"','", monitor[-i][1],"'), type = sampler_record_wrapper, control = list(sampler_function = 'sampler_RW_block', control = list("
          )
      } else if (DefaultSamplerList[[i]]$type == 'RW') {
        str1 <-
          paste0(
            str1,"\n mcmcConf$addSampler(target = '",monitor[i],"', type = sampler_record_wrapper, control = list(sampler_function = 'sampler_RW', control = list("
          )
      } else{
        str1 <-
          paste0(
            str1,"\n mcmcConf$addSampler(target = '",monitor[i],"', type = sampler_record_wrapper, control = list(sampler_function = '",DefaultSamplerList[[i]]$type,"', control = list("
          )
      }
      if (is.null(conjInfo)){
      m = length(DefaultSamplerList[[i]]$control)
      if (m == 1) {
        str1 <-
          paste0(
            str1, names(DefaultSamplerList[[i]]$control[1]),"=",DefaultSamplerList[[i]]$control[1]
          )
      } else if (m > 1) {
        str1 <-
          paste0(
            str1, names(DefaultSamplerList[[i]]$control[1]),"=",DefaultSamplerList[[i]]$control[1]
          )
        for (j in 1:(m - 1)) {
          str1 <-
            paste0(
              str1,",", names(DefaultSamplerList[[i]]$control[j + 1]),"=",DefaultSamplerList[[i]]$control[j +
                                                                                                            1]
            )
          
        }
        
      }
      
        str1 <-
        paste0(str1,")), name = '",DefaultSamplerList[[i]]$name,"')\n")
      }
    }
  }
  str1 <- paste0(str1,str2)
  
  MCMCdefs[[1]] = parse(text = str1)
  
  return (MCMCdefs)
}  



GroupOfLeastMixing <- function(Samples, leastMixing){
  empCov <- cov(Samples)
  empCor <<- cov2cor(empCov)
  distMatrix <- as.dist(1 - abs(empCor))
  hTree <- hclust(distMatrix, method = 'complete')
  h=mean(hTree$height)
  cutreeList <- cutree(hTree, h=h)
  uniqueCutreeList <- unique(cutreeList)
  candidateGroupsList <- lapply(uniqueCutreeList, function(ct) {
    return(names(which(cutreeList==ct)))}) 
  for(i in seq_along(candidateGroupsList)){
    if(leastMixing %in% candidateGroupsList[[i]])
      return (candidateGroupsList[[i]])
  }
} 

## Construct MCMCconf from DefaultSamplerList
BuildDefaultConf <- function(DefaultSamplerList, monitor){
  n = length(monitor)
  str2<-"conf$addMonitors(c('"
    if (n==1){
      str2 <-paste0(str2, monitor[1])
    } else if (n>1){
      str2 <-paste0(str2, monitor[1]) 
      for(k in 2:n){
        str2 <-paste0(str2,"','",monitor[k]) 
      }
     }
    
    str2<- paste0(str2,"'))\n")
   
    str1 <- "conf <- configureMCMC(Rmodel)\n"
  
  for (i in 1:n){
    str1<- paste0(str1,"conf$removeSamplers('",monitor[i],"')\n")
  }      
  for (i in 1:n){
    if(regexpr('conjugate', DefaultSamplerList[[i]]$type)>0){
      str1 <-
              paste0(
                str1, "\n conf$addSampler(target = '",monitor[i],"', type = 'sampler_conjugate_wrapper', control=list(), name = '",DefaultSamplerList[[i]]$name,"') \n"
              )
    } else if(DefaultSamplerList[[i]]$type=='sampler_RW_block'){
    str1 <- paste0(str1,"\n conf$addSampler(target = c('", DefaultSamplerList[[i]]$target[1])
          if (length(DefaultSamplerList[[i]]$target) > 1) {
            for (j in 2:length(DefaultSamplerList[[i]]$target))
              str1 <- paste0(str1,"','", DefaultSamplerList[[i]]$target[j])
          }
          str1 <-
            paste0(str1,"'), type = 'sampler_AF_slice', control = list(")
    } else{
      str1<-paste0(str1,"conf$addSampler(target = '",monitor[i],"', type ='", DefaultSamplerList[[i]]$type,"', control = list(")
    }
    
    if(!regexpr('conjugate', DefaultSamplerList[[i]]$type)>0){
    m = length(DefaultSamplerList[[i]]$control)
    if (m==1){
      str1 <-paste0(str1, names(DefaultSamplerList[[i]]$control[1]),"=",DefaultSamplerList[[i]]$control[1])
    } else if (m>1){
      str1 <-paste0(str1, names(DefaultSamplerList[[i]]$control[1]),"=",DefaultSamplerList[[i]]$control[1])
      for(j in 1:(m-1)){
        str1 <-paste0(str1,",", names(DefaultSamplerList[[i]]$control[j+1]),"=",DefaultSamplerList[[i]]$control[j+1])
        
      }
      
    }
    
      str1 <-paste0(str1,"), name = '",DefaultSamplerList[[i]]$name,"')\n")
    }
    
  }
  str1<-paste0(str1,str2)
  
  conf =parse(text=str1) 

  return (conf)
}


ImproveMixing1 <- function(code, constants, data, inits, niter, burnin, tuning, monitors, makePlot, calculateEfficiency, setSeed, DefaultSamplerList, CandidateSamplerList, verbose) {
  
  bestEfficiency = c(a=0)


  repeat {
  
    
    #### Build and run the default sampler.
    conf <-BuildDefaultConf(DefaultSamplerList, monitor= monitor)
    eval(conf)
    Rmcmc <- buildMCMC(conf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    
    timeResult <- system.time({ Cmcmc$run(niter) })
    print(timeResult[3])
    
    #### Get the least mixing index.
    Samples <- as.matrix(Cmcmc$mvSamples)
    output<-apply(Samples, 2, effectiveSize)/timeResult[3]
    lindex <- which.min(output)
    
    leastMixing <- names(output[lindex])
    if(verbose){
      print(output)
      print("Current efficiency:")
      print(bestEfficiency)
      print("New least efficiency:")
      print(output[lindex])
      print("Least mixing variable:")
      print(leastMixing)
    }
    
      
    if(output[lindex] < bestEfficiency){
      print("Can not improve. Stop iteration.")
      return (list(DefaultSamplerList, monitor))
      
    } else{
      bestEfficiency <- output[lindex]
    }
    leastMixing <- names(output[lindex])
    if(verbose){
      print(output)
      print("Least mixing variable:")
      print(leastMixing)
    }
    
      leastIndex = which(monitor==leastMixing)
      ### Cluster the least mixing variable using distance matrix
      GroupLM<-GroupOfLeastMixing(Samples, leastMixing)
      #### Choose the target as the least mixing variable, generate a list of cadidate samplers.
      targetNames= list(monitor[leastIndex])
      print("Block includes:")
      if(length(GroupLM)>1){
        print(GroupLM)
      } else { ## Just block anything to test
        print(c(monitor[leastIndex],monitor[-leastIndex][1])) 
      }  
      MCMCdefs <-BuildCombinedConf(DefaultSamplerList=DefaultSamplerList, CandidateSamplerList=CandidateSamplerList, targetNames=targetNames, GroupLM=GroupLM, monitor= monitor)
      MCMCs = names(MCMCdefs)
      #### Then we run codess to identify the best mixing sampler for the worst mixing variable.
    test1<-MCMC_CODESS(code = code, constants = constants, data=data, inits=inits, MCMCs = MCMCs, MCMCdefs = MCMCdefs, targetNames= targetNames,  niter = niter, burnin = burnin, tuning=tuning, monitors = monitor,
    makePlot = makePlot, calculateEfficiency = calculateEfficiency, setSeed = setSeed)
  
    #### Find the best mixing sampler
    if(verbose)
      print(test1$summary)
    output<-test1$summary[[1]]
    n=nrow(output)
    bestIndex <- which.max(output[(length(monitor)+1):n,'efficiency'])
    #### Replace the worst mixing sampler with the best mixing sampler found.
    DefaultSamplerList[[leastIndex]]<-CandidateSamplerList$target[[bestIndex]]
if(CandidateSamplerList$target[[bestIndex]]$type=="sampler_RW_block")
{
  if(length(GroupLM)>1)
    DefaultSamplerList[[leastIndex]]$target=GroupLM
  else
    DefaultSamplerList[[leastIndex]]$target=c(monitor[leastIndex],monitor[-leastIndex][1])
} else {
  DefaultSamplerList[[leastIndex]]$target=monitor[leastIndex]

} 

    if(verbose){
      print("The best mixing index:")
      print(bestIndex)
      print("The sampler of the least mixing now is:")
      print(DefaultSamplerList[[leastIndex]])
      print("End of iteration")
      print("###############################################################")
    }
  
    
    
  }  

}


ImproveMixing <- function(code, constants, data, inits, niter, burnin, tuning, monitors, makePlot, calculateEfficiency, setSeed, DefaultSamplerList, CandidateSamplerList, verbose) {
  
  bestEfficiency = c(a=0)
  j = 0

  initialMCMCconf <<- configureMCMC(Rmodel)
            nInitialSamplers <- length(initialMCMCconf$samplerConfs)
            initialMCMCconf$addSampler(type = 'RW',       target = scalarNodeVector[1], print=FALSE)  ## add one RW sampler
            initialMCMCconf$addSampler(type = 'RW_block', target = scalarNodeVector[1], print=FALSE)  ## add one RW_block sampler
            initialMCMCconf$addSampler(type = 'sampler_conjugate', target = scalarNodeVector[1], print=FALSE)  ## add one conjugate sampler
            initialMCMCconf$addSampler(type = 'slice_sampler', target = scalarNodeVector[1], print=FALSE)  ## add one conjugate sampler
            initialMCMCconf$addMonitors(monitors, print=FALSE)
            RinitialMCMC <- buildMCMC(initialMCMCconf)
            Cmodel <<- compileNimble(Rmodel)
            CinitialMCMC <- compileNimble(RinitialMCMC, project = Rmodel)   ## (new version) yes, we need this compileNimble call -- this is the whole point!
            initialMCMCconf$setSamplers(1:nInitialSamplers, print=FALSE)  ## important for new version: removes all news samplers added to initial MCMC conf
            
            
            

  repeat {
  
    j= j+1
    #### Build and run the default sampler.

conf <- configureMCMC(oldConf = initialMCMCconf)  ## new version
            conf$setSamplers()  ## new version -- removes all the samplers from initalMCMCconf
            
            for (i in 1:nInitialSamplers){
    if(regexpr('conjugate', DefaultSamplerList[[i]]$type)>0){
      conf$addSampler(target = monitor[i],type = 'sampler_conjugate', control=list(), name = DefaultSamplerList[[i]]$name)
    } else if(DefaultSamplerList[[i]]$type=='sampler_RW_block'){
      conf$addSampler(target = DefaultSamplerList[[i]]$target, type = sampler_RW_block, control = DefaultSamplerList[[i]]$control, name = DefaultSamplerList[[i]]$name)
    } else{
      conf$addSampler(target = monitor[i], type =DefaultSamplerList[[i]]$type, control = DefaultSamplerList[[i]]$control, name = DefaultSamplerList[[i]]$name)
    }
    

    
  }
  conf$addMonitors(monitors, print=FALSE)
  


    Rmcmc <- buildMCMC(conf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    
    timeResult <- system.time({ Cmcmc$run(niter) })
    print(timeResult[3])
    
    #### Get the least mixing index.
    Samples <- as.matrix(Cmcmc$mvSamples)
    output<-apply(Samples, 2, effectiveSize)/timeResult[3]
    lindex <- which.min(output)
    
    
    leastMixing <- names(output[lindex])
    if(verbose){
      print(output)
      print("Current efficiency:")
      print(bestEfficiency)
      print("New least efficiency:")
      print(output[lindex])
      print("Least mixing variable:")
      print(leastMixing)
    }
    
      
    if(j>2){
    #(output[lindex] < bestEfficiency){
      print("Can not improve. Stop iteration.")
      return (list(DefaultSamplerList, monitor))
      break
    } else{
      bestEfficiency <- output[lindex]
    }
    leastMixing <- names(output[lindex])
    if(verbose){
      print(output)
      print("Least mixing variable:")
      print(leastMixing)
    }
    
      leastIndex = which(monitor==leastMixing)
      ### Cluster the least mixing variable using distance matrix
      GroupLM<-GroupOfLeastMixing(Samples, leastMixing)
      #### Choose the target as the least mixing variable, generate a list of cadidate samplers.
      targetNames= list(monitor[leastIndex])
      print("Block includes:")
      if(length(GroupLM)>1){
        print(GroupLM)
      } else { ## Just block anything to test
        print(c(monitor[leastIndex],monitor[-leastIndex][1])) 
      }  
      MCMCdefs <-BuildCombinedConf(DefaultSamplerList=DefaultSamplerList, CandidateSamplerList=CandidateSamplerList, targetNames=targetNames, GroupLM=GroupLM, monitor= monitor)
      MCMCs = names(MCMCdefs)
      #### Then we run codess to identify the best mixing sampler for the worst mixing variable.
    test1<-MCMC_CODESS(code = code, constants = constants, data=data, inits=inits, MCMCs = MCMCs, MCMCdefs = MCMCdefs, targetNames= targetNames,  niter = niter, burnin = burnin, tuning=tuning, monitors = monitor,
    makePlot = makePlot, calculateEfficiency = calculateEfficiency, setSeed = setSeed)
  
    #### Find the best mixing sampler
    if(verbose)
      print(test1$summary)
    output<-test1$summary[[1]]
    n=nrow(output)
    bestIndex <- which.max(output[(length(monitor)+1):n,'efficiency'])
    #### Replace the worst mixing sampler with the best mixing sampler found.
    DefaultSamplerList[[leastIndex]]<-CandidateSamplerList$target[[bestIndex]]
if(CandidateSamplerList$target[[bestIndex]]$type=="sampler_RW_block")
{
  if(length(GroupLM)>1)
    DefaultSamplerList[[leastIndex]]$target=GroupLM
  else
    DefaultSamplerList[[leastIndex]]$target=c(monitor[leastIndex],monitor[-leastIndex][1])
} else {
  DefaultSamplerList[[leastIndex]]$target=monitor[leastIndex]

} 

    if(verbose){
      print("The best mixing index:")
      print(bestIndex)
      print("The sampler of the least mixing now is:")
      print(DefaultSamplerList[[leastIndex]])
      print("End of iteration")
      print("###############################################################")
    }
  
    
    
  }   

}

autoCodess <- function(code, constants=list(), data=list(), inits=list(),
                      niter = 10000,
                      CandidateGroupList = CandidateGroupList,
                      setSeed0 = TRUE,
                      verbose = FALSE,
                      saveSamples = FALSE,
                      round = TRUE ) {
    ab <- autoCodessClass_oldClass(code, constants, data, inits,
                         control = list(niter=niter, setSeed0=setSeed0, verbose=verbose, saveSamples=saveSamples))
 #   if(!'auto' %in% run) run <- c(run, 'auto')  ## always use 'autoCodess' routine
    ab$run(CandidateGroupList)
#     lastAutoInd <- max(grep('^auto', ab$naming))   ## index of final 'auto' iteration
#     lastAutoGrouping <- ab$grouping[[lastAutoInd]]  ## grouping of final 'auto' iteration
#     nonTrivialGroups <- lastAutoGrouping[unlist(lapply(lastAutoGrouping, function(x) length(x)>1))]
#     abList <- list(ab)
#     names(abList)[1] <- 'model'
#     df <- createDFfromABlist(abList, niter)
#     dfmin <- reduceDF(df, round = round)
#     cat('\nAuto-Codess summary:\n')
#     print(dfmin)
#     if(length(nonTrivialGroups) > 0) {
#         cat('\nAuto-Codess converged on the node groupings:\n')
#         for(i in seq_along(nonTrivialGroups)) {
#             group <- nonTrivialGroups[[i]]
#             cat(paste0('[', i, '] '))
#             cat(paste0(group, collapse = ', '))
#             cat('\n')
#         }
#     } else cat('\nAuto-Codess converged on all scalar (univariate) sampling\n')
#     cat('\n')
#     ## create a new Rmodel and conf with the autoBlock groupings:
#     Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
#     conf <- configureMCMC(Rmodel, nodes = NULL)
#     for(nodeGroup in lastAutoGrouping) addSamplerToConf(Rmodel, conf, nodeGroup)
    
       
    retList <- list(ess=ab$ess, efficiency=ab$essPT)
    if(saveSamples) retList$samples <- ab$samples
    return(invisible(retList))
}

autoCodessModel_oldClass <- setRefClass(
    Class = 'autoCodessModel_oldClass',
    fields = list(
        code = 'ANY',
        constants = 'list',
        data = 'list',
        inits = 'list',
        md = 'ANY',
        Rmodel = 'ANY',
        Cmodel = 'ANY',
        scalarNodeVector = 'character',
        latentNodeVector = 'character',
        nodeGroupScalars = 'list',
        nodeGroupAllBlocked = 'list',
        monitorsVector = 'character',
        initialMCMCconf = 'ANY'
        ),
    methods = list(
        initialize = function(code, constants, data, inits) {
            library(nimble)
            code <<- code
            constants <<- if(missing(constants)) list() else constants
            data <<- if(missing(data)) list() else data
            inits <<- if(missing(inits)) list() else inits
            md <<- nimbleModel(code=code, constants=constants, returnDef=TRUE)
            Rmodel <<- md$newModel(data=data, inits=inits)
            scalarNodeVector <<- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, returnScalarComponents=TRUE)
            latentNodeVector <<- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, latent=TRUE)
            nodeGroupScalars <<- lapply(scalarNodeVector, function(x) x)
            nodeGroupAllBlocked <<- list(scalarNodeVector)
            stochNodeVector <- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, returnScalarComponents=FALSE)
            monitorsVector <<- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE)
        },
        ## here is where the initial MCMC conf is created, for re-use -- for new version
        createInitialMCMCconf = function() {
            initialMCMCconf <<- configureMCMC(Rmodel)
            nInitialSamplers <- length(initialMCMCconf$samplerConfs)
            initialMCMCconf$addSampler(type = 'sampler_RW',       target = latentNodeVector[1], print=FALSE)  ## add one RW sampler
            initialMCMCconf$addSampler(type = 'sampler_RW_block', target = latentNodeVector[1], print=FALSE)  ## add one RW_block sampler
            initialMCMCconf$addSampler(type = 'sampler_conjugate', target = latentNodeVector[1], print=FALSE)  ## add one RW_block sampler
            initialMCMCconf$addSampler(type = 'sampler_slice', target = latentNodeVector[1], print=FALSE)  ## add one RW_block sampler
            initialMCMCconf$addSampler(type = 'sampler_AF_slice', target = latentNodeVector[1:2], control = list(sliceWidths = c(.5,.5), sliceFactorBurnInIters = 500,sliceFactorAdaptInterval = 50, sliceSliceAdaptIters = 50), print=FALSE)  ## add one RW_block sampler
            
            #addCustomizedSamplersToInitialMCMCconf()
            initialMCMCconf$addMonitors(monitorsVector, print=FALSE)
            RinitialMCMC <- buildMCMC(initialMCMCconf)
            Cmodel <<- compileNimble(Rmodel)
            CinitialMCMC <- compileNimble(RinitialMCMC, project = Rmodel)   ## (new version) yes, we need this compileNimble call -- this is the whole point!
            initialMCMCconf$setSamplers(1:nInitialSamplers, print=FALSE)  ## important for new version: removes all news samplers added to initial MCMC conf
        },
        addCustomizedSamplersToInitialMCMCconf = function(runListCode) {
            if(is.list(runListCode)) { lapply(runListCode, function(el) addCustomizedSamplersToInitialMCMCconf(el)); return() }
            if(is.call(runListCode)) {
                if(is.call(runListCode[[1]]) && length(runListCode[[1]])==3 && runListCode[[1]][[3]]=='addSampler') {
                    runListCode[[1]][[2]] <- as.name('initialMCMCconf')
                    eval(substitute(RUNLISTCODE, list(RUNLISTCODE=runListCode)))
                    return()
                }
                lapply(runListCode, function(el) addCustomizedSamplersToInitialMCMCconf(el))
                return()
            }
        },
        createGroups = function(listOfBlocks = list()) {
            listOfBlocks <- lapply(listOfBlocks, function(blk) Rmodel$expandNodeNames(blk, returnScalarComponents=TRUE))
            nodes <- scalarNodeVector
            nodes <- setdiff(nodes, unlist(listOfBlocks))
            nodeList <- lapply(nodes, function(x) x)
            for(ng in listOfBlocks) nodeList[[length(nodeList)+1]] <- ng
            return(nodeList)
        },
        newModel = function() {
            newRmodel <- md$newModel(data=data, inits=inits)
            return(newRmodel)
        }
    )
)



autoCodessParamDefaults <- function() {
    list(
        cutree_heights = c(0.5),#seq(0, 0.2, by=0.1),
        makePlots = FALSE,
        niter = 50000,
        saveSamples = FALSE,
        setSeed0 = TRUE,
        verbose = TRUE
        )
}



autoCodessClass_oldClass <- setRefClass(
    Class = 'autoCodessClass_oldClass',
    fields = list(
        ## special
        abModel = 'ANY',
        it = 'numeric',
        ## overall control
        cutree_heights = 'numeric',
        makePlots = 'logical',
        niter = 'numeric',
        saveSamples = 'logical',
        setSeed0 = 'logical',
        verbose = 'logical',
        ## persistant lists of historical data
        keepTrack = 'list',
        naming = 'list',
        candidateGroups = 'list',
        grouping = 'list',
        groupSizes = 'list',
        groupIDs = 'list',
        samplers = 'list',
        Cmcmcs = 'list',
        timing = 'list',
        samples = 'list',
        means = 'list',
        sds = 'list',
        ess = 'list',
        essPT = 'list',
        burnedSamples = 'list',
        empCov = 'list',
        empCor = 'list',
        distMatrix = 'list',
        LeastIndex = 'list',
        hTree = 'list',
        keepTrackTemp = 'matrix'
        ),
    methods = list(
        initialize = function(code, constants=list(), data=list(), inits=list(), control=list(),DefaultSamplerList=list(), CandidateSamplerList=list()) {
            library(lattice)
            library(coda)
            library(nimble)
            abModel <<- autoCodessModel_oldClass(code=code, constants=constants, data=data, inits=inits)
            defaultsList <- autoCodessParamDefaults()
            for(i in seq_along(defaultsList)) if(is.null(control[[names(defaultsList)[i]]])) control[[names(defaultsList)[i]]] <- defaultsList[[i]]
            for(i in seq_along(control)) eval(substitute(verbose <<- VALUE, list(verbose=as.name(names(control)[i]), VALUE=control[[i]])))
            it <<- 0
            keepTrackTemp <<-matrix(0, ncol=length(CandidateSamplerList)+1, nrow= length(DefaultSamplerList))
	          colnames(keepTrackTemp) <<-c(names(CandidateSamplerList),'Include')
	          rownames(keepTrackTemp) <<-names(DefaultSamplerList)
	          
        },
        run = function(candidateGroupsList) {
            abModel$createInitialMCMCconf()  ## here is where the initial MCMC conf is created, for re-use -- for new version
            
	          oldConf = abModel$initialMCMCconf
	    
	    
	          autoIt <- 1
                                     
            confList <- list(createConfFromGroups(DefaultSamplerList))
              runConfListAndSaveBest(confList, paste0('auto',i), auto=TRUE)
                                         
                                 
            names(candidateGroups) <<- naming
            names(grouping) <<- naming
            names(groupSizes) <<- naming
            names(groupIDs) <<- naming
            names(samplers) <<- naming
            names(Cmcmcs) <<- naming
            names(timing) <<- naming
            if(saveSamples) names(samples) <<- naming
            names(means) <<- naming
            names(sds) <<- naming
            names(ess) <<- naming
            names(essPT) <<- naming
        },
determineCandidateGroupsFromCurrentSample = function() {
            cutreeList <- lapply(cutree_heights, function(height) cutree(hTree[[it]], h = height))
            names(cutreeList) <- paste0('cut', cutree_heights)
            uniqueCutreeList <- unique(cutreeList)
            
            candidateGroupsList <- lapply(uniqueCutreeList, function(ct) determineGroupsFromCutree(ct))
            return(candidateGroupsList)
        },
        determineGroupsFromCutree = function(ct) {
            return(lapply(unique(ct), function(x) names(ct)[ct==x]))
        },
        runConfListAndSaveBest = function(confList, name, auto=FALSE) {
            RmcmcList <- timingList <- samplesList <- meansList <- sdsList <- essList <- essPTList <- essPTminList <- list()
            for(i in seq_along(confList)) {
                checkOverMCMCconf(confList[[i]])
                ##confList[[i]]$addMonitors(abModel$monitorsVector, print=FALSE) ## original version
                RmcmcList[[i]] <- buildMCMC(confList[[i]])
            }
            ##toCompileList <- c(list(Rmodel), RmcmcList) ## original version
            ##compiledList <- compileNimble(toCompileList)
            ##Cmodel <- compiledList[[1]]
            ##CmcmcList <- compiledList[-1]
            Cmodel <- abModel$Cmodel
            CmcmcList <- compileNimble(RmcmcList, project = abModel$Rmodel)
            if(!is.list(CmcmcList)) CmcmcList <- list(CmcmcList)  ## make sure compileNimble() returns a list...
            for(i in seq_along(CmcmcList)) {
                Cmodel$setInits(abModel$inits)
                if(setSeed0) set.seed(0)
                timingList[[i]] <- as.numeric(system.time(CmcmcList[[i]]$run(niter))[3])
                ## slight hack here, to remove samples of any deterministic nodes...
                samplesTEMP <- as.matrix(CmcmcList[[i]]$mvSamples)
                namesToKeep <- setdiff(dimnames(samplesTEMP)[[2]], abModel$Rmodel$getNodeNames(determOnly=TRUE, returnScalarComponents=TRUE))
                samplesList[[i]] <- samplesTEMP[, namesToKeep]
                ## end of slight hack...
                meansList[[i]] <- apply(samplesList[[i]], 2, mean)
                sdsList[[i]]   <- apply(samplesList[[i]], 2, sd)
                essList[[i]]   <- apply(samplesList[[i]], 2, effectiveSize)
                if(!saveSamples) samplesList[[i]] <- NA
                essPTList[[i]] <- essList[[i]] / timingList[[i]]
                essPTminList[[i]] <- sort(essPTList[[i]])[1]
            }
            bestInd <- as.numeric(which(unlist(essPTminList) == max(unlist(essPTminList))))
            if(!is.null(names(confList))) name <- paste0(name, '-', names(confList)[bestInd])
            it <<- it + 1
            naming[[it]] <<- name
            candidateGroups[[it]] <<- lapply(confList, function(conf) determineGroupsFromConf(conf))
            grouping[[it]] <<- candidateGroups[[it]][[bestInd]]
            groupSizes[[it]] <<- determineNodeGroupSizesFromGroups(grouping[[it]])
            groupIDs[[it]] <<- determineNodeGroupIDsFromGroups(grouping[[it]])
            samplers[[it]] <<- determineSamplersFromGroupsAndConf(grouping[[it]], confList[[bestInd]])
            Cmcmcs[[it]] <<- CmcmcList[[bestInd]]
            timing[[it]] <<- timingList[[bestInd]]
            samples[[it]] <<- samplesList[[bestInd]]
            means[[it]] <<- meansList[[bestInd]]
            sds[[it]] <<- sdsList[[bestInd]]
            ess[[it]] <<- essList[[bestInd]]
            essPT[[it]] <<- sort(essPTList[[bestInd]])
            if(TRUE) {
                ## slight hack here, to remove samples of any deterministic nodes...
                samplesTEMP <- as.matrix(CmcmcList[[bestInd]]$mvSamples)
                namesToKeep <- setdiff(dimnames(samplesTEMP)[[2]], abModel$Rmodel$getNodeNames(determOnly=TRUE, returnScalarComponents=TRUE))
                burnedSamples[[it]] <<- samplesTEMP[(floor(niter/2)+1):niter, namesToKeep]
                ## end of slight hack...
                empCov[[it]] <<- cov(burnedSamples[[it]])
                empCor[[it]] <<- cov2cor(empCov[[it]])
                distMatrix[[it]] <<- as.dist(1 - abs(empCor[[it]]))
                hTree[[it]] <<- hclust(distMatrix[[it]], method = 'complete')
                print(hTree[[it]])
            }
            if(!saveSamples) burnedSamples[[it]] <<- NA
            if(verbose) printCurrent(name, confList[[bestInd]])
            if(makePlots && auto) makeCurrentPlots(name)
        },
        determineGroupsFromConf = function(conf) {
            groups <- list()
            for(ss in conf$samplerConfs) {
                if(ss$name == 'crossLevel') {
                    topNodes <- ss$target
                    lowNodes <- conf$model$getDependencies(topNodes, self=FALSE, stochOnly=TRUE, includeData=FALSE)
                    nodes <- c(topNodes, lowNodes)
                } else {
                    nodes <- ss$target
                }
                groups[[length(groups)+1]] <- conf$model$expandNodeNames(nodes, returnScalarComponents=TRUE)
            }
            return(groups)
        },
        determineNodeGroupSizesFromGroups = function(groups) {
            groupSizeVector <- numeric(0)
            for(gp in groups) for(node in gp) groupSizeVector[[node]] <- length(gp)
            return(groupSizeVector)
        },
        determineNodeGroupIDsFromGroups = function(groups) {
            groupIDvector <- numeric(0)
            for(i in seq_along(groups)) for(node in groups[[i]]) groupIDvector[[node]] <- i
            return(groupIDvector)
        },
        determineSamplersFromGroupsAndConf = function(groups, conf) {
            samplerConfs <- conf$samplerConfs
            if(length(groups) != length(samplerConfs)) stop('something wrong')
            samplerVector <- character(0)
            for(i in seq_along(groups)) for(node in groups[[i]]) samplerVector[[node]] <- samplerConfs[[i]]$name
            return(samplerVector)
        },
        createConfFromGroups = function(groups) {
            ##conf <- configureMCMC(Rmodel, nodes=NULL, monitors=character(0)) ## original version
            conf <- configureMCMC(oldConf = abModel$initialMCMCconf)  ## new version
            conf$setSamplers()  ## new version -- removes all the samplers from initalMCMCconf
            
            n = length(groups)
            for (i in 1:n){
              if(regexpr('conjugate', groups[[i]]$type)>0){
                if(i<0){
                  conf$addSampler(target = groups[[i]]$target,type = 'sampler_conjugate', control=list())
                } else {
                  conf$addSampler(target = groups[[i]]$target,type = 'sampler_RW')
                }
      
              } else if(groups[[i]]$type=='sampler_RW_block' | length(groups[[i]]$target)>1 ){
      
                conf$addSampler(target = groups[[i]]$target, type = 'sampler_AF_slice', control = list(sliceWidths = rep(.5, length(groups[[i]]$target) ), sliceFactorBurnInIters = 50,sliceFactorAdaptInterval = 50, sliceSliceAdaptIters = 50))
              } else{
                conf$addSampler(target = groups[[i]]$target, type = groups[[i]]$type , control = groups[[i]]$control)
              }
            }
            
            
            return(conf)
        },
        checkOverMCMCconf = function(conf) {
            warn <- FALSE
            for(ss in conf$samplerConfs) {
                ## if(ss$name == 'end') {
                ##     msg <- 'using \'end\' sampler may lead to results we don\'t want'
                ##     cat(paste0('\nWARNING: ', msg, '\n\n')); warning(msg)
                ## }
                if(grepl('^conjugate_', ss$name) && nimbleOptions('verifyConjugatePosteriors')) {
                    ##msg <- 'conjugate sampler running slow due to checking the posterior'
                    ##cat(paste0('\nWARNING: ', msg, '\n\n')); warning(msg)
                    warn <- TRUE
                }
            }
            if(warn) {
                msg <- 'Conjugate sampler functions in \'default\' conf are running slow due to verifying the posterior;\nThis behaviour can be changed using a NIMBLE package option.'
                warning(msg, call. = FALSE)
            }
        },
        printCurrent = function(name, conf) {
            cat(paste0('\n################################\nbegin iteration ', it, ': ', name, '\n################################\n'))
            if(length(candidateGroups[[it]]) > 1) { cat('\ncandidate groups:\n'); cg<-candidateGroups[[it]]; for(i in seq_along(cg)) { cat(paste0('\n',names(cg)[i],':\n')); printGrouping(cg[[i]]) } }
            cat('\ngroups:\n'); printGrouping(grouping[[it]])
            cat('\nsamplers:\n'); conf$getSamplers()
            cat(paste0('\nMCMC runtime: ', round(timing[[it]], 1), ' seconds\n'))
            cat('\nESS:\n'); print(round(ess[[it]], 0))
            cat('\nESS/time:\n'); print(round(essPT[[it]], 1))
            cat(paste0('\n################################\nend iteration ', it, ': ', name, '\n################################\n\n'))
        },
        makeCurrentPlots = function(name) {
            dev.new()
            if(inherits(try(plot(as.dendrogram(hTree[[it]]), ylim=c(0,1), main=name), silent=TRUE), 'try-error')) dev.off()
        },
        printGrouping = function(g) {
            for(i in seq_along(g)) cat(paste0('[', i, '] ', paste0(g[[i]], collapse=', '), '\n'))
        },
        groupingsEquiv = function(grouping1, grouping2) {
            grouping1 <- lapply(grouping1, sort)
            grouping2 <- lapply(grouping2, sort)
            while(length(grouping1) > 0) {
                grp1 <- grouping1[[1]]
                found <- FALSE
                for(i in seq_along(grouping2)) {
                    grp2 <- grouping2[[i]]
                    if(identical(grp1, grp2)) {
                        found <- TRUE
                        grouping1[1] <- grouping2[i] <- NULL
                        break
                    }
                }
                if(!found) return(FALSE)
            }
            if(length(grouping2) == 0) return(TRUE) else return(FALSE)
        }
        )
)



#' Class \code{MCMCsuiteClass}
#'
#' @aliases MCMCsuiteClass-class
#'
#' @description
#' Objects of this class create, run, and organize output from a suite of MCMC algorithms, all applied to the same model, data, and initial values.
#' This can include WinBUGS, OpenBUGS, JAGS and Stan MCMCs, as well as NIMBLE MCMC algorithms.
#' Trace plots and density plots for the MCMC samples may also be generated and saved.
#'
#' @seealso \link{MCMCsuite}
#' 
#' @author Daniel Turek
#' @export
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#' })
#' output <- MCMC_CODESS(code,
#'                     data = list(x=3),
#'                     inits = list(mu=0),
#'                     niter = 10000,
#'                     monitors = 'mu',
#'                     MCMCs = c('nimble', 'nimble_RW'),
#'                     summaryStats = c('mean', 'sd', 'max', 'function(x) max(abs(x))'),
#'                     makePlot = FALSE)
#' }
#' 
MCMC_CODESSClass <- setRefClass(
  
  Class = 'MCMC_CODESSClass',
  
  fields = list(
    ## set in initialize()
    code = 'ANY',   ## parsed expression for the model code; must be contained in { ... }    --- ORIGINAL ARGUMENT
    constants = 'list',   ## list of the constants (ORIGINAL ARGUMENT)
    data = 'list',   ## list of the data    --- ORIGINAL ARGUMENT
    inits = 'list',  ## named list of initial values used for all MCMC algorithms    --- ORIGINAL ARGUMENT
    constantsAndData = 'list',   ## data list used for WinBUGS, OpenBUGS, JAGS.  is equal to c(constantList, dataList)
    Rmodel = 'ANY',   ## Rmodel object
    
    ## setMonitors()
    monitors = 'character',    ## the original character vector argument to initialize()    --- ORIGINAL ARGUMENT --- SLIGHTLY MODIFIED
    
    monitorVars = 'character',    ## character vector of VARIABLE names of parameters to save
    monitorVars1 = 'character',    ## character vector of VARIABLE names of parameters to save
    
    monitorNodesNIMBLE = 'character',  ## character vector of the monitor node names, with spaces as in nimble: 'y[1, 1]'
    monitorNodesBUGS = 'character',    ## same as monitorNodes, except for WinBUGS and OpenBUGS: no spaces in node names: 'y[1,1]'
    nMonitorNodes = 'numeric',   ## number of monitorNodes
    targetNames = 'list',        
    
    ## set in initialize()
    niter = 'numeric',    ## number of MCMC iterations to run    --- ORIGINAL ARGUMENT
    burnin = 'numeric',   ## burn-in period, the number of initial samples to discard, prior to thinning    --- ORIGINAL ARGUMENT
    thin = 'numeric',   ## thinning interval    --- ORIGINAL ARGUMENT
    tuning = 'list',	    ## tuning parameters for bkde2d
    nkeep = 'numeric',   ## number of samples we'll keep. equal to (niter/thin - burnin)
    burninFraction = 'numeric',  ## fraction of total sampling effort spent on burnin (burnin / (nkeep + burnin))
    
    ## setSummaryStats()
    summaryStats = 'character',    ## character vector of parseable summary statistic functions    --- ORIGINAL ARGUMENT
    calculateEfficiency = 'logical',   ## logical specifying whether to calculate ESS and Efficiency    --- ORIGINAL ARGUMENT
    summaryStatFunctions = 'list',  ## list of the function objects for summary statistics
    summaryStatDimNames = 'character',   ## character vector of the dimension names in output$summary
    nSummaryStats = 'numeric',   ## the number of summary statistics
    
    ## setMCMCs()
    MCMCs = 'character',   ## character vector of the MCMC analyses.  'winbugs', 'openbugs', 'jags', 'stan', or anything else is nimble    --- ORIGINAL ARGUMENT
    winbugsMCMCflag = 'logical',   ## whether 'winbugs' is in MCMCs
    openbugsMCMCflag = 'logical',   ## whether 'openbugs' is in MCMCs
    jagsMCMCflag = 'logical',   ## whether 'jags' is in MCMCs
    stanMCMCflag = 'logical',   ## whether 'stan' is in MCMCs
    nimbleMCMCs = 'character',    ## the names of the remaining (presumably nimble) MCMCs
    nNimbleMCMCs = 'numeric',    ## the number of remaining (nimble) MCMCs
    nimbleMCMCflag = 'logical',   ## a flag indicating whether there are any remaining (nimble) MCMCs
    nMCMCs = 'numeric',   ## the number of MCMC algorithms being ran
    
    ## setMCMCdefs()
    MCMCdefs = 'list',   ## named list of {} expression code blocks, corresponding the setup for nimble MCMCs    --- ORIGINAL ARGUMENT --- SLIGHTLY MODIFIED
    MCMCdefNames = 'character',   ## names of the MCMCdefs list
    
    ## set in initialize()
    winbugs_directory = 'character',    ## directory for WinBUGS program    --- ORIGINAL ARGUMENT
    winbugs_program = 'character',     ## program for WinBUGS    --- ORIGINAL ARGUMENT
    openbugs_directory = 'character',    ## directory for OpenBUGS program    --- ORIGINAL ARGUMENT
    openbugs_program = 'character',     ## program for OpenBUGS    --- ORIGINAL ARGUMENT
    stan_model = 'character',     ## *.stan model file    --- ORIGINAL ARGUMENT
    makePlot = 'logical',    ## whether to generate plots    --- ORIGINAL ARGUMENT
    savePlot = 'logical',   ## whether or not to save plot PDFs    --- ORIGINAL ARGUMENT
    plotName = 'character',     ## name of the file where we save density and trace plots    --- ORIGINAL ARGUMENT
    setSeed = 'logical',   ## whether to setSeed(0) prior to running each algorithm    --- ORIGINAL ARGUMENT
    debug = 'logical',   ## whether to enter browser() before running each algorithm    --- ORIGINAL ARGUMENT
    modelFileName = 'character',     ## name of the text file where we write the model code, set to a fixed value
    
    ## Maps with possible transformations from Stan to WinBUGS/OpenBUGS
    ## e.g. for blocker: StanNameMaps <- list(tau = list(StanSourceName = 'sigmasq_delta', transform = function(x) 1/x)) ## transform can be omitted
    StanNameMaps = 'ANY',
    
    ## set in run()
    Cmodel = 'ANY',   ## compiled Cmodel object
    RmcmcTargetList = 'list',    ## list of the R (nimble) MCMC target
    RmcmcNamesList = 'list',    ## list of the R (nimble) MCMC samplers' names
    
    RmcmcFunctionList = 'list',    ## list of the R (nimble) MCMC functions
    CmcmcFunctionList = 'list',    ## list of the C (nimble) MCMC functions
    output = 'list'   ## list of numeric outputs: samples, summary, timing
    
  ),
  
  methods = list(
    
    initialize = function(
      code,
      constants           = list(),
      data                = list(),
      inits               = list(),
      monitors            = character(),
      niter               = 10000,
      burnin              = 2000,
      thin                = 1,
      tuning		  = list(),
      summaryStats        = c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp'),
      calculateEfficiency = FALSE,
      MCMCs               = 'nimble',
      MCMCdefs            = list(),
      targetNames		= list(),
      
      winbugs_directory   = 'C:/WinBUGS14',
      winbugs_program     = 'WinBUGS',
      openbugs_directory  = 'C:/OpenBUGS323',
      openbugs_program    = 'OpenBUGS',
      stan_model          = '',
      stan_inits          = NULL,
      stan_data           = NULL,
      stanNameMaps        = list(),
      makePlot            = TRUE,
      savePlot            = TRUE,
      plotName            = 'MCMC_CODESS',
      setSeed             = TRUE,
      check               = getNimbleOption('checkModel'),
      debug               = FALSE) {
      
      if(debug) browser()
      code <<- code
      constants <<- constants
      data <<- data
      inits <<- inits
      constantsAndData <<- c(constants, data)
      Rmodel <<- nimbleModel(code=code, constants=constants, data=data, inits=inits, check=check)
      niter <<- niter
      burnin <<- burnin
      thin <<- thin
      tuning <<- tuning
      nkeep <<- floor(niter/thin) - burnin
      burninFraction <<- burnin / (nkeep + burnin)
      setMonitors(monitors)
      targetNames <<- targetNames
      setSummaryStats(summaryStats, calculateEfficiency)
      setMCMCs(MCMCs)
      setMCMCdefs(MCMCdefs)
      winbugs_directory <<- winbugs_directory
      winbugs_program <<- winbugs_program
      openbugs_directory <<- openbugs_directory
      openbugs_program <<- openbugs_program
      stan_model <<- stan_model
      if(is.null(stan_inits)) stan_inits <- gsub('stan$', 'init.R', stan_model)
      if(is.null(stan_data))  stan_data  <- gsub('stan$', 'data.R', stan_model)
      StanNameMaps <<- stanNameMaps
      makePlot <<- makePlot
      savePlot <<- savePlot
      plotName <<- plotName
      setSeed <<- setSeed
      debug <<- debug
      modelFileName <<- 'model.txt'
      
      ## run
      checkMCMCdefNames()
      init_output()
      writeModelFile()
      if(debug)              browser()
      if(winbugsMCMCflag)    run_winbugs()
      if(openbugsMCMCflag)   run_openbugs()
      if(jagsMCMCflag)       run_jags()
      if(stanMCMCflag)       run_stan(stan_data, stan_inits)
      if(nimbleMCMCflag)     run_nimble()
      unlink(modelFileName)
      if(makePlot)           generate_plots()
    },
    
    setMonitors = function(newMonitors) {
      if(length(newMonitors) == 0) newMonitors <- Rmodel$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
      newMonitors <- Rmodel$expandNodeNames(newMonitors, returnScalarComponents = TRUE)
      dataFlags <- unlist(lapply(newMonitors, function(mon) eval(parse(text=mon, keep.source=FALSE)[[1]], envir=Rmodel$isDataEnv)))
      newMonitors <- newMonitors[!dataFlags]
      monitors <<- newMonitors
      monitorVars <<- unique(removeIndexing(monitors))
      monitorVars1 <<- unique(monitors)
      monitorNodesNIMBLE <<- monitors
      monitorNodesBUGS <<- gsub(' ', '', monitorNodesNIMBLE)
      nMonitorNodes <<- length(monitorNodesNIMBLE)
    },
    
    setSummaryStats = function(summaryStats_arg, calculateEfficiency) {
      calculateEfficiency <<- calculateEfficiency
      if(calculateEfficiency) {
        n <- length
        ess <- effectiveSize
        efficiency <- function(x) return(0)   ## placeholder; calculation done in addToOutput()
        codamethod <- function(x) return(1)
        summaryStats_arg <- c(summaryStats_arg, 'n', 'ess', 'efficiency','codamethod')
      }
      summaryStats <<- summaryStats_arg
      CI95_low <- function(x) quantile(x, probs = 0.025)
      CI95_upp <- function(x) quantile(x, probs = 0.975)
      summaryStatFunctions <<- lapply(summaryStats, function(txt) eval(parse(text=txt)[[1]]))
      summaryStatDimNames <<- gsub('function *\\(.*?\\)', '', summaryStats)
      summaryStatDimNames <<- gsub('^ *', '', summaryStatDimNames)
      nSummaryStats <<- length(summaryStats)
    },
    
    setMCMCs = function(MCMCs) {
      MCMCs <<- unique(MCMCs)
      winbugsMCMCflag <<- 'winbugs' %in% MCMCs
      openbugsMCMCflag <<- 'openbugs' %in% MCMCs
      jagsMCMCflag <<- 'jags' %in% MCMCs
      stanMCMCflag <<- 'stan' %in% MCMCs
      nimbleMCMCs <<- setdiff(MCMCs, c('winbugs', 'openbugs', 'jags', 'stan'))
      nNimbleMCMCs <<- length(nimbleMCMCs)
      nimbleMCMCflag <<- if(nNimbleMCMCs > 0) TRUE else FALSE
      nMCMCs <<- length(MCMCs)
    },
    
    setMCMCdefs = function(newMCMCdefs) {
      MCMCdefs <<- list(nimble        = quote(configureMCMC(Rmodel)),
                        nimble_noConj = quote(configureMCMC(Rmodel, useConjugacy = FALSE)),
                        nimble_RW     = quote(configureMCMC(Rmodel, onlyRW       = TRUE)),
                        nimble_slice  = quote(configureMCMC(Rmodel, onlySlice    = TRUE)),
                        autoBlock     = quote(configureMCMC(Rmodel, autoBlock    = TRUE)))
      MCMCdefs[names(newMCMCdefs)] <<- newMCMCdefs
      MCMCdefNames <<- names(MCMCdefs)
    },
    
    init_output = function() {
      samples <- list()
      summary <- list()
      monitor <-list()
      timing <- rep(NA, nMCMCs+1)
      names(timing) <- c(MCMCs, 'nimble_compile')
      if(stanMCMCflag) timing['stan_compile'] <- NA
      runParams <- c(niter = niter, burnin = burnin, thin = thin, nkeep = nkeep, burninFraction = burninFraction) 
      initialOutput <- list(samples=samples, summary=summary, monitor=monitor, timing=timing, runParams = runParams)
      if(calculateEfficiency){
        initialOutput$efficiency <- list()
        initialOutput$codamethod <- list()
      }
      output <<- initialOutput
    },
    
    run_winbugs = function() {
      if(setSeed) set.seed(0)
      if(requireNamespace('R2WinBUGS', quietly = TRUE)) {
        timeResult <- system.time({
          winbugs_out <- R2WinBUGS::bugs(data=constantsAndData, inits=list(inits), parameters.to.save=monitorVars, model.file=modelFileName,
                                         n.chains=1, n.iter=niter, n.burnin=0, n.thin=thin, bugs.directory=winbugs_directory, program=winbugs_program)
        })
        tempArray <- winbugs_out$sims.array[, 1, ]        ## must use sims.array
        samplesArray <- tempArray[(burnin+1):floor(niter/thin), monitorNodesBUGS, drop=FALSE]
        addToOutput('winbugs', samplesArray, timeResult)
      } else warning("run_winbugs: R2WinBUGS package is required for 'winbugs' option.")
    },
    
    run_openbugs = function() {
      if(requireNamespace('R2WinBUGS', quietly = TRUE)) {
        if(setSeed) set.seed(0)
        timeResult <- system.time({
          openbugs_out <- R2WinBUGS::bugs(data=constantsAndData, inits=list(inits), parameters.to.save=monitorVars, model.file=modelFileName,
                                          n.chains=1, n.iter=niter, n.burnin=0, n.thin=thin, bugs.directory=openbugs_directory, program=openbugs_program)
        })
        tempArray <- openbugs_out$sims.array[, 1, ]        ## must use sims.array
        samplesArray <- tempArray[(burnin+1):floor(niter/thin), monitorNodesBUGS, drop=FALSE]
        addToOutput('openbugs', samplesArray, timeResult)
      } else warning("run_openbugs: R2WinBUGS package is required for 'openbugs' option.")
    },
    
    run_jags = function() {
      if(setSeed) set.seed(0)
      if(requireNamespace('rjags', quietly = TRUE)) {
        jags_mod <- rjags::jags.model(file=modelFileName, data=constantsAndData, inits=inits, n.chains=1, quiet=FALSE)
        timeResult <- system.time({
          jags_out <- rjags::coda.samples(model=jags_mod, variable.names=monitorVars, n.iter=niter, thin=thin)
        })
        samplesArray <- jags_out[[1]][(burnin+1):floor(niter/thin), monitorNodesBUGS, drop=FALSE]
        addToOutput('jags', samplesArray, timeResult)
      } else warning("run_jags: rjags package is required for 'jags' option.")
    },
    
    run_stan = function(dataFile, initFile) {
      if(setSeed) set.seed(0)
      if(require('rstan', quietly = TRUE)) {
        if(stan_model == '') stop('must provide \'stan_model\' argument to run Stan MCMC')
        ##            dataFile <- gsub('stan$', 'data.R', stan_model)
        ##            initFile <- gsub('stan$', 'init.R', stan_model)
        if(!is.list(dataFile)) 
          constantsAndDataStan <- fileToList(dataFile)
        else
          constantsAndDataStan <- dataFile
        
        if(!is.list(initFile)) {
          if(file.exists(initFile))
            initsStan <- fileToList(initFile)
          else
            initsStan <- NULL
        } else
          initsStan <- initFile
        
        
        timeResult <- system.time(stan_mod <- rstan::stan_model(file = stan_model))
        addTimeResult('stan_compile', timeResult)
        
        if(is.null(initsStan)) {
          ## missing model.init.R file (stan inits file)
          timeResult <- system.time(stan_out <- rstan::sampling(stan_mod, data=constantsAndDataStan, chains=1, iter=niter, thin=thin))
        } else {
          ## we have the model.init.R file
          ## this one includes inits = ...
          timeResult <- system.time(stan_out <- rstan::sampling(stan_mod, data=constantsAndDataStan, chains=1, iter=niter, thin=thin, init=list(initsStan)))
        }
        
        tempArray <- rstan::extract(stan_out, permuted = FALSE, inc_warmup = TRUE)[, 1, ]
        for(BUGSname in names(StanNameMaps)) {
          iCol <- which(StanNameMaps[[BUGSname]]$StanSourceName == colnames(tempArray))
          if(length(iCol)==1) {
            if(!is.null(StanNameMaps[[BUGSname]]$transform))
              tempArray[,iCol] <- StanNameMaps[[BUGSname]]$transform(tempArray[,iCol])
            colnames(tempArray)[iCol] <- BUGSname
          }
        }
        dimnames(tempArray)[[2]] <- gsub('_', '.', dimnames(tempArray)[[2]])
        if(!all(monitorNodesBUGS %in% dimnames(tempArray)[[2]])) {
          missingNames <- setdiff(monitorNodesBUGS, dimnames(tempArray)[[2]])
          warning(paste0('Stan output is missing values for: ', paste0(missingNames,collapse=', ')))
        }
        samplesArray <- array(0, dim = c(nkeep, length(monitorNodesBUGS)))
        dimnames(samplesArray)[[2]] <- monitorNodesBUGS
        monitorsWeHave <- intersect(monitorNodesBUGS, dimnames(tempArray)[[2]])
        samplesArray[, monitorsWeHave] <- tempArray[(burnin+1):floor(niter/thin), monitorsWeHave, drop=FALSE]
        addToOutput('stan', samplesArray, timeResult)
      } else warning("run_stan: rstan package is required for 'stan' option.")
    },
    
    run_nimble = function() {
      for(iMCMC in seq_along(nimbleMCMCs)) {
        mcmcTag <- nimbleMCMCs[iMCMC]
        mcmcDef <- MCMCdefs[[mcmcTag]]
        mcmcConf <- eval(mcmcDef)
        TargetIndex <-c()
        SamplerIndex <-c()
        for ( i in 1: length(monitorVars1))
        {
          if(monitorVars1[i] == targetNames[[1]]){
            TargetIndex <- c(TargetIndex, mcmcConf$findSamplersOnNodes(monitorVars1[i]))
          }
          
        }
        
        RmcmcTargetList[[iMCMC]] <<- TargetIndex
        Nmonitor <-length(monitorVars1)
        RmcmcNamesList[[iMCMC]] <<-rep(NA, length(TargetIndex)+Nmonitor)
        Samplers <- mcmcConf$getSamplers()
        for (j in 1:(length(TargetIndex)+Nmonitor)){
          if (j <= Nmonitor)
            RmcmcNamesList[[iMCMC]][j] <<-monitorVars1[j]
          else{ 
            if(iMCMC==1){
              RmcmcNamesList[[iMCMC]][j] <<- Samplers[[TargetIndex[j-Nmonitor]]]$name 
                            
            }
            else{
              RmcmcNamesList[[iMCMC]][j] <<- MCMCs[j-Nmonitor] 

            }                
            
          }
        }
        
        mcmcConf$addMonitors(monitorVars, print = FALSE)
        mcmcConf$setThin(thin, print = FALSE)
        RmcmcFunctionList[[mcmcTag]] <<- buildMCMC(mcmcConf)
      }
      timeResult <- system.time({
        Cmodel <<- compileNimble(Rmodel)
        CmcmcFunctionList_temp <- compileNimble(RmcmcFunctionList, project = Rmodel)
        if(nNimbleMCMCs == 1) { CmcmcFunctionList[[nimbleMCMCs[1]]] <<- CmcmcFunctionList_temp
        } else                { CmcmcFunctionList                   <<- CmcmcFunctionList_temp }
      })
      addTimeResult('nimble_compile', timeResult)
      
      for(iMCMC in seq_along(nimbleMCMCs)) {
        Cmodel$setInits(inits);     calculate(Cmodel)
        mcmcTag <- nimbleMCMCs[iMCMC]
        Cmcmc <- CmcmcFunctionList[[mcmcTag]]
        if(setSeed) set.seed(0)
        if (length(RmcmcTargetList[[iMCMC]])>0){
          monitorVars2 <- monitorVars1
          for (i in 1 : length(RmcmcTargetList[[iMCMC]])){
            
            monitorVars2 <- c(monitorVars2, paste0(mcmcTag, RmcmcTargetList[[iMCMC]][i]))
            
          }
          monitorNodesNIMBLE <<- monitorVars2
          nMonitorNodes <<- length(monitorVars2)
          
          
        }
        
        Cmcmc$run(niter, time = TRUE)
        timeResults <-Cmcmc$getTimes()
        timeEach  <- rep(timeResult[3], length(monitorVars1))
        timeOthers <- 0
        CmvSamples <- Cmcmc$mvSamples
        samplesArray <- as.matrix(CmvSamples)
        colnames(samplesArray) = monitorVars1
        
        for(i in 1: (length(timeResults)-length(RmcmcTargetList[[iMCMC]]))){
          timeOthers = timeOthers + timeResults[i] 
        }
	  	
	

        if (length(RmcmcTargetList[[iMCMC]])>0){
          for (i in 1 : length(RmcmcTargetList[[iMCMC]])){
            
            beforeSamples <- Cmcmc$samplerFunctions$contentsList[[RmcmcTargetList[[iMCMC]][i]]]$before
            afterSamples <- Cmcmc$samplerFunctions$contentsList[[RmcmcTargetList[[iMCMC]][i]]]$after
            x = cbind(beforeSamples, afterSamples)
            samplesArray <- cbind(samplesArray, codess(x=x, tuning=tuning))
            timeEach <- c(timeEach, timeOthers+timeResults[[RmcmcTargetList[[iMCMC]][i]]])
            
            
          }
          
          
        }
        samplesArray <- samplesArray[(burnin+1):floor(niter/thin), , drop=FALSE]
        addToOutput(mcmcTag, samplesArray, timeEach, monitorNodesNIMBLE,iMCMC)
      }
      
    },
    
    addToOutput = function(MCMCtag, samplesArray, timeEach, monitorNodesNIMBLE,iMCMC) {
      output$samples[[MCMCtag]] <<- t(samplesArray) ## makes dim1:monitors, and dim2:iter
      rownames(output$samples[[MCMCtag]]) <<- RmcmcNamesList[[iMCMC]]
      addTimeResult(MCMCtag, timeEach)
      summaryArray <- array(NA, c(nSummaryStats, nMonitorNodes))
      dimnames(summaryArray) <- list(summaryStatDimNames, RmcmcNamesList[[iMCMC]])
      
      for(iStat in seq_along(summaryStats)) {
        summaryArray[iStat, ] <- apply(samplesArray, 2, summaryStatFunctions[[iStat]])
        
      }
      if(calculateEfficiency) {
        for (i in 1 : nMonitorNodes){
          essDim <- which(summaryStatDimNames == 'ess')
          effDim <- which(summaryStatDimNames == 'efficiency')
          codamethodDim <- which(summaryStatDimNames == 'codamethod')
          if(i> length(monitorVars1))
            summaryArray[codamethodDim,i ] <- 0  
          thisTime <- timeEach[i]
          summaryArray[effDim,i ] <- summaryArray[essDim,i] / thisTime
          
          
          
        }
        
      }
      output$summary[[MCMCtag]] <<- t(summaryArray)
      output$monitor[[MCMCtag]] <<- monitorNodesNIMBLE
      if(calculateEfficiency) {
        output$efficiency[[MCMCtag]]$min  <<- min(summaryArray[effDim,] )
        output$efficiency[[MCMCtag]]$mean <<- mean(summaryArray[effDim,] )
      }
    },
    
    addTimeResult = function(MCMCtag, timeEach) {
      output$timing[MCMCtag] <<- timeEach[1]
    },
    
    generate_plots = function() {
      cols <- c(2:6, 8:9)
      #if(nMCMCs > length(cols))    { message('too many MCMCs to plot'); return() }
      
      ## for each monitorNode, generate traceplot for each MCMC
      #for(monitorNode in monitorNodesNIMBLE) {
      for(iMCMC in seq_along(nimbleMCMCs)){
        
        # dev.new()
        nSamplers <- nrow(output$samples[[nimbleMCMCs[iMCMC]]])
        npage= 0
        if (nSamplers>3)
          npage = floor(nSamplers/3)
        nlast = nSamplers - npage*3
      if(npage>0){
        for (j in 0:(npage-1)){
        par(mfrow=c(3,1), mar=c(3,3,2,1), mgp=c(0,0.6,0), tcl=-0.3)
        for(i in 1:3) {
          plot(x=1:nkeep, y=output$samples[[nimbleMCMCs[iMCMC]]][j*3+i, ],
               main=paste0(names(output$summary)[iMCMC], ' traceplot:  ', RmcmcNamesList[[iMCMC]][j*3+i]),
               type='l', col=cols[j*3+i], , xlab='', ylab='', xaxt='n', bty='l') }
        
        }
       }
        
        
        
        #if(savePlot)   { dev.print(device = pdf, file = filename) }
      
      if(nlast>0){
	par(mfrow=c(3,1), mar=c(3,3,2,1), mgp=c(0,0.6,0), tcl=-0.3)
        for(i in 1:nlast) {
          plot(x=1:nkeep, y=output$samples[[nimbleMCMCs[iMCMC]]][npage*3+i, ],
               main=paste0(names(output$summary)[iMCMC], ' traceplot:  ', RmcmcNamesList[[iMCMC]][npage*3+i]),
               type='l', col=cols[npage*3+i], , xlab='', ylab='', xaxt='n', bty='l') }
        
      }
      filename <- paste0(names(output$summary)[iMCMC], '_traceplots_','.pdf')
      ## density plots
      #dev.new()
      }
      npage1= 0
      if (length(output$samples)>3)
        npage1 = floor(length(output$samples)/3)
      nlast1 = length(output$samples)- npage1*3
      if(npage1>0){
        for (j in 0:(npage1-1)){
          par(mfrow = c(3,1), mar=c(3,3,2,1), mgp=c(0,0.6,0), tcl=-0.3)
          for(iMCMC in 1:3){
            nSamplers <- nrow(output$samples[[j*3+iMCMC]])
            densityList <- apply(output$samples[[j*3+iMCMC]][ ,drop=FALSE], 1, density)
            xlim <- range(unlist(lapply(densityList, function(d) d$x)))
            xlim <- mean(xlim) + (xlim-mean(xlim)) * 1.1
            ymax <- max(unlist(lapply(densityList, function(d) d$y))) * 1.1
            plot(-100, -100, xlim=xlim, ylim=c(0,ymax),
                 main=paste0('posterior density:  ', names(output$summary)[j*3+iMCMC]),
                 xlab='', ylab='', yaxt='n', bty='n')
            legend(x='topleft', legend=rownames(output$samples[[j*3+iMCMC]]), lty=1, lwd=2, col=cols[1:nSamplers], bty='n')
            for(i in 1:nSamplers)     polygon(densityList[[i]], border=cols[i])
            abline(h=0, col='white')
          }
          
          #  if(savePlot)   { dev.print(device = pdf, file = filename) }
        }
      }
      if(nlast1>0){
        par(mfrow = c(3,1), mar=c(3,3,2,1), mgp=c(0,0.6,0), tcl=-0.3)
        for(iMCMC in 1:nlast1){
          nSamplers <- nrow(output$samples[[npage1*3+iMCMC]])
          densityList <- apply(output$samples[[npage1*3+iMCMC]][ ,drop=FALSE], 1, density)
          xlim <- range(unlist(lapply(densityList, function(d) d$x)))
          xlim <- mean(xlim) + (xlim-mean(xlim)) * 1.1
          ymax <- max(unlist(lapply(densityList, function(d) d$y))) * 1.1
          plot(-100, -100, xlim=xlim, ylim=c(0,ymax),
               main=paste0('posterior density:  ', names(output$summary)[npage1*3+iMCMC]),
               xlab='', ylab='', yaxt='n', bty='n')
          legend(x='topleft', legend=rownames(output$samples[[npage1*3+iMCMC]]), lty=1, lwd=2, col=cols[1:nSamplers], bty='n')
          for(i in 1:nSamplers)     polygon(densityList[[i]], border=cols[i])
          abline(h=0, col='white')
        }
        
        #  if(savePlot)   { dev.print(device = pdf, file = filename) }
      }
      
    },
    
    checkMCMCdefNames = function() {
      if(!all(nimbleMCMCs %in% MCMCdefNames)) stop(paste0('missing MCMCdefs for: ', paste0(setdiff(nimbleMCMCs, MCMCdefNames), collapse=', ')))
    },
    
    writeModelFile = function() {
      writeLines(paste0('model\n', paste0(deparse(code), collapse='\n')), con=modelFileName)
    },
    
    fileToList = function(file) {
      if(!file.exists(file)) {
        warning(paste0('missing Stan input file: \'', file, '\''))
        return(NULL)
      }
      env <- new.env()
      source(file, local = env)
      lst <- list()
      for(name in ls(env))   lst[[name]] <- get(name, env)
      return(lst)
    },
    
    show = function() {
      cat(paste0('MCMC_CODESS object\n',
                 'algorithms:  ', paste0(MCMCs, collapse=', '), '\n',
                 'monitors:  ', paste0(monitorNodesNIMBLE, collapse=', '), '\n',
                 'model code:\n',
                 paste0(deparse(model), collapse='\n')))
    }
  )
)




