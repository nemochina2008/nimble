#' Auto adapt procedure using Auto Block for efficient MCMC sampling
#' 
#' This is an Auto Adapt version of Auto Block algorithm. Runs NIMBLE's automated blocking procedure for a given model object, to dynamically determine a blocking scheme of the continuous-valued model nodes.  This blocking scheme is designed to produce efficient MCMC sampling (defined as number of effective samples generated per second of algorithm runtime).  However, instead of using Adaptive RW scalar and block sampler, this also (optionally) uses all possible candidate samplers.
#' 
#' 
#' @author Dao Nguyen
#'
#' @param Rmodel A NIMBLE model object, created from \link{nimbleModel}.
#'
#' @param autoIt The number of MCMC iterations to run intermediate MCMC algorithms, through the course of the procedure.  Default 20,000.
#'
#' @param run List of additional MCMC algorithms to compare against the automated blocking MCMC.  These may be specified as: the character string 'all' to denote blocking all continuous-valued nodes; the character string 'default' to denote NIMBLE's default MCMC configuration; a named list element consisting of a quoted code block, which when executed returns an MCMC configuration object for comparison; a custom-specificed blocking scheme, specified as a named list element which itself is a list of character vectors, where each character vector specifies the nodes in a particular block.  Default is c('all', 'default').
#'
#' @param verbose Logical specifying whether to output considerable details of the automated block procedure, through the course of execution.  Default FALSE.
#' 
#' @param setSeed Logical specificying whether to call set.seed(0) prior to beginning the blocking procedure.  Default TRUE.
#'
#' @param makePlots Logical specifying whether to plot the hierarchical clustering dendrograms, through the course of execution.  Default FALSE.
#'
#' @param round Logical specifying whether to round the final output results to two decimal places.  Default TRUE.
#' 
#' @return Returns a named list containing elements:
#' \itemize{
#' \item \code{summary}: A data frame containing a numerical summary of the performance of all MCMC algorithms (including that from automated blocking)
#' \item \code{autoGroups}: A list specifying the parameter blockings converged on by the automated blocking procedure
#' \item \code{conf}: A NIMBLE MCMC configuration object corresponding to the results of the automated blocking procedure
#' }
#' 
#' @references
#'
#'  
#'
#' @export
autoAdapt <- function(Rmodel,
                      autoIt = 20000,
                      run = list('all', 'default'),
                      setSeed = TRUE,
                      verbose = FALSE,
                      makePlots = FALSE,
                      round = TRUE ) {

    control <- list(niter=niter, setSeed=setSeed, verbose=verbose, makePlots=makePlots)
    ab <- autoAdaptClass(Rmodel, control)
    
    ab$run(DefaultSamplerList, Candidates, monitor, iterations)
#     abList <- list(ab)
#     names(abList)[1] <- 'model'
#     df <- createDFfromABlist(abList, autoIt)
#     dfmin <- reduceDF(df, round = round)
#     cat('\nAuto-Adapting summary:\n')
#     print(dfmin)
#     lastAutoInd <- max(grep('^auto', ab$naming))   ## index of final 'auto' iteration
#     lastAutoGrouping <- ab$grouping[[lastAutoInd]]  ## grouping of final 'auto' iteration
#     nonTrivialGroups <- lastAutoGrouping[unlist(lapply(lastAutoGrouping, function(x) length(x)>1))]
#     if(length(nonTrivialGroups) > 0) {
#         cat('\nAuto-Adapting converged on the node groupings:\n')
#         for(i in seq_along(nonTrivialGroups)) {
#             group <- nonTrivialGroups[[i]]
#             cat(paste0('[', i, '] '))
#             cat(paste0(group, collapse = ', '))
#             cat('\n')
#         }
#     } else cat('\nAuto-Adapting converged on all scalar (univariate) sampling\n')
#     cat('\n')
# ## create a new MCMC conf with the autoAdapt groupings:
#     conf <- configureMCMC(Rmodel, nodes = NULL)
#     for(nodeGroup in lastAutoGrouping) addSamplerToConf(Rmodel, conf, nodeGroup)
#     retList <- list(summary=dfmin, autoGroups=nonTrivialGroups, conf=conf)
	retList <- list(ess=ab$ess, efficiency=ab$essPT)
  	retList$samples <- ab$samples
  	return(invisible(retList))
#     return(invisible(retList))
}



autoAdaptModel <- setRefClass(
    Class = 'autoAdaptModel',
    fields = list(
        Rmodel_orig = 'ANY',
        Rmodel = 'ANY',
        Cmodel = 'ANY',
        md = 'ANY',
        scalarNodeVector = 'character',
        scalarNodeVectorCont = 'character',
        scalarNodeVectorDisc = 'character',
        nodeGroupScalars = 'list',
        nodeGroupAllAdapted = 'list',
        monitorsVector = 'character',
        initialMCMCconf = 'ANY'
    ),
    methods = list(
        initialize = function(Rmodel_orig) {
            Rmodel_orig <<- Rmodel_orig
            md <<- Rmodel_orig$modelDef
            Rmodel <<- Rmodel_orig$newModel(replicate = TRUE, check = FALSE)
            ##nimCopy(from = Rmodel_orig, to = Rmodel, logProb = TRUE)
            ##for(var in ls(Rmodel_orig$isDataEnv)) Rmodel$isDataEnv[[var]] <<- Rmodel_orig$isDataEnv[[var]]  ## copies data flags to the new model
            scalarNodeVector <<- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, returnScalarComponents=TRUE)
			latentNodeVector <<- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, latent=TRUE,returnScalarComponents=TRUE)
            discreteInd <- sapply(scalarNodeVector, function(n) Rmodel$isDiscrete(n), USE.NAMES=FALSE)
            scalarNodeVectorCont <<- scalarNodeVector[!discreteInd]   ## making work with discrete nodes
            scalarNodeVectorDisc <<- scalarNodeVector[ discreteInd]   ## making work with discrete nodes
            if(length(scalarNodeVectorCont) == 0) stop('autoAdapting only works with one or more continuous-valued model nodes')   ## making work with discrete nodes
            nodeGroupScalars <<- c(unique(lapply(scalarNodeVectorDisc, Rmodel$expandNodeNames)), scalarNodeVectorCont)   ## making work with discrete nodes, and also with dmulti distributions
            ##nodeGroupAllAdapted <<- list(scalarNodeVector)   ## making work with discrete nodes
            ##nodeGroupAllAdapted <<- c(lapply(scalarNodeVectorDisc, function(x) x), list(scalarNodeVectorCont))   ## making work with discrete nodes
            nodeGroupAllAdapted <<- c(unique(lapply(scalarNodeVectorDisc, Rmodel$expandNodeNames)), list(scalarNodeVectorCont))   ## making work with discrete nodes, and also with dmulti distributions
            monitorsVector <<- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE)
        },
        ## here is where the initial MCMC conf is created, for re-use -- for new version
        createInitialMCMCconf = function() {
            initialMCMCconf <<- configureMCMC(Rmodel)
            nInitialSamplers <- length(initialMCMCconf$samplerConfs)
            initialMCMCconf$addSampler(type = 'sampler_RW',       target = latentNodeVector[1], print=FALSE)  ## add one RW sampler
      initialMCMCconf$addSampler(type = 'sampler_RW_block', target = latentNodeVector[1], print=FALSE)  ## add one RW_block sampler
      #initialMCMCconf$addSampler(type = 'sampler_conjugate', target = latentNodeVector[1], print=FALSE)  ## add one RW_block sampler
      initialMCMCconf$addSampler(type = 'sampler_slice', target = latentNodeVector[1], print=FALSE)  ## add one RW_block sampler
      initialMCMCconf$addSampler(type = 'sampler_AF_slice', target = latentNodeVector[1:2], control = list(sliceWidths = c(.5,.5), sliceFactorBurnInIters = 500,sliceFactorAdaptInterval = 50, sliceSliceAdaptIters = 50), print=FALSE)  ## add one RW_block sampler     
      initialMCMCconf$addSampler(type = 'RW_rotated_block',target = latentNodeVector[1:2], control = list(factorAdaptInterval = 10, ## eigenVectors will be recalculated every 1,000 iterations
                                                                                                          factorBurnInIters = 100, ## eigenVectors will stop being recalculated after 50,000 iterations
                                                                                                          coordinateProportion = .8,  ## the highest 80% of eigenVector directions will be sampled
                                                                                                          scaleVector = rep(.1, 2), ## initial scale values for the univariate random walk sampler
                                                                                                          scaleAdaptInterval = 20))  ## interval on which 
      
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
        createGroups = function(listOfAdapts = list()) {
            listOfAdapts <- lapply(listOfAdapts, function(blk) Rmodel$expandNodeNames(blk, returnScalarComponents=TRUE))
            if(any(unlist(listOfAdapts) %in% scalarNodeVectorDisc)) stop('cannot put Adapt sampler on discrete-valued model nodes')
            nodes <- scalarNodeVector
            nodes <- setdiff(nodes, unlist(listOfAdapts))
            nodeList <- lapply(nodes, function(x) x)
            for(ng in listOfAdapts) nodeList[[length(nodeList)+1]] <- ng
            return(nodeList)
        },
        resetCmodelInitialValues = function() {
            nimCopy(from = Rmodel_orig, to = Cmodel, logProb = TRUE)
            calculate(Cmodel)
        }
    )
)



autoAdaptParamDefaults <- function() {
	list(
            makePlots = FALSE,
            niter = 1000,
            setSeed = TRUE,
            verbose = FALSE
        )
}


autoAdaptClass <- setRefClass(

    Class = 'autoAdaptClass',
    
    fields = list(

           ## special
    abModel = 'ANY',
    it = 'numeric',
    monitor = 'character',
    ## overall control
    cutree_heights = 'numeric',
    allIndex = 'numeric',
    makePlots = 'logical',
    niter = 'numeric',
    saveSamples = 'logical',
    setSeed0 = 'logical',
    verbose = 'logical',
    ## persistant lists of historical data
    NodeInfo = 'list',
	currentSamplerIndex = 'list',
	ScaleList = 'list',
	AcceptanceRateList = 'list',
	StepSizeList = 'list',
	BestSamplers = 'list',
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
    oldIndex='numeric',
    burnedSamples = 'list',
    empCov = 'list',
    empCor = 'list',
    distMatrix = 'list',
    LeastIndex = 'list',
    hTree = 'list',
    keepTrackTemp = 'matrix',
	nCandidates = 'numeric'

    ),
    
    methods = list(
        
	initialize = function(Rmodel, control=list(),DefaultSamplerList=list(), CandidateSamplerList=list(), monitor= character()) {
            abModel <<- autoAdaptModel(Rmodel)
            defaultsList <- autoAdaptParamDefaults()
            for(i in seq_along(defaultsList)) if(is.null(control[[names(defaultsList)[i]]])) control[[names(defaultsList)[i]]] <- defaultsList[[i]]
            for(i in seq_along(control)) eval(substitute(verbose <<- VALUE, list(verbose=as.name(names(control)[i]), VALUE=control[[i]])))
            it <<- 0
			allIndex <<- length(CandidateSamplerList)
		  	keepTrackTemp <<-matrix(0, ncol=length(CandidateSamplerList)+1, nrow= length(DefaultSamplerList))
		  	nCandidates <<- length(CandidateSamplerList) 
		  	colnames(keepTrackTemp) <<-c(names(CandidateSamplerList),'Include')
		  	rownames(keepTrackTemp) <<-names(DefaultSamplerList)
        },

        run = function(DefaultSamplerList, Candidates, monitor, iterations) {
            
          abModel$createInitialMCMCconf()  ## here is where the initial MCMC conf is created, for re-use -- for new version
          
          oldConf = abModel$initialMCMCconf
          n = length(DefaultSamplerList)
          Indices <- 1:allIndex 
          NodeInfo <<- vector(mode="list", length=n)
          monitor <<-monitor
          methods <- list("conjugate","sampler_RW","RW_log", "sampler_slice","sampler_RW_block","sampler_AF_slice","RW_rotated_block")
          currentSamplerIndex <<-list()
          
          for(i1 in 1: n){
            names(NodeInfo)[i1] <<- names(DefaultSamplerList)[i1]
            NodeInfo[[i1]]$type <<- DefaultSamplerList[[i1]]$type
            NodeInfo[[i1]]$target <<- DefaultSamplerList[[i1]]$target
            NodeInfo[[i1]]$currentIndices <<- c()
            NodeInfo[[i1]]$scale <<- rep(1,nCandidates)  
            NodeInfo[[i1]]$acceptanceRate <<- rep(1,nCandidates)
            NodeInfo[[i1]]$nTry <<- rep(1,nCandidates)
            
            NodeInfo[[i1]]$Nodes <<- unique(c(NodeInfo[[i1]]$Nodes, NodeInfo[[i1]]$target))    
            #random walk log scale
            if ((DefaultSamplerList[[i1]]$type == 'sampler_RW' | DefaultSamplerList[[i1]]$type == 'RW') & !is.null(DefaultSamplerList[[i1]]$control$log)){
              if(DefaultSamplerList[[i1]]$control$log){
                j <- which(Candidates==3)	
                NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
                NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "RW_log"))
                currentSamplerIndex[[i1]]<<-j
                
              } else {
                j <- which(Candidates==2)	        
                NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
                
                NodeInfo[[i1]]$scale[j] <<- 1
                
                
                NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "RW"))
                currentSamplerIndex[[i1]]<<-j
                
              }
            } else if(regexpr('conjugate', NodeInfo[[i1]]$type)>0){
              j <- which(Candidates==1)	
              NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
              NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "conjugate"))
              currentSamplerIndex[[i1]]<<-j
              
            } else if(NodeInfo[[i1]]$type=='sampler_slice'){
              j <- which(Candidates==4)	
              NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
              NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "slice"))
              currentSamplerIndex[[i1]]<<-j
              
            } else if(NodeInfo[[i1]]$type=='sampler_RW_block'){
              j <- which(Candidates==5)	
              NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
              NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "RW_block"))
              currentSamplerIndex[[i1]]<<-j
              
            }  else if(NodeInfo[[i1]]$type=='sampler_AF_slice'){
              j <- which(Candidates==6)	
              NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
              NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "sampler_AF_slice"))
              
            }  else if(NodeInfo[[i1]]$type=='RW_rotated_block'){
              j <- which(Candidates==7)	
              NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
              NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "RW_rotated_block"))
              currentSamplerIndex[[i1]]<<-j
              
            } else {
              j <- which(Candidates==2)	  
              NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
              NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "RW"))
              currentSamplerIndex[[i1]]<<-j
              
            }
            
          } 
          
          BestSamplers <<- DefaultSamplerList
          
          
          bestEfficiency = c(a=0)
          count =0
          bestIndex=0
          index <-2
          
          
          
          
          for(i in 1 : iterations){
            
            confList <- list(createConfFromGroups(DefaultSamplerList))
            print("Config List:")
            confList[[1]]$printSamplers()
            runConfListAndSaveBest(confList, paste0('auto',i), auto=FALSE)
            leastMixing <- names(LeastIndex[[it]])
            print(leastMixing)
            print(LeastIndex[[it]])
            print("Least ESS:")
            print(essPT[[it]][LeastIndex[[it]]])
            if(bestEfficiency < essPT[[it]][LeastIndex[[it]]]){
              bestEfficiency <- essPT[[it]][LeastIndex[[it]]]
              bestEssPT <<- list(essPT[[it]])
              oldIndex <<- LeastIndex[[it]]
              BestSamplers <<- DefaultSamplerList
              
            } else {
              DefaultSamplerList <- BestSamplers
              #LeastIndex[[it]]<- oldIndex
            }
            
            ## NodeInfo store information of the 
            avail <- NodeInfo[[LeastIndex[[it]]]]$currentIndices
            print("Avail:")
            print(avail)
            if(length(Indices[-avail]>0)){
              index <- sample(Indices[-avail])[1]
              print(index)
            } else {
              NodeInfo[LeastIndex[[it]]]$currentIndices <<- c()
              index <- sample(Indices)[1]
              print(index)
            }
            print("index:")
            #DefaultSamplerList[[LeastIndex[[it]]]]$type <- CandidateSamplerList[[index]]$type
            
            
            print(DefaultSamplerList[[LeastIndex[[it]]]]$type)
            NodeInfo[[LeastIndex[[it]]]]$currentIndices <<- unique(c(NodeInfo[[LeastIndex[[it]]]]$currentIndices, index))
            print("current indices:")
            print(NodeInfo[[LeastIndex[[it]]]]$currentIndices)
            
            
            
            if(DefaultSamplerList[[LeastIndex[[it]]]]$type %in% c('sampler_RW_block', 'sampler_AF_slice','RW_rotated_block')  & length(DefaultSamplerList[[LeastIndex[[it]]]]$target)<2){
              
              h = sample(seq(0.5,0.9, by=0.1))[1]
              GroupLM<-try(GroupOfLeastMixing(samples[[it]], leastMixing, h= h))
              if (class(GroupLM) == 'try-error') {
                GroupLM <- leastMixing
              } else if(length(GroupLM)>1){
                print(GroupLM)
                DefaultSamplerList[[LeastIndex[[it]]]]$target <- GroupLM
                
              } else { ## Just block anything to test
                print(c(names(essPT[[it]][LeastIndex[[it]]]),names(essPT[[it]][-LeastIndex[[it]]][1])))
                DefaultSamplerList[[LeastIndex[[it]]]]$target <- c(names(essPT[[it]][LeastIndex[[it]]]),names(essPT[[it]][-LeastIndex[[it]]][1])) 
              }  
              
              
              
            }
            if(length(DefaultSamplerList[[LeastIndex[[it]]]]$target)>1){
              if(DefaultSamplerList[[LeastIndex[[it]]]]$type=='sampler_RW_block'){
                
              } else if (DefaultSamplerList[[LeastIndex[[it]]]]$type=='sampler_AF_slice'){
                DefaultSamplerList[[LeastIndex[[it]]]]$control=list(sliceWidths = rep(.5, length(DefaultSamplerList[[LeastIndex[[it]]]]$target) ), sliceFactorBurnInIters = 5000,sliceFactorAdaptInterval = 250, sliceSliceAdaptIters = 20, sliceMaxSteps=20)
                
              } else if (DefaultSamplerList[[LeastIndex[[it]]]]$type=='RW_rotated_block'){
                DefaultSamplerList[[LeastIndex[[it]]]]$control=list(factorAdaptInterval =1000, ## eigenVectors will be recalculated every 1,000 iterations
                                                                    factorBurnInIters = 50000, ## eigenVectors will stop being recalculated after 50,000 iterations
                                                                    coordinateProportion = .6,  ## the highest 80% of eigenVector directions will be sampled
                                                                    scaleVector = rep(.3, length(DefaultSamplerList[[LeastIndex[[it]]]]$target)), ## initial scale values for the univariate random walk sampler
                                                                    scaleAdaptInterval = 250)
              } else {
                
                DefaultSamplerList[[LeastIndex[[it]]]]$type <- 'RW_rotated_block'
              }
            } 
            
            for(i1 in 1: n){
              
              #random walk log scale
              if ((DefaultSamplerList[[i1]]$type == 'sampler_RW' | DefaultSamplerList[[i1]]$type == 'RW') & !is.null(DefaultSamplerList[[i1]]$control$log)){
                if(DefaultSamplerList[[i1]]$control$log){
                  j <- which(Candidates==3)	
                  NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
                  NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "RW_log"))
                  currentSamplerIndex[[i1]]<<-j
                  NodeInfo[[i1]]$scale[j] <<- 1
                } else {
                  j <- which(Candidates==2)	        
                  NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
                  
                  NodeInfo[[i1]]$scale[j] <<- 1
                  
                  
                  NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "RW"))
                  currentSamplerIndex[[i1]]<<-j
                  
                }
              } else if(regexpr('conjugate', NodeInfo[[i1]]$type)>0){
                j <- which(Candidates==1)	
                NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
                NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "conjugate"))
                currentSamplerIndex[[i1]]<<-j
                
              } else if(NodeInfo[[i1]]$type=='sampler_slice'){
                j <- which(Candidates==4)	
                NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
                NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "slice"))
                currentSamplerIndex[[i1]]<<-j
                
              } else if(NodeInfo[[i1]]$type=='sampler_RW_block'){
                j <- which(Candidates==5)	
                NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
                NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "RW_block"))
                currentSamplerIndex[[i1]]<<-j
                
              }  else if(NodeInfo[[i1]]$type=='sampler_AF_slice'){
                j <- which(Candidates==6)	
                NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
                NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "sampler_AF_slice"))
                
              }  else if(NodeInfo[[i1]]$type=='RW_rotated_block'){
                j <- which(Candidates==7)	
                NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
                NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "RW_rotated_block"))
                currentSamplerIndex[[i1]]<<-j
                
              } else {
                j <- which(Candidates==2)	  
                NodeInfo[[i1]]$currentIndices <<- unique(c(NodeInfo[[i1]]$currentIndices, j))
                NodeInfo[[i1]]$Samplers <<- unique(c(NodeInfo[[i1]]$Samplers, "RW"))
                currentSamplerIndex[[i1]]<<-j
                
              }
              
            } 
            print("Now, the sampler of the least mixing is:")
            print(DefaultSamplerList[[LeastIndex[[it]]]]$type)
            print(DefaultSamplerList[[LeastIndex[[it]]]]$target)
            print("Efficiency is:")
            print(bestEfficiency)
            print("Current sampler indices:")
            print(currentSamplerIndex)
            
            print("NodeInfo:")
            print(NodeInfo)
            
            print("Scale List:")
            print(ScaleList)
            
            
            
            
            
          }  
          
          
          
          
          
          #names(candidateGroups) <<- naming
          #names(grouping) <<- naming
          #names(groupSizes) <<- naming
          #names(groupIDs) <<- naming
          #names(samplers) <<- naming
          
          
        },

        determineCandidateGroupsFromCurrentSample = function() {
            cutree_heights <- c(0.5)
            cutreeList <- lapply(cutree_heights, function(height) cutree(hTree[[it]], h = height))
            names(cutreeList) <- paste0('cut', cutree_heights)
            uniqueCutreeList <- unique(cutreeList)
            for(i in seq_along(uniqueCutreeList)) { for(j in seq_along(cutreeList)) { if(all(uniqueCutreeList[[i]]==cutreeList[[j]])) { names(uniqueCutreeList)[i] <- names(cutreeList)[j]; break } } }
            candidateGroupsList <- lapply(uniqueCutreeList, function(ct) determineGroupsFromCutree(ct))
            return(candidateGroupsList)
        },
        
        determineGroupsFromCutree = function(ct) {
            groupsContOnly <- lapply(unique(ct), function(x) names(ct)[ct==x])   ## making work with discrete nodes
            ##groupsAllNodes <- c(lapply(abModel$scalarNodeVectorDisc, function(x) x), groupsContOnly)   ## making work with discrete nodes
            groupsAllNodes <- c(unique(lapply(abModel$scalarNodeVectorDisc, abModel$Rmodel$expandNodeNames)), groupsContOnly)   ## making work with discrete nodes and dmulti distribution
            return(groupsAllNodes)   ## making work with discrete nodes
        },
        
        runConfListAndSaveBest = function(confList, name, auto=FALSE) {
            #lapply(confList, function(conf) checkOverMCMCconf(conf))
            RmcmcList <- lapply(confList, function(conf) buildMCMC(conf))
            CmcmcList <- compileNimble(RmcmcList, project = abModel$Rmodel)
            if(!is.list(CmcmcList)) CmcmcList <- list(CmcmcList)  ## make sure compileNimble() returns a list...
            timingList <- essList <- essPTList <- essPTminList <- list()
            for(i in seq_along(CmcmcList)) {
                if(setSeed) set.seed(0)
                abModel$resetCmodelInitialValues()
                timingList[[i]] <- as.numeric(system.time(CmcmcList[[i]]$run(niter))[3])
                 burnedSamples <- extractAndBurnSamples(CmcmcList[[i]])
                 essList[[i]] <- apply(burnedSamples, 2, effectiveSize)
                 essList[[i]] <- essList[[i]][essList[[i]] > 0]  ## exclude nodes with ESS=0 -- for discrete nodes which are fixed to a certain value; making work with discrete nodes
                 essPTList[[i]] <- essList[[i]] / timingList[[i]]
                 essPTminList[[i]] <- sort(essPTList[[i]])[1]
            }
            bestInd <- as.numeric(which(unlist(essPTminList) == max(unlist(essPTminList))))
             if(length(bestInd) > 1) stop('there should never be an exact tie for the best...')
             if(!is.null(names(confList))) name <- paste0(name, '-', names(confList)[bestInd])
             
             it <<- it + 1
             naming[[it]] <<- name
             candidateGroups[[it]] <<- lapply(confList, function(conf) determineGroupsFromConf(conf))
             grouping[[it]] <<- candidateGroups[[it]][[bestInd]]
             groupSizes[[it]] <<- determineNodeGroupSizesFromGroups(grouping[[it]])
             groupIDs[[it]] <<- determineNodeGroupIDsFromGroups(grouping[[it]])
             samplers[[it]] <<- determineSamplersFromGroupsAndConf(grouping[[it]], confList[[bestInd]])
             timing[[it]] <<- timingList[[bestInd]]
             ess[[it]] <<- essList[[bestInd]]
             essPT[[it]] <<- sort(essPTList[[bestInd]])
             
			 essPT[[it]] <<- essPTList[[bestInd]]
      		 LeastIndex[[it]] <<- which.min(essPT[[it]])
      		 print(LeastIndex[[it]])
      		 if(it>1)
        		samples[[it-1]] <<- NA

             if(auto) {
                 burnedSamples <- extractAndBurnSamples(CmcmcList[[bestInd]])
                 burnedSamples <- burnedSamples[, abModel$scalarNodeVectorCont]   ## making work with discrete nodes
             
                 ##empCov[[it]] <<- cov(burnedSamples)
                 e <- try(cov(burnedSamples))
                 if(inherits(e, 'try-error')) {
                     message('try-error, going into browser'); browser(); 1; 2
                 } else empCov[[it]] <<- e
             
                 ##empCor[[it]] <<- cov2cor(empCov[[it]])
                 e <- try(cov2cor(empCov[[it]]))
                 if(inherits(e, 'try-error')) {
                     message('try-error, going into browser'); browser(); 3; 4
                 } else empCor[[it]] <<- e
             
                 distMatrix[[it]] <<- as.dist(1 - abs(empCor[[it]]))
             
                 ##hTree[[it]] <<- hclust(distMatrix[[it]])
                 e <- try(hclust(distMatrix[[it]]))
                 if(inherits(e, 'try-error')) {
                     message('try-error, going into browser'); browser(); 5; 6
                 } else hTree[[it]] <<- e
             }
             
             if(verbose) printCurrent(name, confList[[bestInd]])
             if(makePlots && auto) makeCurrentPlots(name)
        },

        extractAndBurnSamples = function(Cmcmc) {
            samples <- as.matrix(Cmcmc$mvSamples)
            ## make sure we don't keep samples from deterministic nodes
            namesToKeep <- setdiff(dimnames(samples)[[2]], abModel$Rmodel$getNodeNames(determOnly=TRUE, returnScalarComponents=TRUE))
            burnedSamples <- samples[(floor(niter/2)+1):niter, namesToKeep]
            burnedSamples
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
          conf <- configureMCMC(oldConf = abModel$initialMCMCconf)  ## new version
          conf$setSamplers()  ## new version -- removes all the samplers from initalMCMCconf
          
          n = length(groups)
          for (i in 1:n){
            if(regexpr('conjugate', groups[[i]]$type)>0){
              if(it<=8){
                if(i>0){
                  conf$addSampler(target = groups[[i]]$target,type = 'sampler_RW', control=list())
                } else {
                  conf$addSampler(target = groups[[i]]$target,type = 'sampler_slice')
                }
              } else {
                print("Try to change from conjugate sampler to slice sampler")
                conf$addSampler(target = groups[[i]]$target,type = 'sampler_RW')
                
              }
              
            } else if(length(groups[[i]]$target)>1 ){
              
              if(groups[[i]]$type=='sampler_RW_block'){
                conf$addSampler(target = groups[[i]]$target, type = 'sampler_RW_block',control =  (list(adaptive = TRUE,popCov=1e-12*diag(length(groups[[i]]$target)), adaptInterval = 20)))
              }
              else if (groups[[i]]$type=='sampler_AF_slice'){ 
                conf$addSampler(target = groups[[i]]$target, type = 'sampler_AF_slice', control = list(sliceWidths = rep(.5, length(groups[[i]]$target) ), sliceFactorBurnInIters = 5000,sliceFactorAdaptInterval = 50, sliceSliceAdaptIters = 20, sliceMaxSteps=20))
              } else {
                conf$addSampler(target = groups[[i]]$target, type ='RW_rotated_block', control = list(factorAdaptInterval =1000, ## eigenVectors will be recalculated every 1,000 iterations
                                                                                                      factorBurnInIters = 50000, ## eigenVectors will stop being recalculated after 50,000 iterations
                                                                                                      coordinateProportion = .6,  ## the highest 80% of eigenVector directions will be sampled
                                                                                                      scaleVector = rep(.3, length(groups[[i]]$target)), ## initial scale values for the univariate random walk sampler
                                                                                                      scaleAdaptInterval = 200))  ## interval on which   
              }
            } else{
              conf$addSampler(target = groups[[i]]$target, type = groups[[i]]$type, control = groups[[i]]$control )  # 
            }
          }
          #conf$addMonitors(monitor)
          
          return(conf)
          
        },
        
        sortGroups = function(groups) {
            eachGroupSorted <- lapply(groups, sort)
            groupsAsStrings <- lapply(eachGroupSorted, function(grp) paste0(grp, collapse = '_'))
            sortedInd <- sort(unlist(groupsAsStrings), index.return = TRUE)$ix
            sortedGroups <- eachGroupSorted[sortedInd]
            return(sortedGroups)
        },
        
        checkOverMCMCconf = function(conf) {
            warn <- FALSE
            for(ss in conf$samplerConfs) {
                ## if(ss$name == 'posterior_predictive') {
                ##     msg <- 'using \'posterior_predictive\' sampler may lead to results we don\'t want'
                ##     cat(paste0('\nWARNING: ', msg, '\n\n')); warning(msg)
                ## }
                if(grepl('^conjugate_', ss$name) && getNimbleOption('verifyConjugatePosteriors')) {
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

addSamplerToConf <- function(Rmodel, conf, nodeGroup) {
    if(length(nodeGroup) > 1) {
        conf$addSampler(target = nodeGroup, type = 'RW_block', print = FALSE); return()
    }
    if(!(nodeGroup %in% Rmodel$getNodeNames()) && !Rmodel$isDiscrete(nodeGroup)) {
        conf$addSampler(target = nodeGroup, type = 'RW', print = FALSE); return()
    }
    if(nodeGroup %in% Rmodel$getMaps('nodeNamesEnd')) {
        ##cat(paste0('warning: using \'posterior_predictive\' sampler for node ', nodeGroup, ' may lead to results we don\'t want\n\n'))
        conf$addSampler(target = nodeGroup, type = 'posterior_predictive', print = FALSE); return()
    }
    ## conjugacyResult <- Rmodel$checkConjugacy(nodeGroup)
    ## if((!is.null(conjugacyResult)) && conjOveride) {
    ##     conf$addSampler(target = ??????, type = conjugacyResult$samplerType, control = conjugacyResult$control, print = FALSE); return()
    ## }
    if(Rmodel$isBinary(nodeGroup)) {
        conf$addSampler(target = nodeGroup, type = 'binary', print = FALSE); return()
    }
    if(Rmodel$isDiscrete(nodeGroup)) {
        if(Rmodel$getNodeDistribution(nodeGroup) == 'dmulti') {
            conf$addSampler(target = nodeGroup, type = 'RW_multinomial', print = FALSE); return()
        }
        conf$addSampler(target = nodeGroup, type = 'slice', print = FALSE); return()
    }
    if(length(Rmodel$expandNodeNames(nodeGroup, returnScalarComponents = TRUE)) > 1) {
        conf$addSampler(target = nodeGroup, type = 'RW_block', print = FALSE); return()
    }
    conf$addSampler(target = nodeGroup, type = 'RW', print = FALSE); return()
}

createDFfromABlist <- function(lst, niter) {
    df <- data.frame(model=character(), Adapting=character(), timing=numeric(), node=character(), groupSize = numeric(), groupID = numeric(), sampler = character(), ess=numeric(), essPT=numeric(), stringsAsFactors=FALSE)
    for(iAB in seq_along(lst)) {
        ab <- lst[[iAB]]
        abName <- names(lst)[iAB]
        for(iAdapt in seq_along(ab$naming)) {
            Adapting <- ab$naming[[iAdapt]]
            timing <- ab$timing[[iAdapt]]
            ess <- ab$ess[[iAdapt]]
            nodes <- names(ess)
            essPT <- ab$essPT[[iAdapt]][nodes]            ## sort
            groupSizes <- ab$groupSizes[[iAdapt]][nodes]  ##
            groupIDs <- ab$groupIDs[[iAdapt]][nodes]      ##
            samplers <- ab$samplers[[iAdapt]][nodes]      ##
            newIndDF <- (1:length(nodes)) + dim(df)[1]
            df[newIndDF,] <- NA
            df[newIndDF,]$model <- abName
            df[newIndDF,]$Adapting <- Adapting
            df[newIndDF,]$timing <- timing
            df[newIndDF,]$node <- nodes
            df[newIndDF,]$groupSize <- groupSizes
            df[newIndDF,]$groupID <- groupIDs
            df[newIndDF,]$sampler <- samplers
            df[newIndDF,]$ess <- ess
            df[newIndDF,]$essPT <- essPT
        }
    }
    df$timePer10k <- df$timing * 10000/niter
    df$essPer10k  <- df$ess    * 10000/niter * 2
    df$Efficiency <- df$essPer10k / df$timePer10k
    df$mcmc <- gsub('-.+', '', df$Adapting)
    return(df)
}



plotABS <- function(df, xlimToMin=FALSE, together) {
    models <- unique(df$model)
    nModels <- length(models)
    if(missing(together)) together <- if(nModels <= 5) TRUE else FALSE
    nVertPlots <- if(together) nModels*2 else nModels
    xVarNames <- c('ess', 'essPT')
    parCmd <- quote(par(mfrow=c(nVertPlots,1),mar=c(1,0,1,0),tcl=-.1,mgp=c(3,0,0),cex.axis=.7))
    if(together) { eval(parCmd) }
    for(xVarName in xVarNames) {
        if(!together) { eval(parCmd) }
        maxMinXVar<-0; for(mod in models) {dfMod<-df[df$model==mod,]; blks<-unique(dfMod$Adapting); for(blk in blks) {maxMinXVar<-max(maxMinXVar,min(dfMod[dfMod$Adapting==blk,xVarName]))}}
        maxXVar <- if(xlimToMin) maxMinXVar else max(df[, xVarName])
        xlim <- c(maxXVar*-0.05, maxXVar)
        maxTiming <- max(df[, 'timing'])
        for(mod in models) {
            dfMod <- df[df$model==mod,]
            Adaptings <- unique(dfMod$Adapting)
            nAdaptings <- length(Adaptings)
            bestBlk<-''; bestEssPT<-0; for(blk in Adaptings) { if(min(dfMod[dfMod$Adapting==blk,'essPT'])>bestEssPT) {bestEssPT<-min(dfMod[dfMod$Adapting==blk,'essPT']); bestBlk<-blk} }
            plot(-100,-100,xlim=xlim,ylim=c(0,nAdaptings+1),xlab='',ylab='',main=paste0(xVarName, ' for model ', mod))
            for(iAdapting in 1:nAdaptings) {
                Adapting <- Adaptings[iAdapting]
                dfModAdapt <- dfMod[dfMod$Adapting==Adapting,]
                xVarValues <- dfModAdapt[,xVarName]
                groupSizes <- dfModAdapt[,'groupSize']
                timing <- dfModAdapt[,'timing'][1]   # only first element
                timingOnXaxis <- timing/maxTiming * xlim[2]
                yCoord <- nAdaptings+1-iAdapting
                lines(x=c(0,timingOnXaxis), y=rep(yCoord,2), lty=1, lwd=2, col='lightgrey')
                col <- if(Adapting == bestBlk) 'green' else 'black'
                text(x=xVarValues, y=yCoord, labels=groupSizes, cex=0.7, col=col)
                col <- if(Adapting == bestBlk) 'green' else 'blue'
                text(x=xlim[1], y=yCoord, labels=Adapting, col=col)
                if(timing==maxTiming) text(xlim[2], yCoord+1, paste0('t = ',round(timing,1)))
            }
        }
    }
}


printMinTimeABS <- function(df, round=TRUE, addAutoMax=TRUE, sortOutput=FALSE) {
    namesToRemove <- intersect(c('groupID', 'sampler'), names(df))
    for(name in namesToRemove) { ind <- which(names(df)==name); df <- df[, -ind] }
    models <- unique(df$model)
    cat('\n')
    dfReturn <- data.frame()
    for(mod in models) {
        dfMod <- df[df$model == mod, ]
        Adaptings <- unique(dfMod$Adapting)
        dfOut <- dfMod[numeric(0), ]
        for(blk in Adaptings) {
            dfModBlk <- dfMod[dfMod$Adapting == blk, ]
            ind <- which(dfModBlk$essPT == min(dfModBlk$essPT))[1]
            dfOut[dim(dfOut)[1] + 1, ] <- dfModBlk[ind, ]
        }
        if(sortOutput) dfOut <- dfOut[sort(dfOut$essPT,index.return=TRUE)$ix, ]
        dimnames(dfOut)[[1]] <- 1:(dim(dfOut)[1])
        if(round) {
            dfOut$timing     <- round(dfOut$timing, 2)
            dfOut$timePer10k <- round(dfOut$timePer10k, 2)
            dfOut$ess        <- round(dfOut$ess, 1)
            dfOut$essPer10k  <- round(dfOut$essPer10k, 1)
            dfOut$essPT      <- round(dfOut$essPT, 1)
            dfOut$Efficiency <- round(dfOut$Efficiency, 1)
        }
        if(addAutoMax && ('auto0' %in% Adaptings)) {
            autoAdaptings <- Adaptings[grepl('^auto', Adaptings)]
            dfAuto <- dfOut[dfOut$Adapting %in% autoAdaptings,]
            maxEffInd <- which(dfAuto$Efficiency == max(dfAuto$Efficiency))
            nextInd <- dim(dfOut)[1] + 1
            dfOut[nextInd,] <- dfAuto[maxEffInd,]
            dfOut[nextInd, 'Adapting'] <- dfOut[nextInd, 'mcmc'] <- 'autoMax'
        }
        print(dfOut)
        cat('\n')
        dfReturn <- rbind(dfReturn, dfOut)
    }
    return(invisible(dfReturn))
}


reduceDF <- function(df, addAutoMax=TRUE, sortOutput=TRUE, round=TRUE) {
    df = data.frame(mcmc=df$mcmc, node=df$node, S=df$essPer10k, C=df$timePer10k, Efficiency=df$Efficiency, stringsAsFactors=FALSE)
    dfOut <- df[numeric(), ]
    mcmcs <- unique(df$mcmc)
    for(mcmc in mcmcs) {
        dfBlk <- df[df$mcmc == mcmc, ]
        ind <- which(dfBlk$Efficiency == min(dfBlk$Efficiency))[1]
        dfOut[dim(dfOut)[1]+1, ] <- dfBlk[ind, ]
    }
    dfOut[dfOut$mcmc=='auto0', 'mcmc'] <- 'All Scalar'
    dfOut[dfOut$mcmc=='all', 'mcmc'] <- 'All Block'
    dfOut[dfOut$mcmc=='default', 'mcmc'] <- 'Default'
    if(addAutoMax) {
        autoAdaptings <- dfOut$mcmc[grepl('^auto', dfOut$mcmc)]
        autoLast <- autoAdaptings[length(autoAdaptings)]
        ## replace autoLast with 'autoMax'
        dfOut[dfOut$mcmc==autoLast, 'mcmc'] <- 'Auto-Adapting'
        ## remove any remaining 'auto#' entries
        dfOut <- dfOut[!dfOut$mcmc %in% autoAdaptings,]
    }
    if(sortOutput) dfOut <- dfOut[sort(dfOut$Efficiency,index.return=TRUE)$ix, ]
    dimnames(dfOut)[[1]] <- 1:(dim(dfOut)[1])
    if(round) {
        dfOut$S          <- round(dfOut$S, 2)
        dfOut$C          <- round(dfOut$C, 2)
        dfOut$Efficiency <- round(dfOut$Efficiency, 2)
    }
    return(dfOut)
}



