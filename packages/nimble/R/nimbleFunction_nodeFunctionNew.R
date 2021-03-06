nodeFunctionNew <- function(LHS, RHS, name = NA, altParams, bounds, parentsSizeAndDims, logProbNodeExpr, type, setupOutputExprs, dynamicIndexInfo = NULL, evaluate = TRUE, where = globalenv()) {
    if(!(type %in% c('stoch', 'determ')))       stop(paste0('invalid argument to nodeFunction(): type = ', type))
    setupOutputLabels <- nndf_makeNodeFunctionIndexLabels(setupOutputExprs) ## should perhaps move to the declInfo for preservation
    LHSrep <- nndf_replaceSetupOutputsWithIndexedNodeInfo(LHS, setupOutputLabels)
    RHSrep <- nndf_replaceSetupOutputsWithIndexedNodeInfo(RHS, setupOutputLabels)
    altParamsRep <- lapply(altParams, nndf_replaceSetupOutputsWithIndexedNodeInfo, setupOutputLabels)
    boundsRep <- lapply(bounds, nndf_replaceSetupOutputsWithIndexedNodeInfo, setupOutputLabels)
    logProbNodeExprRep <- nndf_replaceSetupOutputsWithIndexedNodeInfo(logProbNodeExpr, setupOutputLabels)
    
    if(nimbleOptions()$allowDynamicIndexing) {
        if(length(dynamicIndexInfo)) {
            for(i in seq_along(dynamicIndexInfo))
                dynamicIndexInfo[[i]]$indexCode <- nndf_replaceSetupOutputsWithIndexedNodeInfo(dynamicIndexInfo[[i]]$indexCode, setupOutputLabels)
            dynamicIndexLimitsExpr <- nndf_generateDynamicIndexLimitsExpr(dynamicIndexInfo)
        } else dynamicIndexLimitsExpr <- NULL
    } else dynamicIndexLimitsExpr <- NULL
    
    if(nimbleOptions('experimentalEnableDerivs')){
      parents <- names(parentsSizeAndDims)
      parentIndexInfoList <- nndf_extractNodeIndices(LHSrep, parents)
      parentIndexInfoList <- nndf_extractNodeIndices(RHSrep, parents, indexExprList = parentIndexInfoList)
      for(i in seq_along(parentIndexInfoList)){
        for(j in seq_along(parentIndexInfoList[[i]])){
          parentsSizeAndDims[[names(parentIndexInfoList)[i]]][[j]]$indexExpr <- parentIndexInfoList[[i]][[j]]$indexExpr 
        }
      }
    }

    nodeFunctionTemplate <-
        substitute(
            nimbleFunction(##contains      = CONTAINS,
                           setup         = SETUPFUNCTION,
                           methods       = METHODS,
                           name          = name,
                           check         = FALSE,
                           enableDerivs  = CALCAD_LIST,
                           where = where)
          ,
            list(##CONTAINS      = nndf_createContains(RHS, type), ## this was used for intermediate classes for get_scale style parameter access, prior to getParam
                 SETUPFUNCTION = nndf_createSetupFunction(),  ##nndf = new node function
                 METHODS       = nndf_createMethodList(LHSrep, RHSrep, parentsSizeAndDims, altParamsRep, boundsRep, logProbNodeExprRep, type, dynamicIndexLimitsExpr, RHS),
                 CALCAD_LIST   = if(nimbleOptions('experimentalEnableDerivs')) list(getCalcADFunName()) else list(),
                 where         = where)
        )
    if(evaluate){
      returnFunc <- eval(nodeFunctionTemplate)
      assign('parentsSizeAndDims', parentsSizeAndDims, envir = environment(returnFunc))      
      return(returnFunc)
    }
    else       return(nodeFunctionTemplate)
}

nndf_makeNodeFunctionIndexLabels <- function(setupOutputExprs) {
    setupOutputLabels <- 1:length(setupOutputExprs)
    names(setupOutputLabels) <- names(setupOutputExprs)
    setupOutputLabels
}

nndf_makeNodeFunctionIndexAccessCall <- function(index) {
    ## use name INDEXEDNODEINFO_
    substitute(getNodeFunctionIndexedInfo(INDEXEDNODEINFO_, INDEX), list(INDEX = index)) ## still needs unity decrement for c++
}

nndf_extractNodeIndices <- function(code, nodesToExtract, indexExprList = list()){
  if(is.call(code)){
    if(deparse(code[[1]]) == '[') {
      if(deparse(code[[2]]) %in% nodesToExtract){
        thisIndexExpr <- list()
        for(i in 1:(length(code)-2)){
          if(is.call(code[[i + 2]]) && deparse(code[[i+2]][[1]]) == ':'){
            thisIndexExpr <- c(thisIndexExpr, code[[i+2]][[2]])
          }
          else{
            thisIndexExpr <- c(thisIndexExpr, code[[i+2]])
          }
        }
        if(is.null(indexExprList[[deparse(code[[2]])]])) indexExprList[[deparse(code[[2]])]][[1]]$indexExpr <- thisIndexExpr
        else indexExprList[[deparse(code[[2]])]][[length(indexExprList[[deparse(code[[2]])]]) + 1]]$indexExpr <- thisIndexExpr
        return(indexExprList)
      }
    }
    if(length(code) > 1){
      for(i in 2:length(code)){
        indexExprList <- nndf_extractNodeIndices(code[[i]], nodesToExtract, indexExprList)
      }
    }
  }
  return(indexExprList)
}

nndf_replaceSetupOutputsWithIndexedNodeInfo <- function(code, setupOutputLabels) {
    cLength <- length(code)
    if(cLength == 1) {
        if(is.name(code)) {
            varName <- as.character(code)
            if(varName %in% names(setupOutputLabels))
                return(nndf_makeNodeFunctionIndexAccessCall(setupOutputLabels[[varName]]))
        }
        return(code)
    }
    if(is.call(code)) {
        for(i in 2:cLength)
            code[[i]] <- nndf_replaceSetupOutputsWithIndexedNodeInfo(code[[i]], setupOutputLabels)
        return(code)
    }
    return(code)
}

# need function to be defined to pass CRAN but setupOutputs is
# never called - it is processed out of nimbleFunction setup code
setupOutputs <- function(...) NULL

## creates a function object for use as setup argument to nimbleFunction()
nndf_createSetupFunction <- function() {
    setup <- function(model, BUGSdecl) {
        indexedNodeInfoTable <- indexedNodeInfoTableClass(BUGSdecl)
        setupOutputs(indexedNodeInfoTable)
        invisible(NULL)
    }
    return(setup)
}

indexedNodeInfoTableClass <- function(BUGSdecl) {
    structure(
        list(unrolledIndicesMatrix = BUGSdecl$unrolledIndicesMatrix),
             class = 'indexedNodeInfoTableClass')
}

## creates a list of the methods calculate, simulate, getParam, getBound, and getLogProb, corresponding to LHS, RHS, and type arguments
nndf_createMethodList <- function(LHS, RHS, parentsSizeAndDims, altParams, bounds, logProbNodeExpr, type, dynamicIndexLimitsExpr, RHSnonReplaced) {
    if(type == 'determ') {
        methodList <- eval(substitute(
            list(
                simulate   = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {  DETERMSIM                                                },
                calculate  = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { simulate(INDEXEDNODEINFO_ = INDEXEDNODEINFO_);    returnType(double());   return(invisible(0)) },
                calculateDiff = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {simulate(INDEXEDNODEINFO_ = INDEXEDNODEINFO_);  returnType(double());   return(invisible(0)) },
                getLogProb = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {                returnType(double());   return(0)            }
            ),
            list(LHS=LHS, # no longer used but kept for reference
                 RHS=RHS, # no longer used but kept for reference
                 DETERMSIM = ndf_createDetermSimulate(LHS, RHS, dynamicIndexLimitsExpr = dynamicIndexLimitsExpr, RHSnonReplaced = RHSnonReplaced)
                 )))
        if(nimbleOptions('experimentalEnableDerivs')){
          methodList[['CALCADFUNNAME']]  <- eval(substitute(
            function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { LHS <- RHS;    returnType(double(THISDIM));   return(THISNAME) },
            list(LHS=LHS,
                 RHS=RHS,
                 THISDIM =   as.numeric(parentsSizeAndDims[[1]][[1]]$nDim),
                 THISNAME =  as.name(names(parentsSizeAndDims)[1])
            )))
        }

    }
    if(type == 'stoch') {
        methodList <- eval(substitute(
            list(
                simulate   = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { STOCHSIM                                                         },
                calculate  = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { STOCHCALC_FULLEXPR;   returnType(double());   return(invisible(LOGPROB)) },
                calculateDiff = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {STOCHCALC_FULLEXPR_DIFF; LocalAns <- LocalNewLogProb - LOGPROB;  LOGPROB <<- LocalNewLogProb;
                                            returnType(double());   return(invisible(LocalAns))},
                getLogProb = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {                       returnType(double());   return(LOGPROB)            }
            ),
            list(LHS       = LHS, # no longer used but kept for reference
                 LOGPROB   = logProbNodeExpr,
                 STOCHSIM  = ndf_createStochSimulate(LHS, RHS, dynamicIndexLimitsExpr = dynamicIndexLimitsExpr, RHSnonReplaced = RHSnonReplaced),
                 STOCHCALC_FULLEXPR = ndf_createStochCalculate(logProbNodeExpr, LHS, RHS, dynamicIndexLimitsExpr = dynamicIndexLimitsExpr, RHSnonReplaced = RHSnonReplaced),
                 STOCHCALC_FULLEXPR_DIFF = ndf_createStochCalculate(logProbNodeExpr, LHS, RHS, diff = TRUE, dynamicIndexLimitsExpr = dynamicIndexLimitsExpr, RHSnonReplaced = RHSnonReplaced))))
        if(nimbleOptions('experimentalEnableDerivs')){
          methodList[['CALCADFUNNAME']]  <- eval(substitute(
            function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { STOCHCALC_FULLEXPR_AD;   returnType(double());   return(invisible(LOGPROB_AD)) },
            list(LOGPROB_AD = as.name('logProb'),
                 STOCHCALC_FULLEXPR_AD = ndf_createStochCalculate(as.name('logProb'), LHS, RHS, ADFunc = TRUE, dynamicIndexLimitsExpr = dynamicIndexLimitsExpr, RHSnonReplaced = RHSnonReplaced))))
        }
        if(FALSE) {
        if(nimbleOptions()$compileAltParamFunctions) {
            distName <- as.character(RHS[[1]])
            ## add accessor function for node value; used in multivariate conjugate sampler functions
            type = getType(distName)
            nDim <- getDimension(distName)
            methodList[['get_value']] <- ndf_generateGetParamFunction(LHS, type, nDim)
            ## add accessor functions for stochastic node distribution parameters
            for(param in names(RHS[-1])) {
                if(!param %in% c("lower", "upper")) {
                    type = getType(distName, param)
                    nDim <- getDimension(distName, param)
                    methodList[[paste0('get_',param)]] <- ndf_generateGetParamFunction(RHS[[param]], type, nDim)
                }
            }
            for(i in seq_along(altParams)) {
                altParamName <- names(altParams)[i]
                type = getType(distName, altParamName)
                nDim <- getDimension(distName, altParamName)
                methodList[[paste0('get_',altParamName)]] <- ndf_generateGetParamFunction(altParams[[altParamName]], type, nDim)
            }
        }
        } ## if(FALSE) to cut out old get_XXX param system
            ## new for getParam, eventually to replace get_XXX where XXX is each param name
            ## TO-DO: unfold types and nDims more thoroughly (but all types are implemented as doubles anyway)
            ## understand use of altParams vs. all entries in typesListAllParams (i.e., getDistributionInfo(distName)$types
        ## need a value Entry
        distName <- as.character(RHS[[1]])

        allParams <- c(list(value = LHS), as.list(RHS[-1]), altParams)

        typesNDims <- getDimension(distName, includeParams = TRUE)
        typesTypes <- getType(distName, includeParams = TRUE)
        paramIDs <- getParamID(distName, includeParams = TRUE)
        ## rely on only double for now
        for(nDimSupported in c(0, 1, 2)) {
            boolThisCase <- typesNDims == nDimSupported ## & typesTypes == 'double' ## until (if ever) we have separate handling of integer params, these should be folded in with doubles.  We don't normally have any integer params, because we handle integers as doubles
            paramNamesToUse <- getParamNames(distName)[boolThisCase]
            caseName <- paste0("getParam_",nDimSupported,"D_double")

            ## Special handling needed for dinterval
            ## Second parameter is a vector but is expected to handle a scalar
            ## It goes in nDimSupported == 1, but we need it's return type cast to vector
            ## if it is in fact only a scalar. We do so by wrapping in c() if
            ## it doesn't any `:` in it
            if(nDimSupported == 1) {
                allParams[paramNamesToUse] <- lapply(allParams[paramNamesToUse],
                                                     function(x) {
                                                         if(':' %in% all.names(x))
                                                             x
                                                         else
                                                             substitute(c(X), list(X = x))
                                                     })
            }
            
            if(length(paramNamesToUse) > 0)
                methodList[[caseName]] <- nndf_generateGetParamSwitchFunction(allParams[paramNamesToUse], paramIDs[paramNamesToUse], type = 'double', nDim = nDimSupported)
        }
        nDimSupported <- 0  # even multivar nodes have single lower and single upper since truncation not supported for mv distributions, so lower and upper come from distribution 'range'
        caseName <- paste0("getBound_",nDimSupported,"D_double")
        methodList[[caseName]] <- nndf_generateGetBoundSwitchFunction(bounds, seq_along(bounds), type = 'double', nDim = nDimSupported)
    }
    parentsArgs <-c()
    if(nimbleOptions('experimentalEnableDerivs')){
      names(methodList)[names(methodList) == 'CALCADFUNNAME'] <-  getCalcADFunName() ## replace CALCADFUNNAME with real name
      parentsArgs <- if(length(parentsSizeAndDims) > 0) list() else NULL
      for(i in seq_along(parentsSizeAndDims)){
        for(j in seq_along(parentsSizeAndDims[[i]])){
          parentsArgs[[paste0(names(parentsSizeAndDims)[i], '_', j)]] <- substitute(double(PARDIM, PARSIZES), 
                                                                                    list(PARDIM = as.numeric(parentsSizeAndDims[[i]][[j]]$nDim), 
                                                                                         PARSIZES = nndf_makeParentSizeExpr(parentsSizeAndDims[[i]][[j]])))  
          body(methodList[[getCalcADFunName()]]) <- nndf_addArgInfoToCalcAD(body(methodList[[getCalcADFunName()]]), names(parentsSizeAndDims)[i], j)
        }
      }
      if(type == 'determ') body(methodList[[getCalcADFunName()]]) <- nndf_addArgInfoToCalcAD(body(methodList[[getCalcADFunName()]]), names(parentsSizeAndDims)[1], 1)
      formals(methodList[[getCalcADFunName()]]) <- c(formals(methodList[[getCalcADFunName()]]), parentsArgs)
    }
    ## add model$ in front of all names, except the setupOutputs
  
    methodList <- nndf_addModelDollarSignsToMethods(methodList, exceptionNames = c("LocalAns", "LocalNewLogProb","PARAMID_","PARAMANSWER_", "BOUNDID_", "BOUNDANSWER_", "INDEXEDNODEINFO_"), 
                                                    ADexceptionNames = c(names(parentsArgs), 'logProb'))
    return(methodList)
}

nndf_addArgInfoToCalcAD <- function(code, argName, argNum){
  for(i in seq_along(code)){
    if(is.name(code[[i]]) && (deparse(code[[i]]) == argName)){
      code[[i]] <- parse(text = paste0(deparse(code[[i]]), '_', argNum))[[1]]
      return(code)
    }
    else if(length(code[[i]]) > 1){
      didArgAdd <- nndf_addArgInfoToCalcAD(code[[i]], argName, argNum)
      if(!is.null(didArgAdd)){ 
        code[[i]] <- didArgAdd
        return(code)
      }
    }
  }
  return(NULL)
}

nndf_makeParentSizeExpr <- function(sizeInfoList){
  if(sizeInfoList$nDim == 0) return(parse(text = 'c(1)')[[1]])
  else{
    dimInds <- sizeInfoList$lengths[sizeInfoList$lengths > 1]
    return(parse(text = paste0('c(', paste(dimInds, collapse = ', '), ')'))[[1]])
  }
}

nndf_addModelDollarSignsToMethods <- function(methodList, exceptionNames = character(), ADexceptionNames = character()) {
    for(i in seq_along(methodList)) {
      if(names(methodList)[i] == getCalcADFunName()){
        body(methodList[[i]]) <- removeIndices(body(methodList[[i]]))
        body(methodList[[i]]) <- addModelDollarSign(body(methodList[[i]]), exceptionNames = c(exceptionNames, ADexceptionNames))
      }
      else  body(methodList[[i]]) <- addModelDollarSign(body(methodList[[i]]), exceptionNames = c(exceptionNames))
    }
    return(methodList)
}

nndf_generateGetParamSwitchFunction <- function(typesListAll, paramIDs, type, nDim) {
    if(any(unlist(lapply(typesListAll, is.null)))) stop(paste('problem creating switch function for getParam from ', paste(paste(names(typesListAll), as.character(typesListAll), sep='='), collapse=',')))
    paramIDs <- as.integer(paramIDs)
    answerAssignmentExpressions <- lapply(typesListAll, function(x) substitute(PARAMANSWER_ <- ANSEXPR, list(ANSEXPR = x)))
    switchCode <- as.call(c(list(quote(nimSwitch), quote(PARAMID_), paramIDs), answerAssignmentExpressions))
    # avoid arg name mismatch based on R partial arg name matching
    names(answerAssignmentExpressions) <- NULL
    names(switchCode)[2:3] <- c('paramID', 'IDoptions')
    if(nDim == 0) {
        answerInitCode <- quote(PARAMANSWER_ <- 0)  ## this avoids a Windows compiler warning about a possibly unassigned return variable

        ans <- try(eval(substitute(
            function(PARAMID_ = integer(), INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {
                returnType(TYPE(NDIM))
                ANSWERINITCODE
                SWITCHCODE
                return(PARAMANSWER_)
            },
            list(TYPE = as.name(type), NDIM=nDim, ANSWERINITCODE = answerInitCode, SWITCHCODE = switchCode)
        )))
    } else {
        ans <- try(eval(substitute(
            function(PARAMID_ = integer(), INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {
                returnType(TYPE(NDIM))
                SWITCHCODE
                return(PARAMANSWER_)
            },
            list(TYPE = as.name(type), NDIM=nDim, SWITCHCODE = switchCode)
        )))
    }
    if(inherits(ans, 'try-error')) browser()
    attr(ans, 'srcref') <- NULL
    ans
}

nndf_generateGetBoundSwitchFunction <- function(typesListAll, boundIDs, type, nDim) {
    if(any(unlist(lapply(typesListAll, is.null)))) stop(paste('problem creating switch function for getBound from ', paste(paste(names(typesListAll), as.character(typesListAll), sep='='), collapse=',')))
    boundIDs <- as.integer(boundIDs)
    answerAssignmentExpressions <- lapply(typesListAll, function(x) substitute(BOUNDANSWER_ <- ANSEXPR, list(ANSEXPR = x)))
    switchCode <- as.call(c(list(quote(nimSwitch), quote(BOUNDID_), boundIDs), answerAssignmentExpressions))
    if(nDim == 0) {
        answerInitCode <- quote(BOUNDANSWER_ <- 0)  ## this avoids a Windows compiler warning about a possibly unassigned return variable

        ans <- try(eval(substitute(
            function(BOUNDID_ = integer(), INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {
                returnType(TYPE(NDIM))
                ANSWERINITCODE
                SWITCHCODE
                return(BOUNDANSWER_)
            },
            list(TYPE = as.name(type), NDIM=nDim, ANSWERINITCODE = answerInitCode, SWITCHCODE = switchCode)
        )))
    } else {
        ans <- try(eval(substitute(
            function(BOUNDID_ = integer(), INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {
                returnType(TYPE(NDIM))
                SWITCHCODE
                return(BOUNDANSWER_)
            },
            list(TYPE = as.name(type), NDIM=nDim, SWITCHCODE = switchCode)
        )))
    }
    if(inherits(ans, 'try-error')) browser()
    attr(ans, 'srcref') <- NULL
    ans
}

nndf_generateDynamicIndexLimitsExpr <- function(dynamicIndexInfo) {
    indivLimits <- lapply(dynamicIndexInfo, function(x) substitute(VAREXPR >= LOWER & VAREXPR <= UPPER,
                              list(VAREXPR =  x$indexCode,
                             LOWER = x$lower,
                             UPPER = x$upper)))
    dynamicIndexLimitsExpr <- indivLimits[[1]]
    if(length(indivLimits) > 1)
        for(i in 2:length(indivLimits))
            dynamicIndexLimitsExpr <- substitute(FIRST & SECOND,
                                                 list(FIRST = dynamicIndexLimitsExpr,
                                                      SECOND = indivLimits[[i]]))
    ## FIXME: that puts extra () in the expression; could potentially also construct as
    ## tmp <- quote(1 & 1)
    ## tmp[[2]] <- dynamicIndexLimitsExpr
    ## tmp[[3]] <- indivLimits[[i]]
    ## check if () are stripped out in C++; if so then it doesn't matter anyway
    return(dynamicIndexLimitsExpr)
}

