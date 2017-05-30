


## returns a list of all deterministic dependents up to the first stochastic dependent,
## omitting any nodes in 'omit', or down the path of 'omit' nodes
## works only in terms of vertex IDs, as in the igraph object.
gd_getDependencies_IDs <- function(graph, maps, nodes, omit, downstream) {
  # nonStochNodes <- which(maps$types != 'stoch') 		We should be able to speed things up by looking up by graphID instead of intersecting...

    nodes <- setdiff(nodes, omit)
    newNodes <- if(length(nodes) > 0)    unlist(maps$edgesFrom2To[nodes]) else integer(0)  ## first set of dependencies, including from LHSinferred
    newNodes <- setdiff(newNodes, omit)
    
    boolLHSinferred <- maps$types[nodes] == 'LHSinferred'
    LHSinferredNodes <- nodes[boolLHSinferred]
    nodes <- setdiff(nodes[!boolLHSinferred], omit) ## filter out LHSinferred
    
    if(length(LHSinferredNodes)>0) {
        ## something like x[1], an inferred piece of x[1:10] because x[1] appeared somewhere on its own
        ## include the node it is from (x[1:10]) and well as *non-inferred* dependencies of that node
        fullNodes <- unique(maps$vertexID_2_nodeID[ LHSinferredNodes ]) ## get the x[1:10]
        fullNodes <- setdiff(fullNodes, omit)  ## filter omits
        nodes <- c(nodes, fullNodes)           ## add to nodes

        fullNodesForRecursion <- fullNodes
      ##  fullNodesForRecursion <- if(downstream)   fullNodes   else  fullNodes[maps$types[fullNodes] != 'stoch'] ## find recursion nodes

        fullNodesDeps <- if(length(fullNodesForRecursion) > 0) unlist(maps$edgesFrom2To[fullNodesForRecursion]) else integer(0) ## get dependencies of x[1:10]
        
        fullNodesDepsLHSinferred <- maps$types[fullNodesDeps] == 'LHSinferred' ## filter out LHSinferred dependencies, e.g. x[2]
        fullNodesDeps <-  fullNodesDeps[!fullNodesDepsLHSinferred]
        fullNodesDeps <- setdiff(fullNodesDeps, omit) ## filter omits
        
        newNodes <- c(newNodes, fullNodesDeps)
    }
    
    while(length(newNodes) > 0) {
        nodes <- c(nodes, newNodes)
        newNodesForRecursion <- if(downstream)   newNodes   else  newNodes[maps$types[newNodes] != 'stoch'] 
        newNodes <- if(length(newNodesForRecursion) > 0)  unlist(maps$edgesFrom2To[newNodesForRecursion]) else integer(0) 
        newNodes <- setdiff(newNodes, omit)
    }
    nodes <- unique(nodes)
    nodes <- sort(nodes)    # topological sort
    return(nodes)

    ## old
    ## nodes <- setdiff(nodes, omit)
    ## newNodes <- if(length(nodes) > 0)    unlist(maps$edgesFrom2To[nodes]) else integer(0) 
    ## newNodes <- setdiff(newNodes, omit)
    ## while(length(newNodes) > 0) {
    ##     nodes <- c(nodes, newNodes)
    ##     newNodesForRecursion <- if(downstream)   newNodes   else  newNodes[maps$types[newNodes] != 'stoch'] 
    ##     newNodes <- if(length(newNodesForRecursion) > 0)  unlist(maps$edgesFrom2To[newNodesForRecursion]) else integer(0) 
    ##     newNodes <- setdiff(newNodes, omit)
    ## }
    ## nodes <- unique(nodes)
    ## nodes <- sort(nodes)    # topological sort
    ## return(nodes)

}


gd_allNeighbors <- function(graph, nodes) stop("shouldn't be calling gd_allNeighbors any more")

enhanceDepsForDerivs <- function(inputNodes, deps, model) {
    ## exploratory idea for enhancing an ordered set of nodes with parents, by input order, indexed by which input node they correspond to.
    ## Example:
    ## say deps has [2, 4, 6]
    ## 

    ## path to info: idea 1: go up first
    ## graphID_2_declID
    ## from the declID we 

    ## idea 2: 
    ## from input nodes, see edgesFrom2To and see edgesFrom2ParentExprID
    ## see if any of the To nodes are in deps and provide the parentExprID

    if(!is.integer(inputNodes)) inputNodes <- model$modelDef$nodeName2GraphIDs(inputNodes)
    depIDs <- if(!is.integer(deps)) model$modelDef$nodeName2GraphIDs(deps) else deps

    maps <- model$modelDef$maps
    
    declIDs <- maps$graphID_2_declID[depIDs]
    
    depIndex_2_parentDepIndices <- lapply(declIDs, integer)

    for(i in seq_along(depIDs)) {
        thisNode <- depIDs[i]
        if(thisNode %in% inputNodes) {
            depIndex_2_parentDepIndices[[i]] <- -which(inputNodes == thisNode) ## e.g. set -2 for 2nd input node
        }
        toNodes <- maps$edgesFrom2To[[ thisNode ]]
        parentExprIDs <- maps$edgesFrom2ParentExprID[[ thisNode ]]
        for(iTo in seq_along(toNodes)) {
            thisToNode <- toNodes[iTo]
            if(thisToNode %in% depIDs) {
                iThisNodeInDeps <- which(depIDs == thisToNode)
                thisParentExprID <- parentExprIDs[iTo]
                depIndex_2_parentDepIndices[[iThisNodeInDeps]][ thisParentExprID ] <- i
            }
        }
    }
    list(deps, depIndex_2_parentDepIndices)
}

explainDerivContent <- function(enhancedDeps, model) {
    writeLines('The following calculations would be done from this input:')
    deps <- enhancedDeps[[1]]
    depIDs <- if(!is.integer(deps)) model$modelDef$nodeName2GraphIDs(deps) else deps
    declIDs <- model$modelDef$maps$graphID_2_declID[depIDs]
    derivInfo <- enhancedDeps[[2]]
    for(i in seq_along(depIDs)) {
        thisDerivInfo <- derivInfo[[i]]
        if(length(thisDerivInfo) > 0) {
            if(all(thisDerivInfo < 0)) { ## initiating node
                if(length(thisDerivInfo) > 1) stop('Problem with thisDerivInfo')
                description <- paste0('Input parameter ', -thisDerivInfo)
            } else {
                used <- thisDerivInfo > 0
                if(sum(used) > 0) {
                    argumentNames <- lapply( m$modelDef$declInfo[[ declIDs[i] ]]$symbolicParentNodesReplaced, deparse)
                    description <- paste0('Argument ', seq_along(thisDerivInfo)[used], ' (', argumentNames[used],') comes from calculation ', thisDerivInfo[used])
                }
                else
                    description <- "(No arguments are from previous calculations)\n"
            }
        } else
            description <- ""
        BUGSline <- deparse(m$modelDef$declInfo[[ declIDs[i] ]]$codeReplaced)
        output <- paste0(i,': ', deps[i], ' (from ', BUGSline, ')\n', paste0('\t', description, collapse = '\n'))
        writeLines(output)
    }
}
