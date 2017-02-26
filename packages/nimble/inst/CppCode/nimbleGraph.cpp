#include "nimble/nimbleGraph.h"

graphNode::graphNode(int inputCgraphID, NODETYPE inputType, const string &inputName ) :
  role(UNKNOWNROLE),
  type(inputType),
  CgraphID(inputCgraphID),
  name(inputName),
  touched(false),
  touched2(false),
  numChildren(0) {
  RgraphID = CgraphID + 1;
}

void graphNode::addChild( graphNode *toNode, int childParentExpressionID) {
#ifdef _DEBUGNIMGRAPH
  PRINTF("adding child %s to parent %s\n", toNode->name.c_str(), name.c_str());
#endif
  children.push_back(toNode);
  childrenParentExpressionIDs.push_back(childParentExpressionID);
  numChildren++;
  toNode->addParent(this);
}

void graphNode::addParent(graphNode *fromNode) {
#ifdef _DEBUGNIMGRAPH
  PRINTF("adding parent %s to child %s\n", fromNode->name.c_str(), name.c_str());
#endif
  parents.push_back(fromNode);
}

void SEXP_2_nodeType(SEXP Stypes, vector<NODETYPE> &ans) {
  //  enum NODETYPE {UNKNOWNTYPE, STOCH, DETERM, RHSONLY};
  if(!isString(Stypes)) {
    PRINTF("Error:  called for SEXP that is not a string!\n");
    return;
  }
  int nn = LENGTH(Stypes);
  ans.resize(nn);
  string oneString;
  for(int i = 0; i < nn; i++) {
    oneString.assign(CHAR(STRING_ELT(Stypes, i)), LENGTH(STRING_ELT(Stypes, i)));
    if(oneString == "stoch")
      ans[i] = STOCH;
    else if(oneString == "determ")
      ans[i] = DETERM;
    else if(oneString == "RHSonly")
      ans[i] = RHSONLY;
    else if(oneString == "LHSinferred")
      ans[i] = LHSINFERRED;
    else if(oneString == "unknown")
      ans[i] = UNKNOWNTYPE;
    else {
      ans[i] = UNKNOWNTYPE;
      PRINTF("In SEXP_2_nodeType: unknown string type label %s\n", oneString.c_str());
    }
  }
}

SEXP setGraph(SEXP SedgesFrom, SEXP SedgesTo, SEXP SedgesFrom2ParentExprIDs, SEXP Stypes, SEXP Snames, SEXP SnumNodes) {
  vector<int> edgesFrom = SEXP_2_vectorInt(SedgesFrom, -1); // -1 subtracted here
  vector<int> edgesTo = SEXP_2_vectorInt(SedgesTo, -1); // -1 substracted here
  vector<int> edgesFrom2ParentExprIDs = SEXP_2_vectorInt(SedgesFrom2ParentExprIDs);
  vector<NODETYPE> types;
  SEXP_2_nodeType(Stypes, types);
  vector<string> names;
  STRSEXP_2_vectorString(Snames, names);
  int numNodes = SEXP_2_int(SnumNodes);
  nimbleGraph *newGraph = new nimbleGraph;
  newGraph->setNodes(edgesFrom, edgesTo, edgesFrom2ParentExprIDs, types, names, numNodes);
  SEXP SextPtrAns;
  PROTECT(SextPtrAns = R_MakeExternalPtr(newGraph, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(SextPtrAns, &nimbleGraphFinalizer, TRUE);
  UNPROTECT(1);
  return(SextPtrAns);
}

void nimbleGraphFinalizer(SEXP SgraphExtPtr) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr));
  delete graphPtr;
}

nimbleGraph::~nimbleGraph() {
  int n = graphNodeVec.size();
  for(int i = 0; i < n; i++) {
    delete graphNodeVec[i];
  }
}

SEXP anyStochDependencies(SEXP SgraphExtPtr) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr));
  vector<int> ans(graphPtr->anyStochDependencies());
  SEXP Sans;
  PROTECT(Sans = allocVector(LGLSXP, ans.size()));
  int *SansPtr = INTEGER(Sans);
  for(unsigned int i = 0; i < ans.size(); i++) {
    if(ans[i] == 0) PRINTF("Element %i was not processed\n", i);
    SansPtr[i] = ans[i]==2 ? 1 : 0;
  }
  UNPROTECT(1);
  return(Sans);
}

SEXP getDependencies(SEXP SgraphExtPtr, SEXP Snodes, SEXP Somit, SEXP Sdownstream) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr));
  vector<int> nodes = SEXP_2_vectorInt(Snodes, -1); // subtract 1 index for C
  vector<int> omit = SEXP_2_vectorInt(Somit, -1);
  bool downstream = SEXP_2_bool(Sdownstream);
  vector<int> ans = graphPtr->getDependencies(nodes, omit, downstream);
  return(vectorInt_2_SEXP(ans, 1)); // add 1 index for R
}

SEXP topologicalSortOrder(SEXP SgraphExtPtr) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr));
  vector<int> sortOrder = graphPtr->topologicalSortOrder();
  return(vectorInt_2_SEXP(sortOrder, 1));
}

SEXP anyStochParents(SEXP SgraphExtPtr) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr));
  vector<int> ans(graphPtr->anyStochParents());
  SEXP Sans;
  PROTECT(Sans = allocVector(LGLSXP, ans.size()));
  int *SansPtr = INTEGER(Sans);
  for(unsigned int i = 0; i < ans.size(); i++) {
    if(ans[i] == 0) PRINTF("Element %i was not processed\n", i);
    SansPtr[i] = ans[i]==2 ? 1 : 0;
  }
  UNPROTECT(1);
  return(Sans);
}

void nimbleGraph::setNodes(const vector<int> &edgesFrom, const vector<int> &edgesTo,
		     const vector<int> &edgesFrom2ParentExprIDs,
		     const vector<NODETYPE> &types,
		     const vector<string> &names,
		     int inputNumNodes) {
  if(inputNumNodes < 0) PRINTF("Error in setNodes: inputNumNodes < 0\n");
  numNodes = static_cast<unsigned int>(inputNumNodes);
  unsigned int numEdges = edgesFrom.size();

#ifdef _DEBUGNIMGRAPH
  PRINTF("numNodes %i\n", numNodes);
  PRINTF("numEdges %i\n", numEdges);
#endif
  if((numEdges != edgesTo.size()) | (numEdges != edgesFrom2ParentExprIDs.size()) | (numNodes != types.size()) | (numNodes != names.size())) {
    PRINTF("Something is not the right size\n");
    return;
  }
  graphNodeVec.resize(numNodes);
  for(unsigned int iNode = 0; iNode < numNodes; iNode++) {
    graphNodeVec[iNode] = new graphNode(iNode, types[iNode], names[iNode]);
  }
  for(unsigned int iEdge = 0; iEdge < numEdges; iEdge++) {
    graphNodeVec[ edgesFrom[iEdge]]->addChild( graphNodeVec[edgesTo[iEdge]], edgesFrom2ParentExprIDs[iEdge] );
  }
}

vector<int> nimbleGraph::anyStochDependencies() {
  vector<int> ans(numNodes, 0);
  for(unsigned int i = 0; i < numNodes; i++) {
    anyStochDependenciesOneNode(ans, i);
  }
  return(ans);
}

bool nimbleGraph::anyStochDependenciesOneNode(vector<int> &anyStochDependencies,  int CgraphID) {
  // 0 = untouched, 1 = false, 2 = true
  if(anyStochDependencies[CgraphID] != 0) return(anyStochDependencies[CgraphID] == 2);
  bool thisHasAstochDep(false);
  graphNode *thisGraphNode = graphNodeVec[CgraphID];
  graphNode *thisChildNode;
  int numChildren = thisGraphNode->numChildren;
  /* If no children, answer is false */
  if(numChildren == 0) {
    anyStochDependencies[CgraphID] = 1;
    return(false);
  }
  int i(0);
  /* Check type of children without recursing.  If any are STOCH, answer is true */
  while((i < numChildren) & (!thisHasAstochDep)) {
    if(thisGraphNode->children[i]->type == STOCH) {
      thisHasAstochDep = true;
    }
    i++;
  }
  /* If answer was true, we're done */
  if(thisHasAstochDep) {
    anyStochDependencies[CgraphID] = 2;
    return(true);
  }
  /* all children were not STOCH, so now recurse through children */
  i = 0;
  while((i < numChildren) & (!thisHasAstochDep)) {
    thisChildNode = thisGraphNode->children[i];
    if(anyStochDependenciesOneNode(anyStochDependencies, thisChildNode->CgraphID)) {
      thisHasAstochDep = true;
    }
    i++;
  }
  if(thisHasAstochDep) {
    anyStochDependencies[CgraphID] = 2;
    return(true);
  }
  anyStochDependencies[CgraphID] = 1;
  return(false);
}


vector<int> nimbleGraph::anyStochParents() {
  vector<int> ans(numNodes, 0);
  for(int i = static_cast<int>(numNodes - 1); i >= 0; i--) {
    anyStochParentsOneNode(ans, i);
  }
  return(ans);
}

bool nimbleGraph::anyStochParentsOneNode(vector<int> &anyStochParents,  int CgraphID) {
  // 0 = untouched, 1 = false, 2 = true
  if(anyStochParents[CgraphID] != 0) return(anyStochParents[CgraphID] == 2);
  bool thisHasAstochParent(false);
  graphNode *thisGraphNode = graphNodeVec[CgraphID];
  graphNode *thisParentNode;
  int numParents = thisGraphNode->parents.size();
  if(numParents == 0) {
    anyStochParents[CgraphID] = 1;
    return(false);
  }
  int i(0);
  while((i < numParents) & (!thisHasAstochParent)) {
    if(thisGraphNode->parents[i]->type == STOCH) {
      thisHasAstochParent = true;
    }
    i++;
  }
  if(thisHasAstochParent) {
    anyStochParents[CgraphID] = 2;
    return(true);
  }
  i = 0;
  while((i < numParents) & (!thisHasAstochParent)) {
    thisParentNode = thisGraphNode->parents[i];
    if(anyStochParentsOneNode(anyStochParents, thisParentNode->CgraphID)) {
      thisHasAstochParent = true;
    }
    i++;
  }
  if(thisHasAstochParent) {
    anyStochParents[CgraphID] = 2;
    return(true);
  }
  anyStochParents[CgraphID] = 1;
  return(false);
}

//#define _DEBUG_GETDEPS

vector<int> nimbleGraph::getDependencies(const vector<int> &Cnodes, const vector<int> &Comit, bool downstream) {
  // assume on entry that touched = false on all nodes
  // Cnodes and Comit are C-indices (meaning they start at 0)
  int n = Comit.size();
  int i;
  vector<int> ans;
  // touch omit nodes
#ifdef _DEBUG_GETDEPS
  int iDownstream = static_cast<int>(downstream);
  PRINTF("debugging output for getDependencies with %i nodes, %i omits, and downstream = %i.  C indices (graphIDs) shown are 0-based\n", Cnodes.size(), Comit.size(), iDownstream);
#endif
  for(i = 0; i < n; i++) {
    graphNodeVec[ Comit[i] ]->touched = true;
#ifdef _DEBUG_GETDEPS
    PRINTF("touching %i to omit\n", Comit[i]);
#endif
  }
  n = Cnodes.size();
  graphNode *thisGraphNode;
  int thisGraphNodeID;
  for(i = 0; i < n; i++) {
    thisGraphNodeID = Cnodes[i];
    thisGraphNode = graphNodeVec[ thisGraphNodeID ];
#ifdef _DEBUG_GETDEPS
    PRINTF("Working on input node %i\n", thisGraphNodeID);
#endif
    if(!thisGraphNode->touched) {
#ifdef _DEBUG_GETDEPS
      PRINTF("  Adding node %i to ans and recursing\n", thisGraphNodeID);
#endif
      ans.push_back(thisGraphNodeID);
      thisGraphNode->touched = true;
      getDependenciesOneNode(ans, thisGraphNodeID, downstream, 1);
    } else {
#ifdef _DEBUG_GETDEPS
      PRINTF("  Node %i was already touched.\n", thisGraphNodeID);
#endif
      if((thisGraphNode->type == STOCH) & !downstream) {
	/* In this case the input node was already touched, so it is a dependency */
	/* of something earlier on the input list.  But since it was on the input list */
	/* we still need to get its dependencies.  But if downstream is TRUE (==1), then */
	/* its dependencies will have already been pursued so we don't need to */
#ifdef _DEBUG_GETDEPS
      PRINTF("  But is stochastic and downstream is false, so we are recursing into its dependencies.\n");
#endif
	getDependenciesOneNode(ans, thisGraphNodeID, downstream, 1);
      }
    }
  }

  // untouch nodes and omit
  n = Comit.size();
  for(i = 0; i < n; i++) {
    graphNodeVec[ Comit[i] ]->touched = false;
  }
  n = ans.size();
  for(i = 0; i < n; i++) {
    graphNodeVec[ ans[i] ]->touched = false;
  }
  std::sort(ans.begin(), ans.end());
  return(ans);
}

void nimbleGraph::getDependenciesOneNode(vector<int> &deps, int CgraphID, bool downstream, unsigned int recursionDepth) {
  if(recursionDepth > graphNodeVec.size()) {
    PRINTF("ERROR: getDependencies has recursed too far.  Something must be wrong.\n");
    return;
  }
#ifdef _DEBUG_GETDEPS
  PRINTF("    Entering recursion for node %i\n", CgraphID);
#endif
  graphNode *thisGraphNode = graphNodeVec[CgraphID];
  int numChildren = thisGraphNode->numChildren;
  int i(0);
  graphNode *thisChildNode;
  int thisChildCgraphID;
#ifdef _DEBUG_GETDEPS
  PRINTF("      Starting to iterate through %i children of node %i\n", numChildren, CgraphID);
#endif
  for(; i < numChildren; i++) {
    thisChildNode = thisGraphNode->children[i];
    if(thisChildNode->touched) continue;
    thisChildCgraphID = thisChildNode->CgraphID;
#ifdef _DEBUG_GETDEPS
    PRINTF("        Adding child node %i\n", thisChildCgraphID);
#endif
    deps.push_back(thisChildNode->CgraphID);
    thisChildNode->touched = true;
    if(downstream | (thisChildNode->type != STOCH)) {
#ifdef _DEBUG_GETDEPS
    PRINTF("          Recursing into child node %i\n", thisChildCgraphID);
#endif
      getDependenciesOneNode(deps, thisChildCgraphID, downstream, recursionDepth + 1);
    }
  }
#ifdef _DEBUG_GETDEPS
  PRINTF("      Done iterating through %i children of node %i\n", numChildren, CgraphID);
#endif
}

//#define _DEBUG_TOPSORT

// vector<int> nimbleGraph::topologicalSortOrder() { // first try
//   // newCindices will an ordering so oldNodeIDs[Cindices] is validly sorted

//   // things to check: for(j = 0; j < 0; ++j) behavior.
//   // create some kind of failsafe like n^2 passes over everything,
//   // or no new nodes added from a pass.
  
  
//   // then go through their child nodes, check if all parents are touched.
//   //

//   vector<int> newCindices;
  
//   vector<graphNode*> nodesToCheck, newNodesToCheck;

//   // first take nodes with no parents.
//   int i, j;
// #ifdef _DEBUG_TOPSORT
//   std::cout<<"COLLECTING NO-PARENT NODES:\n";
// #endif

//   for(i = 0; i < numNodes; ++i) {
//     if(graphNodeVec[i]->parents.size() == 0) {
// #ifdef _DEBUG_TOPSORT
//       std::cout<<"  Adding no-parent node "<<graphNodeVec[i]->name<<"\n";
//       std::cout<<"     Adding children for next round:";
// #endif
//       newCindices.push_back(graphNodeVec[i]->CgraphID); // should be same as i
//       graphNodeVec[i]->touched = true;
//       for(j = 0; j < graphNodeVec[i]->numChildren; ++j) {
// 	// impossible for children to be touched, because they must have parents
// 	nodesToCheck.push_back(graphNodeVec[i]->children[j]);
// #ifdef _DEBUG_TOPSORT
// 	std::cout<<" "<<graphNodeVec[i]->children[j]->name;
// #endif
	
//       }
// #ifdef _DEBUG_TOPSORT
//       std::cout<<"\n";
// #endif

//     }
//   }

//   bool done(false), allParentsTouched;

//   int numNodesToCheck, numParentsThisNode;
//   graphNode *thisNodeToCheck;
//   while(!done) {
//     // iterate over nodesToCheck
//     numNodesToCheck = nodesToCheck.size();
// #ifdef _DEBUG_TOPSORT
//   std::cout<<"NEW ITERATION TO CHECK "<<numNodesToCheck<<" NODES\n";
// #endif
//     for(i = 0; i < numNodesToCheck; ++i) {
//       // check one node
//       allParentsTouched = true;
//       thisNodeToCheck = nodesToCheck[i];
// #ifdef _DEBUG_TOPSORT
//       std::cout<<"  "<<i<<": checking "<<thisNodeToCheck->name <<"\n";
// #endif

//       if(thisNodeToCheck->touched) {
// #ifdef _DEBUG_TOPSORT
// 	std::cout<<"    Already touched\n";
// #endif
// 	continue; // conceivable if it was first added to newNodesToCheck and then later touched during last iteration
// 	// actually now this shouldn't happen, but we'll keep it in to be safe and sure.
//       }
//       if(thisNodeToCheck->touched2) {
// #ifdef _DEBUG_TOPSORT
// 	std::cout<<"    Already checked in this iteration\n";      
// #endif
// 	continue;
//       }
//       thisNodeToCheck->touched2 = true; // use touched2 to flag that we've already checked this node in this pass
//       j = 0;
//       numParentsThisNode = thisNodeToCheck->parents.size();
//       // see if all parents were touched
//       while(allParentsTouched && j < numParentsThisNode) {
// 	if(!(thisNodeToCheck->parents[j]->touched)) allParentsTouched = false;
// 	// touched2 means it was touched in this round and shouldn't be counted for this child
// 	// simply because that is the ordering choice we make 
// 	else if(thisNodeToCheck->parents[j]->touched2) allParentsTouched = false;
// 	++j;
//       }
//       if(allParentsTouched) { // ok to add to sorted order and include children for next round
// #ifdef _DEBUG_TOPSORT
// 	std::cout<<"    All parents touched TRUE. Adding to next round:";
// #endif

// 	newCindices.push_back(thisNodeToCheck->CgraphID);
// 	thisNodeToCheck->touched = true; 
// 	for(j = 0; j < thisNodeToCheck->numChildren; ++j) {
// 	  // impossible for children to be touched, because this parent was just touched
// #ifdef _DEBUG_TOPSORT
// 	  std::cout<<" "<<thisNodeToCheck->children[j]->name;
// #endif

// 	  newNodesToCheck.push_back(thisNodeToCheck->children[j]);
// 	}
// #ifdef _DEBUG_TOPSORT
// 	std::cout<<"\n";
// #endif

//       } else { // not all parents touched, so add same node to next round
// #ifdef _DEBUG_TOPSORT
// 	std::cout<<"    All parents touched FALSE. Adding self to next round.\n";
// #endif
// 		newNodesToCheck.push_back(thisNodeToCheck);
// 	// NOT necessary to re-check a node that failed this time,
// 	// unless it is the child of something that just passed.
//       }
//     }
//     for(i = 0; i < numNodesToCheck; ++i) {
//       nodesToCheck[i]->touched2 = false;
//     }
//     if(newNodesToCheck.size() == 0) {
//       done = true;
//     } else {
//       nodesToCheck = newNodesToCheck;
//       newNodesToCheck.clear();
//     }
//   }
  
//   for(i = 0; i < numNodes; ++i) {
//     graphNodeVec[i]->touched = false;
//   }
//   return(newCindices);
// }

// A second try, directly imitating igraph_topological_sorting from structural_properties.c from igraph source code

#include<queue>

vector<int> nimbleGraph::topologicalSortOrder() {
  vector<int> newCindices;
  std::queue<int> indicesToTouch;
  vector<int> degreeIn(numNodes);
  vector<int> sortedChildGraphIDs;
  int i, j, CgraphID;
  // all indexing must use CgraphID
  // typically this will match order of graphNodeVec but we will not assume it
  vector<int> CgraphID2index(numNodes); // perhaps this should be maintained at nimbleGraph class level
  for(i = 0; i < numNodes; ++i) {
    CgraphID =  graphNodeVec[i]->CgraphID;
    CgraphID2index[CgraphID] = i;
    if(CgraphID >= numNodes) std::cout<<"Error in sort: A CgraphID is too big.\n";
    degreeIn[i] = graphNodeVec[i]->parents.size();
    if(degreeIn[i] == 0) indicesToTouch.push(i); 
  }
  while(indicesToTouch.size() > 0) {
    int nextIndexToTouch = indicesToTouch.front();
    indicesToTouch.pop();
    newCindices.push_back( graphNodeVec[nextIndexToTouch]->CgraphID );
    degreeIn[nextIndexToTouch] = -1;
    int numChildren = graphNodeVec[nextIndexToTouch]->numChildren;
    sortedChildGraphIDs.resize(numChildren);
    for(j = 0; j < numChildren; ++j) {
      sortedChildGraphIDs[j] = graphNodeVec[nextIndexToTouch]->children[j]->CgraphID;
      // BUILD A VECTOR OF CgraphIDs and sort them before inspecting and pushing them
      // This should match sorting done in igraph_neighbors in type_indexededgelist.c in igraph source
    }
    std::sort(sortedChildGraphIDs.begin(), sortedChildGraphIDs.end() );

    for(j = 0; j < numChildren; ++j) {
      //      int childIndex = CgraphID2index[ graphNodeVec[nextIndexToTouch]->children[j]->CgraphID ];
      int childIndex = CgraphID2index[ sortedChildGraphIDs[j] ];
      if(--(degreeIn[ childIndex ]) == 0) {
	indicesToTouch.push( childIndex );
      }
    }
  }
  if(newCindices.size() < numNodes) {
    std::cout<<"Error from topological sorting: graph seems to have a cycle.\n";
  }
  return(newCindices);
}
