
## This takes in a code-expander function and packs it up as a recognizable BUGS module
## For now, all it does is make a list with a class attribute.
## It could provide more general features later.
makeBUGSmodule <- function( fun ) {
    
    ans <- structure(list(fun), class = "BUGSmodule")
    ans
}

## Some borrowed and modified code from nimbleEcology repo
## This is an example of a code-expanding function
## We could either allow data to be passed in or provide other ways to label factors and lengths
lm <- function(LHS, RHS, stoch, family=c("Normal"), priors=c("Normal","t","Uniform"), dropbase=TRUE){

    ## The match.call, mf, etc. would belong in R function nimble.lm
  ##   cl <- match.call()
  ## mf <- match.call(expand.dots = FALSE)
  
  ## m <- match(c("mod","dat"), names(mf), 0L)
  
  ## ## construct a call
  ## mf <- mf[c(1L, m)] 
  ## names(mf)[names(mf)=='dat'] <- 'data'
  ## names(mf)[names(mf)=='mod'] <- 'formula'
  ## mf$drop.unused.levels <- TRUE
  
  ## ## change name of function
  ## mf[[1L]] <- quote(stats::model.frame)
  
  ## ## and eval it in calling environment
  ## mf <- eval(mf, parent.frame())
  ## ## This gets us mf with data filled in from dat argument or from calling environment as needed
  
  ## #Convert Character Variables to Factors
  ## for (i in 1:dim(mf)[2]){
  ##   mf[,i] <- makeFac(mf[,i],char.only = TRUE)
  ## }
  
  
  ## #Check link and prior arguments are supported
  ## family = match.arg(family)
  ## prior = match.arg(priors)
  
  ## mod = as.formula(mod)

    LHSvar <- LHS[[2]] ## x of x[1:n]
    LHSlength <- LHS[[3]][[3]] ## n of x[1:n] assuming to be quick that it starts at 1
    
  index.name <- quote(i)
   
  #Likelihood
  if (family=="Normal"){
    meanfcn <- quote(dnorm(mu[i],tau))}  #could change to other disributions later on
    dep <- LHS2BUGSterm(LHSvar, index.name)
    meanfctn <-  substitute(LHSvar ~ RHS, list(LHS = dep, RHS = meanfcn))

  
    #Generate mean function mu[i]
    #Intercept
    int.to.add <- attr(terms(mod), 'intercept')
  
    terms.to.add <- attr(terms(mod), 'term.labels')
    assert_that(length(terms.to.add) > 0)
  
  
  ## Assume Fixed Intercept for Now
  if (int.to.add==1){
      RHS <- quote(b0)
      ## append additional terms
      for(iTerm in seq_along(terms.to.add)) {
        RHS <- addTerm(RHS, formulaTerm2BUGSterm(terms.to.add[iTerm], mf, index.name,id = iTerm))
    }
  }
  
  if(int.to.add==0){
    RHS <- formulaTerm2BUGSterm(terms.to.add[1], mf, index.name, id=1)
    for(iTerm in seq_along(terms.to.add[-1])+1) {
      RHS <- addTerm(RHS, formulaTerm2BUGSterm(terms.to.add[iTerm], mf, index.name,id = iTerm))
    }
  }

  
  #Make Term for mean function, mu  
  LHS <- LHS2BUGSterm(quote(mu), index.name)
  fullLine <- substitute(LHS ~ RHS, list(LHS = LHS, RHS = RHS))
  lines <- list(meanfctn, fullLine) ## could have potentially more
  forLoop <- embedLinesInForLoop(lines, index.name, start = 1, finish = as.numeric(dim(mf)[1])) ## num.data will have to be figured out
  
  
  #Specify Priors
  priorsList <- list()
  
  #Prior for Variance Term, could allow options for this later
  priorsList[[1]] <-  quote(tau ~ dgamma (0.001, 0.001))
  priorsList[[2]] <-  quote(sigma2 <- 1/tau)
  
  
  #Priors for Beta Terms
  for(iTerm in seq_along(terms.to.add)) {
    priorsList[[iTerm+2]] <- formulaTerm2PRIOR(terms.to.add[iTerm], mf, index.name,id = iTerm, dropbase, prior = priors)
  }
  
  #Add Prior for Intercept 
  if(int.to.add==1){
    priorsList[[length(priorsList)+1]] <- make.prior(quote(b0),priors)
  }
  
  
  Priors <- embedLinesInCurlyBrackets(priorsList)
  bugs.out <- embedLinesInCurlyBrackets(list(forLoop,Priors))
  
  bugs.out
}




write.abun.sym(abun = z ~ 1 + xm + A, Counts = sillyData, priors=list(xm = 'vague')
