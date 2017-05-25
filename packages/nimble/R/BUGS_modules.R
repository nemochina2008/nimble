
## This takes in a code-expander function and packs it up as a recognizable BUGS module
## For now, all it does is make a list with a class attribute.
## It could provide more general features later.
makeBUGSmodule <- function( fun ) {
    
    ans <- structure(list(process = fun), class = "BUGSmodule")
    ans
}

## 
## This is an example of a code-expanding function
## We could either allow data to be passed in or provide other ways to label factors and lengths
## For now the contents of this are a proxy -- nothing real is done
lmPred <- makeBUGSmodule(
    function(LHS, RHS) { ##, stoch, family=c("Normal"), link = c('Identity'), priors=c("Normal","t","Uniform"), dropbase=TRUE){ ## RHS will entire RHS code, so expand args by match.call

        ## could provide a utility function or entry step to do following every time a  module is called.  Difficulty would be knowing what is supposed to be code.
        RHSargs <- match.call(function(RHSmodel, stoch, family, link, priors, dropbase=TRUE) {}, RHS)
        link <- match.arg(RHSargs[['link']], c('Identity','log'))
        
    ## as a proxy, suppose LHS is predY[1:n] and RHS is X[1:n]
        ## Then the constructed code would be:
        ## (I am splitting link in a cheap way to illustrate glm later
        ##  but in the real system the link would be constructed while handling LHS.)
        if(link == 'Identity') {
            newCode <- quote({ ## This would be constructed based on LHS, RHS, and other arguments
                X_coeff ~ dnorm(0, sd = 1000)
                for(i in 1:n) {
                    predY[i] <- X_coeff * X[i] 
                }
            })
        } else if (link == 'log') {
            newCode <- quote({ ## This would be constructed based on LHS, RHS, and other arguments
                X_coeff ~ dnorm(0, sd = 1000)
                for(i in 1:n) {
                    log(predY[i]) <- X_coeff * X[i] 
                }
            })
        }
    return(list(code = newCode)) ## In a list so we can add more return content later
})

nim_lm <- makeBUGSmodule(
    function(LHS, RHS) { ## skip match.call on RHS for now

    ## as a proxy, suppose LHS is predY[1:n] and RHS is X[1:n]
    ## Then the constructed code would be:
    newCode <- quote({ ## This would be constructed based on LHS, RHS, and other arguments
        for(i in 1:n) {
            Y[i] ~ dnorm(predY[i], sd = sigmga)
        }
        predY[1:n] ~ lmPred(X[1:n], priors = priors) ## rely on recursion to lmPred
    })
    return(list(code = newCode)) ## In a list so we can add more return content later
})

nim_glm <- makeBUGSmodule(
    function(LHS, RHS) { ## skip match.call on RHS for now

    ## as a proxy, suppose LHS is predY[1:n] and RHS is X[1:n]
    ## Then the constructed code would be:
    newCode <- quote({ ## This would be constructed based on LHS, RHS, and other arguments
        for(i in 1:n) {
            Y[i] ~ dexp(predY[i])
        }
        predY[1:n] ~ lmPred(X[1:n], priors = priors, link = 'log') ## rely on recursion to lmPred
    })
    return(list(code = newCode)) ## In a list so we can add more return content later
})

