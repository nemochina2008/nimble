source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context('Profiling Nimble code')

test_that('Profiling matrix multiplication', {
    nimFun <- nimbleFunction(
        run = function(x = double(2), y = double(2), numIters = integer(0)) {
            gprofile({
                for (i in 1:numIters) {
                    z <- x %*% y
                }
            })
            return(z)
            returnType(double(2))
        })
    
    nimbleOptions(useGooglePerftools = TRUE)
    nimbleOptions(showCompilerOutput = TRUE)
    cFun <- compileNimble(nimFun)
    N <- 100
    iters <- 100
    x <- matrix(rnorm(N * N), N, N)
    y <- matrix(rnorm(N * N), N, N)
    cFun(x, y, iters)
    expect_true(file.exists(file.path(tempdir(), 'nimble.profile')))
})

