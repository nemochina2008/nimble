source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
context('nimbleList() tests in RC functions')

## All necessary nlDefinitions below:
testListDef1 <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
otherTestListDef1 <- testListDef1
testListDef2 <- nimbleList(nlScalar = double(0), nestedNL = testListDef1())
differentNlDef <- nimbleList(differentNlScalar = double(0))
assign('testListDef1', testListDef1, envir = .GlobalEnv)
assign('otherTestListDef1', otherTestListDef1, envir = .GlobalEnv)
assign('testListDef2', testListDef2, envir = .GlobalEnv)
assign('differentNlDef', differentNlDef, envir = .GlobalEnv)

cat("### nimbleList RCfun Test 1 ###\n")
##
##Test functions 1 and 1a create nimbleLists inside of RC functions.  Test 1a creates an otherTestListDef1, but returns a testListDef1.
##
nlTestFunc1 <- nimbleFunction(
  run = function(){
    doubleMatrix <- diag(1)
    doubleScalar <- 1
    doubleVector <- numeric(2, 1)
    newList1 <- testListDef1$new(nlVector = doubleVector, nlMatrix = doubleMatrix)
    newList1$nlScalar <- doubleScalar
    returnType(testListDef1())
    return(newList1)
  }
)

RnimbleList1 <- nlTestFunc1()
CnlTestFunc1 <- compileNimble(nlTestFunc1)
CnimbleList1 <- CnlTestFunc1()

expect_identical(RnimbleList1$nlScalar, 1)
expect_identical(RnimbleList1$nlMatrix, CnimbleList1$nlMatrix)


nlTestFunc1a <- nimbleFunction(
  run = function(){
    outList <- otherTestListDef1$new(nlScalar = 5)
    returnType(testListDef1())
    return(outList)
  })

RnimbleList1a <- nlTestFunc1a()
CnlTestFunc1a <- compileNimble(nlTestFunc1a)
CnimbleList1a <- CnlTestFunc1a()

expect_identical(RnimbleList1a$nlScalar, 5)
expect_identical(RnimbleList1a$nlScalar, CnimbleList1a$nlScalar)

cat("### nimbleList RCfun Test 2 ###\n")
###
## Test functions 2 and 2a use a nimbleList returning RC function within a nimbleFunction.  In 2a, the names of the nimbleList definitions
## differ between the two functions (but refer to the same definition)
###
testListDef1 <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
otherTestListDef1 <- testListDef1
testListDef2 <- nimbleList(nlScalar = double(0), nestedNL = testListDef1())
differentNlDef <- nimbleList(differentNlScalar = double(0))
assign('testListDef1', testListDef1, envir = .GlobalEnv)
assign('otherTestListDef1', otherTestListDef1, envir = .GlobalEnv)
assign('testListDef2', testListDef2, envir = .GlobalEnv)
assign('differentNlDef', differentNlDef, envir = .GlobalEnv)

innerNlTestFunc2 <- nimbleFunction(
  run = function(){
    newList2 <- testListDef1$new(nlScalar = 2)
    returnType(testListDef1())
    return(newList2)
  }
)
assign('innerNlTestFunc2', innerNlTestFunc2, envir = .GlobalEnv)

outerNlTestFunc2 <- nimbleFunction(
  setup = function(){
  },
  run = function(){
    outList <- innerNlTestFunc2()
    returnType(testListDef1())
    return(outList)
  }
)


nlTestFunc2 <- outerNlTestFunc2()
RnimbleList2 <- nlTestFunc2$run()
CnlTestFunc2 <- compileNimble(nlTestFunc2, showCompilerOutput = TRUE)
CnimbleList2 <- CnlTestFunc2$run()

expect_identical(RnimbleList2$nlScalar, 2)
expect_identical(RnimbleList2$nlScalar, CnimbleList2$nlScalar)

cat("### nimbleList RCfun Test 2a ###\n")
testListDef1 <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
otherTestListDef1 <- testListDef1
testListDef2 <- nimbleList(nlScalar = double(0), nestedNL = testListDef1())
differentNlDef <- nimbleList(differentNlScalar = double(0))
assign('testListDef1', testListDef1, envir = .GlobalEnv)
assign('otherTestListDef1', otherTestListDef1, envir = .GlobalEnv)
assign('testListDef2', testListDef2, envir = .GlobalEnv)
assign('differentNlDef', differentNlDef, envir = .GlobalEnv)

innerNlTestFunc2a <- nimbleFunction(
  run = function(){
    newList2a <- testListDef1$new(nlScalar = .5)
    returnType(testListDef1())
    return(newList2a)
  }
)
assign('innerNlTestFunc2a', innerNlTestFunc2a, envir = .GlobalEnv)

outerNlTestFunc2a <- nimbleFunction(
  setup = function(nlGen){
  },
  run = function(){
    outList <- innerNlTestFunc2a()
    returnType(nlGen())
    return(outList)
  }
)

nlTestFunc2a <- outerNlTestFunc2a(testListDef1)
RnimbleList2a <- nlTestFunc2a$run()
CnlTestFunc2a <- compileNimble(nlTestFunc2a, showCompilerOutput = TRUE)
CnimbleList2a <- CnlTestFunc2a$run()


expect_identical(RnimbleList2a$nlScalar, .5)
expect_identical(RnimbleList2a$nlScalar, CnimbleList2a$nlScalar)

cat("### nimbleList RCfun Test 3 ###\n")
###
## Test functions 3 and 3a take a nimbleList as an argument, copy it, modify it, and return a nimbleList of the same type. 
## 3a uses alternate names for the same nimbleList definition.
###
testListDef1 <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
otherTestListDef1 <- testListDef1
testListDef2 <- nimbleList(nlScalar = double(0), nestedNL = testListDef1())
differentNlDef <- nimbleList(differentNlScalar = double(0))
assign('testListDef1', testListDef1, envir = .GlobalEnv)
assign('otherTestListDef1', otherTestListDef1, envir = .GlobalEnv)
assign('testListDef2', testListDef2, envir = .GlobalEnv)
assign('differentNlDef', differentNlDef, envir = .GlobalEnv)

nlTestFunc3 <- nimbleFunction(
  run = function(nlArg = testListDef1()){
    outList <- testListDef1$new()
    outList <- nlArg
    outList$nlScalar <- nlArg$nlVector[1] + 4.3
    returnType(testListDef1())
    return(outList)
  }
)

nlArg3 <- testListDef1$new(nlVector = c(1,2,3,4))
RnimbleList3 <-nlTestFunc3(nlArg3)
CnlTestFunc3 <- compileNimble(nlTestFunc3)
CnimbleList3 <- CnlTestFunc3(nlArg3)

expect_identical(RnimbleList3$nlVector, c(1,2,3,4))
expect_identical(RnimbleList3$nlScalar, CnimbleList3$nlScalar)

cat("### nimbleList RCfun Test 3a ###\n")
testListDef1 <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
otherTestListDef1 <- testListDef1
testListDef2 <- nimbleList(nlScalar = double(0), nestedNL = testListDef1())
differentNlDef <- nimbleList(differentNlScalar = double(0))
assign('testListDef1', testListDef1, envir = .GlobalEnv)
assign('otherTestListDef1', otherTestListDef1, envir = .GlobalEnv)
assign('testListDef2', testListDef2, envir = .GlobalEnv)
assign('differentNlDef', differentNlDef, envir = .GlobalEnv)

nlTestFunc3a <- nimbleFunction(
  run = function(nlArg = testListDef1()){
    outList <- otherTestListDef1$new()
    outList <- nlArg
    outList$nlScalar <- nlArg$nlVector[1] + 4.3
    returnType(testListDef1())
    return(outList)
  }
)

nlArg3 <- testListDef1$new(nlVector = c(1,2,3,4))
RnimbleList3a <- nlTestFunc3a(nlArg3)
CnlTestFunc3a <- compileNimble(nlTestFunc3a)
CnimbleList3a <- CnlTestFunc3a(nlArg3)

expect_identical(RnimbleList3a$nlVector, c(1,2,3,4))
expect_identical(RnimbleList3a$nlScalar, CnimbleList3a$nlScalar)

cat("### nimbleList RCfun Test 4 ###\n")
###
## Test function 4 creates a nimbleList with a nested nimbleList inside and RC function
###
testListDef1 <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
otherTestListDef1 <- testListDef1
testListDef2 <- nimbleList(nlScalar = double(0), nestedNL = testListDef1())
differentNlDef <- nimbleList(differentNlScalar = double(0))
assign('testListDef1', testListDef1, envir = .GlobalEnv)
assign('otherTestListDef1', otherTestListDef1, envir = .GlobalEnv)
assign('testListDef2', testListDef2, envir = .GlobalEnv)
assign('differentNlDef', differentNlDef, envir = .GlobalEnv)

nlTestFunc4 <- nimbleFunction(
  run = function(){
    outList <- testListDef2$new()
    outList$nestedNL <- otherTestListDef1$new()
    outList$nlScalar <- 2
    outList$nestedNL$nlScalar <- outList$nlScalar
    returnType(testListDef2())
    return(outList)
  }
)

RnimbleList4 <- nlTestFunc4()
CnlTestFunc4 <- compileNimble(nlTestFunc4, showCompilerOutput = TRUE)
CnimbleList4 <- CnlTestFunc4()

expect_identical(RnimbleList4$nlScalar,2)
expect_identical(RnimbleList4$nlScalar, CnimbleList4$nlScalar)


cat("### nimbleList RCfun Test 5 ###\n")
###
## Test function 5 takes 2 nimbleLists as arguments - the type of the second nl argument is nested within the first nl argument.
###
testListDef1 <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
otherTestListDef1 <- testListDef1
testListDef2 <- nimbleList(nlScalar = double(0), nestedNL = testListDef1())
differentNlDef <- nimbleList(differentNlScalar = double(0))
assign('testListDef1', testListDef1, envir = .GlobalEnv)
assign('otherTestListDef1', otherTestListDef1, envir = .GlobalEnv)
assign('testListDef2', testListDef2, envir = .GlobalEnv)
assign('differentNlDef', differentNlDef, envir = .GlobalEnv)

nlTestFunc5 <- nimbleFunction(
  run = function(argList = testListDef2(), argList2 = otherTestListDef1()){
    outList <- argList
    outList$nestedNL <- argList2
    outList$nlScalar <- 2
    outList$nestedNL$nlScalar <- argList2$nlScalar
    returnType(testListDef2())
    return(outList)
  }
)

nlArg5 <- testListDef2$new(nlScalar = 3.14)
nlArg5_2 <- testListDef1$new(nlScalar = 7)

RnimbleList5 <- nlTestFunc5(nlArg5, nlArg5_2)
CnlTestFunc5 <- compileNimble(nlTestFunc5)
CnimbleList5 <- CnlTestFunc5(nlArg5, nlArg5_2)

expect_identical(RnimbleList5$nestedNL$nlScalar, 7)
expect_identical(RnimbleList5$nestedNL$nlScalar, CnimbleList5$nestedNL$nlScalar)


cat("### nimbleList RCfun Test 6 ###\n")
###
## Test functions 6 and 7 test the use of eigen (and eigenNimbleLists) and svd (and svdNimbleLists) inside RC functions
###
testListDef1 <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
otherTestListDef1 <- testListDef1
testListDef2 <- nimbleList(nlScalar = double(0), nestedNL = testListDef1())
differentNlDef <- nimbleList(differentNlScalar = double(0))
assign('testListDef1', testListDef1, envir = .GlobalEnv)
assign('otherTestListDef1', otherTestListDef1, envir = .GlobalEnv)
assign('testListDef2', testListDef2, envir = .GlobalEnv)
assign('differentNlDef', differentNlDef, envir = .GlobalEnv)

nlTestFun6 <- nimbleFunction(
  run = function(){
    simpleMat <- diag(4)
    eigenList <- eigen(simpleMat)
    eigenList$values[1] <- 0
    returnType(eigenNimbleList())
    return(eigenList)
  }
)

RnimbleList6 <- nlTestFun6()
CnlTestFunc6 <- compileNimble(nlTestFun6)
CnimbleList6 <- CnlTestFunc6()

expect_identical(RnimbleList6$values[1], 0)
expect_identical(RnimbleList6$vectors, CnimbleList6$vectors)

cat("### nimbleList RCfun Test 7 ###\n")
testListDef1 <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
otherTestListDef1 <- testListDef1
testListDef2 <- nimbleList(nlScalar = double(0), nestedNL = testListDef1())
differentNlDef <- nimbleList(differentNlScalar = double(0))
assign('testListDef1', testListDef1, envir = .GlobalEnv)
assign('otherTestListDef1', otherTestListDef1, envir = .GlobalEnv)
assign('testListDef2', testListDef2, envir = .GlobalEnv)
assign('differentNlDef', differentNlDef, envir = .GlobalEnv)

nlTestFun7 <- nimbleFunction(
  run = function(svdArg = svdNimbleList()){
    simpleMat <- diag(4)
    svdList <- svd(simpleMat)
    svdList$d <- svdArg$d
    returnType(svdNimbleList())
    return(svdList)
  }
)

svdArgList <- svdNimbleList$new(d = c(4, 3, 2, 1))
RnimbleList7 <- nlTestFun7(svdArgList)
CnlTestFunc7 <- compileNimble(nlTestFun7)
CnimbleList7 <- CnlTestFunc7(svdArgList)

expect_identical(RnimbleList7$d, c(4, 3, 2, 1))
expect_identical(RnimbleList7$d, CnimbleList7$d)

cat("### nimbleList RCfun Test 8 ###\n")
###
## Test function 8 has 2 nimbleList arguments coming from different nlDefs.  The second nlDef will not be added to the symbolTable unless
## it is detected as an argument type - this is currently broken.
###
testListDef1 <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
otherTestListDef1 <- testListDef1
testListDef2 <- nimbleList(nlScalar = double(0), nestedNL = testListDef1())
differentNlDef <- nimbleList(differentNlScalar = double(0))
assign('testListDef1', testListDef1, envir = .GlobalEnv)
assign('otherTestListDef1', otherTestListDef1, envir = .GlobalEnv)
assign('testListDef2', testListDef2, envir = .GlobalEnv)
assign('differentNlDef', differentNlDef, envir = .GlobalEnv)

nlTestFun8 <- nimbleFunction(
  run = function(nlArg1 = testListDef1(), nlArg2 = differentNlDef()){
    nlArg1$nlScalar <- nlArg2$differentNlScalar
    returnType(testListDef1())
    return(nlArg1)
  }
)

nlArg8 = testListDef1$new(nlScalar = .1)
nlArg8_2 = differentNlDef$new(differentNlScalar = .2)
RnimbleList8 <- nlTestFun8(nlArg8, nlArg8_2)
CnlTestFunc8 <- compileNimble(nlTestFun8)
CnimbleList8 <- CnlTestFunc8(nlArg8, nlArg8_2)

expect_identical(RnimbleList8$nlScalar, .2)
expect_identical(RnimbleList8$nlScalar, CnimbleList8$nlScalar)

cat("### nimbleList RCfun Test 9 ###\n")
###
## This ninth set of test functions has two levels of nested functions that return nimbleLists generated from the same definition.
## We then compile another function within the same project and use the same nimbleList generator there.  It also reuses the nlArg8 nimbleList.
###

testListDef1 <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
otherTestListDef1 <- testListDef1
testListDef2 <- nimbleList(nlScalar = double(0), nestedNL = testListDef1())
differentNlDef <- nimbleList(differentNlScalar = double(0))
assign('testListDef1', testListDef1, envir = .GlobalEnv)
assign('otherTestListDef1', otherTestListDef1, envir = .GlobalEnv)
assign('testListDef2', testListDef2, envir = .GlobalEnv)
assign('differentNlDef', differentNlDef, envir = .GlobalEnv)

nlTestFun9 <- nimbleFunction(
  run = function(nlArg1 = testListDef1(), doubleArg = double()){
    returnNl <- innerNlTestFun9_1(nlArg1, doubleArg)
    returnNl$nlMatrix <- eigen(diag(2))$vectors
    returnType(testListDef1())
    return(returnNl)
  }
)
assign('nlTestFun9', nlTestFun9, envir = .GlobalEnv)
    
innerNlTestFun9_2 <- nimbleFunction(
  run = function(nlArg1 = testListDef1(), doubleArg = double()){
    nlArg1$nlScalar <- testListDef1$new(nlScalar = -1)$nlScalar + doubleArg
    dummyNl <- testListDef1$new(nlVector = c(1,2,3))
    nlArg1$nlVector <- dummyNl$nlVector
    returnType(testListDef1())
    return(nlArg1)
  }
)
assign('innerNlTestFun9_2', innerNlTestFun9_2, envir = .GlobalEnv)

innerNlTestFun9_1 <- nimbleFunction(
  run = function(nlArg1 = testListDef1(), doubleArg = double()){
    nlArg1$nlScalar <- innerNlTestFun9_2(nlArg1, doubleArg)$nlScalar
    returnType(testListDef1())
    return(nlArg1)
  }
)

assign('innerNlTestFun9_1', innerNlTestFun9_1, envir = .GlobalEnv)


nlArg9 <- testListDef1$new(nlScalar = .1)
#CnlTestList9 <- compileNimble(list(innerNlTestFun9_2, innerNlTestFun9_1,  nlTestFun9))
#CnimbleList9 <- CnlTestList9[[3]](nlArg9, 10)
#expect_identical(CnimbleList9$nlScalar, 9)

wrapperNlTestFun9 <- nimbleFunction(
  setup = TRUE, 
  run = function(doubleArg = double()){
    newNl <- testListDef1$new()
    newNl <- nlTestFun9(newNl, doubleArg)
    returnType(testListDef1())
    return(newNl)
  }
)

wrapperTestFun9Instance <- wrapperNlTestFun9()
CwrapperTestFun9Instance <- compileNimble(wrapperTestFun9Instance)
CnimbleList9_2 <- CwrapperTestFun9Instance$run(4)
expect_identical(CnimbleList9_2$nlScalar, 3)

wrapperTestFun9Instance_2 <- wrapperNlTestFun9()
CwrapperTestFun9Instance_2 <- compileNimble(wrapperTestFun9Instance_2, project = wrapperTestFun9Instance)
CnimbleList9_2_2 <- CwrapperTestFun9Instance_2$run(0)
expect_identical(CnimbleList9_2_2$nlScalar, -1)


cat("### nimbleList RCfun Test 10 ###\n")
###
## Test 10 calls a nimbleFunction that uses the eigen() function.  This function takes an eigenNimbleList as a
## setup argument.  The function is compiled 3 times, with 3 different nl setup arguments.
###
testListDef1 <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
otherTestListDef1 <- testListDef1
testListDef2 <- nimbleList(nlScalar = double(0), nestedNL = testListDef1())
differentNlDef <- nimbleList(differentNlScalar = double(0))
assign('testListDef1', testListDef1, envir = .GlobalEnv)
assign('otherTestListDef1', otherTestListDef1, envir = .GlobalEnv)
assign('testListDef2', testListDef2, envir = .GlobalEnv)
assign('differentNlDef', differentNlDef, envir = .GlobalEnv)

nlTestFun10 <- nimbleFunction(
  setup = function(eigList){}, 
  run = function(matrixArg = double(2)){
    eigList$values <<- eigen(matrixArg, TRUE)$values
    returnType(eigenNimbleList())
    return(eigList)
  }
)

setupEigenList <- eigenNimbleList$new()
nlFunc10Instance <- nlTestFun10(setupEigenList)
CnlFunc10Instance <- compileNimble(nlFunc10Instance)
CnimbleList10 <- CnlFunc10Instance$run(diag(2))
expect_identical(CnimbleList10$values, c(1,1))

setupEigenList_2 <- eigenNimbleList$new()
nlFunc10Instance_2 <- nlTestFun10(setupEigenList_2)
CnlFunc10Instance_2 <- compileNimble(nlFunc10Instance_2, project = nlFunc10Instance) ## same project
CnimbleList10_2 <- CnlFunc10Instance_2$run(2*diag(2))
expect_identical(CnimbleList10_2$values, c(2,2))

setupEigenList_3 <- eigenNimbleList$new()
nlFunc10Instance_3 <- nlTestFun10(setupEigenList_3)
CnlFunc10Instance_3 <- compileNimble(nlFunc10Instance_3)
CnimbleList10_3 <- CnlFunc10Instance_3$run(3*diag(2))
expect_identical(CnimbleList10_3$values, c(3,3))


