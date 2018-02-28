## this is a set of functions that link across sSIR and HawkesN
library(data.table)
library(Matrix)

#' this function the probability distribution of the final size of a diffusion
#' cascade, given an observed initial prefix.
#' the params are HawkesN parameters
get.size.distribution <- function(params = c(K = 5, c = 0.001, theta = 0.2, N = 400), history = NULL, .transition = NULL) {
  params.S <- HAWKES2SIR.transform.paramters(params = params)
  
  ## get out own calculation of the transition matrix. I know that
  ## `get.final.size.distribution`, but I like things I can see.
  if (is.null(.transition)) {
    .transition <- construct.transition.matrix.and.states(params = params.S)
  }
  
  ## if we were given a history, estimate the number of infected individuals at
  ## the end of the observed history
  if (! is.null(history)) {
    history <- as.data.frame(history)
    maxtime <- max(history$time)
    
    Iest <- round(sum(exp(-params.S["gamma"] * (maxtime - history$time[1:nrow(history)]))))
    Sest <- params.S["N"] - nrow(history)
    
    ## call the function to get the distribution
    myres <- get.final.size.distribution(params = params.S, start.from = list(S = Sest, I = Iest), .transition = .transition)
  } else {
    ## call the function to get the distribution, but no starting point.
    myres <- get.final.size.distribution(params = params.S, .transition = .transition)
  }
  
  return(myres)
}

## this function takes HAWKESN EXPN parameters, and return the equivalent SIR paremeters.
HAWKES2SIR.transform.paramters <- function(params = c(K = 5, c = 0.001, theta = 0.2, N = 400)) {
  params <- unlist(params)
  
  ## first construct the equivalent SIR parameters.
  ## We need to build the expectation over event marks mi
  C <- params["K"]
  params.S <- c(N = round(params["N"]), I.0 = 1, gamma = params["theta"], beta = C*params["theta"]) ## note that I.0 is 1, because in HawkesN I start from one immigrant
  names(params.S) <- c("N", "I.0", "gamma", "beta")
  params.S <- unlist(params.S)
  
  if (sum(!is.finite(params.S)) > 0) {
    stop(sprintf("At least one parameter in the equiv. SIR model is not finite: N = %.2f, I.0 = %.2f, gamma = %.2f, beta = %.2f", 
                 params.S["N"], params.S["I.0"], params.S["gamma"], params.S["beta"]))
  }
  
  return(params.S)
}

#' gets the distribution of the final size of a diffusion. Note that start.from can specify where we start from:
#' start.from <- list(S = 20, I = 2), starts from a susceptible population of 20 and infected of 2, even if N = 100
get.final.size.distribution <- function(params = c(N = 100, I.0 = 40, gamma = 1, beta = 2.3), start.from = NULL, .transition = NULL) {
  names(params) <- c("N", "I.0", "gamma", "beta")
  params <- unlist(params)
  
  ## lets get the theoretical size
  R0 <- params["beta"] / params["gamma"]
  
  f <- function(sinf) {
    # val <- log(sinf) / (sinf - 1) - R0
    val <- params["N"] + (params["N"] / R0 ) * log(sinf / (params["N"] - params["I.0"])) - sinf
    return(val)
  }
  
  ## solve the equation f = 0 on the interval [0, N) -- sinf is the # of infected at \infty
  res <- uniroot(f = f, interval = c(0, params["N"]))
  theo.mean <- params["N"] - res$root
  
  #### now let's get the stochastic distribution
  ## first, the transition matrix and the mapping
  if (!is.null(.transition)) {
    transition <- .transition
  } else {
    transition <- construct.transition.matrix.and.states(params = params)
  }
  ## the initial state
  init.state <- rep(x = 0, time = nrow(transition$stateMapping))
  
  if (is.null(start.from)) {
    start.from <- list(S = (params["N"] - params["I.0"]), I = params["I.0"])
  }
  
  init.state[which(transition$stateMapping$I == start.from$I & transition$stateMapping$S == start.from$S )] <- 1
  
  ## compute the final distribution via iterative multiplication 
  state <- init.state
  for (i in 1:(2*params['N'] - 1)) {
    state <- transition$transitionMatrix %*% state
    # if (i %% 20 == 0) {
      # ggsave(filename = paste0("testani/", i, ".pdf"), plot = p)
    # }
    # p <- drawHeatMap(N = params["N"], transition = transition, state = state, init.state = init.state, max_value = max_value)
    # ggsave(filename = paste0("testani/", i, ".png"), plot = p)
  }
  
  ## get final state and the stochastic mean
  final.state <- state[which(transition$stateMapping$I == 0)]
  names(final.state) <- params["N"] - transition$stateMapping$S[which(transition$stateMapping$I == 0)]
  
  stochastic.mean <- sum(as.numeric(names(final.state)) * final.state)
  
  return(list(theo.mean = theo.mean, stochastic.mean = stochastic.mean, final.state = final.state))
}

#' This function returns the transition matrix and the states mapping matrix for
#' an SIR model with the given parameters. Note that this is a transition matrix
#' that only looks at the next transition, not at the timing of the transitions.
construct.transition.matrix.and.states <- function(params = c(N = 20, I.0 = 1, gamma = 1, beta = 2)) {
  names(params) <- c("N", "I.0", "gamma", "beta")
  params <- unlist(params)
  
  ## first construct the state mapping
  S <- params["N"]:0
  I <- 0:params["N"]
  
  # construct the grid of the two stochastic variables
  state.mapping <- expand.grid(S = S, I = I)
  maxNoStates <- nrow(state.mapping)
    
  ## fill in the transition matrix, starting from each state
  f <- function(mystate) {
    # j <- 10000
    # mystate <- unlist(state.mapping[j,])
    retStack <- matrix(data = NA, nrow = 2, ncol = 3)
    colnames(retStack) <- c("newState", "oldState", "prob")
    
    ## check if we are in a feasable state, otherwise return void
    if (mystate["S"] + mystate["I"] > params["N"])
      return(data.frame(retStack))
    
    j <- (params["N"] + 1) * mystate["I"] + params["N"] + 1 - mystate["S"]
    
    if (mystate["I"] > 0) { 
      ## if on this branch, then we are on a non-absorbing state
      
      ## new recovery -- Linda Allen pag. 30
      newS <- mystate["S"]
      newI <- mystate["I"] - 1
      prob <- params["gamma"] / (params["gamma"] + (params["beta"] / params["N"] ) * newS )
      
      ## compute the new state. This trick is from indexing, because S and I take subsequent values.
      ## it is like counting in base N+1 (where N is max value of S and I)
      newState <- (params["N"] + 1) * newI + params["N"] + 1 - newS
      
      ## this is the old way of getting the next state -- very slow
      # newState <- which(state.mapping$S == newS & state.mapping$I == newI)
      
      if (newState > 0 && newState <= maxNoStates &&
          newI >= 0 && newI <= params["N"] &&
          newS >= 0 && newS <= params["N"]) {
        ## if here, means that the new state is valid
        retStack[1,] <- c(newState, j, prob)
      }
      
      ## new infection
      newS <- mystate["S"] - 1 
      newI <- mystate["I"] + 1 
      
      ## compute the new state. This trick is from indexing, because S and I take subsequent values.
      ## it is like counting in base N+1 (where N is max value of S and I)
      newState <- (params["N"] + 1) * newI + params["N"] + 1 - newS
      
      ## this is the old way of getting the next state
      # newState <- which(state.mapping$S == newS & state.mapping$I == newI)
      
      if (newState > 0 && newState <= maxNoStates &&
          newI >= 0 && newI <= params["N"] &&
          newS >= 0 && newS <= params["N"]) {
        ## if here, means that the new state is valid
        retStack[2,] <- c(newState, j, 1 - prob)
      }
    } else {
      ## we are in an absorbing state
      retStack[1,] <- c(j, j, 1)
    }
    
    return(data.frame(retStack))
  }

  .cl <- makeCluster(spec = detectCores(), type = "FORK")
  result <- parApply(cl = .cl, X = state.mapping, MARGIN = 1, FUN = f )
  stopCluster(cl = .cl)

  
  res <- rbindlist(l = result)
  res <- res[!is.na(res$newState),]
  
  # note we construct a sparse matrix, since most of the entries will be zero
  transitionMat <- Matrix(data = 0, nrow = nrow(state.mapping), ncol = nrow(state.mapping), sparse = T)
  transitionMat[cbind(res$newState, res$oldState)] <- res$prob
  
  ## identify and remove impossible states
  trueState <- (state.mapping$S + state.mapping$I <= params["N"])
  state.mapping <- state.mapping[trueState, ]
  transitionMat <- transitionMat[trueState, trueState]
  
  return(list(stateMapping = state.mapping, transitionMatrix = transitionMat ))
}