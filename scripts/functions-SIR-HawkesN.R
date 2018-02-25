# this is a library with the methods for the SIR and the HawkesN models
library(nloptr)

##################### SIMULATING A STOCHASTIC SIR MODEL ######################

#' Generates the stochastic Continuous Time Markov Chain (discreete states) SIR
#' model. In the stochastic version, time is not given in moments at which the
#' populations are to be computed, it evolves by itself. Tmax gives the end time
#' at which simulation halts.
generate.stochastic.sir <- function(params = c(N = 1300, I.0 = 300, gamma = 0.2, beta = 1),  Tmax = 10) {
  params["b"] <- 0 ## b is the birth = death rate. In a fixed population it is fixed to zero, but I'm puting it here for later dev.
  
  ## start at time 0
  state <- data.frame(t = 0, I = params["I.0"], S = params["N"] - params["I.0"], C = params["I.0"])
  rownames(state) <- NULL
  j <- 1
  while (state$I[j]>0 & state$t[j]<Tmax) {
    ## first, compute next state rates
    ps1 <- params["beta"] * state$I[j] * state$S[j] / params["N"] ## new infection
    ps2 <- params["gamma"] * state$I[j] ## new recovery
    ps3 <- params["b"] * state$I[j] ## death in infected AND birth from infected
    ps4 <- params["b"] * ( params["N"] - state$I[j] - state$S[j]) ## birth from recovered
    stay <- 1 - ps1 - ps2 - ps3 - ps4 ## staying in the same state -> will determine time of next event
    
    ## get time of next event
    u1 <- runif(n = 1, min = 0, max = 1) # uniform random number
    a <-  ps1 + ps2 + ps3 + ps4
    newTime <- state$t[j] - log(u1) / a
    
    ## draw the next state. Each of the next states has a chance proportional to its ps1, 2, 3, 4.
    ## they need to amount to 1 (since there is no chance of staying in the current state -- accounted by the next time here above)
    u2 <- runif(n = 1, min = 0, max = 1)  # uniform random number
    
    if (u2 <= ps1/a ) { ## new infection -> S-1, I+1
      nextState <- data.frame(t = newTime, I = state$I[j]+1, S = state$S[j]-1, C = state$C[j]+1)
    } else if (u2 <= (ps1 + ps2)/a) {  ## new recovery -> S, I-1
      nextState <- data.frame(t = newTime, I = state$I[j]-1, S = state$S[j], C = state$C[j])
    } else if (u2 <- (ps1 + ps2 + ps3)/a) {  ## death in infected AND birth from infected -> S+1, I-1
      nextState <- data.frame(t = newTime, I = state$I[j]-1, S = state$S[j]+1, C = state$C[j])
    } else { ## birth from recovered -> S+1, I
      nextState <- data.frame(t = newTime, I = state$I[j], S = state$S[j]+1, C = state$C[j])
    }
    state <- rbind(state, nextState)
    j <- j + 1
    cat(sprintf("\rCurrent simulation time: %.3f / %.3f (S=%d, I=%d, R=%d, C=%d).", state$t[j], Tmax, state$S[j], state$I[j], params["N"]-state$S[j]-state$I[j], state$C[j]))
  }
  
  rownames(state) <- NULL
  state$R <- state$C - state$I
  names(state) <- c("time", "I", "S", "C", "R")
  state <- state[c("time", "S", "I", "R", "C")]
  
  cat(sprintf("\n--> Simulation done!\n"))
  return(state)
}

#' From a stochastic state data.frame (generated using generate.stochastic.sir),
#' this function extracts the event series. For compatibility with Hawkes, it
#' also creates dummy event magnitudes of 1.1.
SIR2HAWKES.stochastic.event.series <- function(state) {
  state <- as.data.frame(state)
  ## generate as many events with time 0 as there are I.0 (infections at time 0)
  history <- as.data.frame(matrix(data = rep(x = c(1.1, 0), times = state$I[1]), nrow = state$I[1], ncol = 2, byrow = T))
  colnames(history) <- c("magnitude", "time")
  
  ## compute where we had a new infection event
  state$new <- F
  state$new[-1] <- state$C[-1] > state$C[-nrow(state)]
  evnt <- data.frame(magnitude = 1.1, time = state$time[state$new])
  
  ## concatenate the two
  history <- rbind(history, evnt)
  return(history)
}

############################# FITTING A STOCHASTIC SIR MODEL ################################


# fit SIR using loglikelhood maximization
fit.stochastic.sir <- function(mysim, params.start) {
  if (params.start["N"] <  max(mysim[, -1])) {
    params.start["N"] <-  max(mysim[, -1])
  }
  
  params.start["I.0"] <- mysim$I[1]
  res <- lbfgs(x0 = unlist(params.start),                            
               fn = stochastic.sir.complete.neg.log.likelihood,          
               lower =  c(N = max(mysim[, -1]), I.0 = mysim$I[1], gamma = 0, beta = 0),
               upper = c(N = Inf, I.0 = mysim$I[1], gamma = Inf, beta = Inf),
               state = mysim)
  names(res$par) <- names(params.start)
  
  return(res)
}

#' Complete log-likelihood function for a stochastic SIR model, computed for a
#' set of params. It assumes that the state matrix is in the SIR format (with
#' "time", "I", "S", "C", "R" headers) and it contains all events in the
#' process: new infection, new recovery and (optionally) birth and death. Such a
#' matrix can be obtained using the "generate.stochastic.sir" function.
stochastic.sir.complete.neg.log.likelihood <- function(params = c(N = 1300, I.0 = 300, gamma = 0.2, beta = 1), state) {
  state <- as.data.frame(state)
  names(params) <- c("N", "I.0", "gamma", "beta")
  params <- unlist(params)
  params["b"] <- 0 ## b is the birth = death rate. In a fixed population it is fixed to zero, but I'm puting it here for later dev.
  
  ## we will sum here the log-likelihood
  ll <- 0
  
  contribs <- c(0, 0, 0, 0, 0); last_infection_time <- 0; last_recovery_time <- 0
  contrib_Lambda_inf <- 0; contrib_Lambda_rec <- 0
  ## go through events (rows) one by one. Note that the first row is the initial
  ## state (cannot compute any probs there)
  for (j in 2:nrow(state)) {
    ## first, compute next state probabilities/rates
    ps1 <- params["beta"] * state$I[j-1] * state$S[j-1] / params["N"] ## new infection
    ps2 <- params["gamma"] * state$I[j-1] ## new recovery
    ps3 <- params["b"] * state$I[j-1] ## death in infected AND birth from infected
    ps4 <- params["b"] * ( params["N"] - state$I[j-1] - state$S[j-1]) ## birth from recovered
    a <-  ps1 + ps2 + ps3 + ps4 ## the inter-event times are exponentially distributed with param a. 
    
    ## the likelihood of having waited as much as we have
    ll <-  sum(ll, dexp(x = (state$time[j] - state$time[j-1]), rate = a, log = T), na.rm = T)
    
    new_event_prob_timing <- a * exp(-(state$time[j] - state$time[j-1])*a )
    contribs[3] <- sum(contribs[3], log(new_event_prob_timing), na.rm = T)
    
    ## compute contributions to Lambda by current component
    contrib_Lambda_inf <- sum(contrib_Lambda_inf, (state$time[j] - state$time[j-1]) * ps1)
    contrib_Lambda_rec <- sum(contrib_Lambda_rec, (state$time[j] - state$time[j-1]) * ps2)
    
    ## and now, how likely to have observed an event such as what we've seen
    if (state$I[j] == state$I[j-1]+1 && state$S[j] == state$S[j-1]-1) { ## new infection -- S-1, I+1
      ll <- sum(ll, log(ps1/a), na.rm = T)
      
      contribs[1] <- sum(contribs[1], log(ps1/a), na.rm = T)
      
      new_inf_prob_timing <- ps1 * exp(-contrib_Lambda_inf)
      contribs[4] <- sum(contribs[4], log(new_inf_prob_timing), na.rm = T)
      contrib_Lambda_inf <- 0
    } else if (state$I[j] == state$I[j-1]-1 && state$S[j] == state$S[j-1]) {  ## new recovery -- S, I-1
      ll <- sum(ll, log(ps2/a), na.rm = T)
      contribs[2] <- sum(contribs[2], log(ps2/a), na.rm = T)

      
      new_rec_prob_timing <- ps2 * exp(-contrib_Lambda_rec)
      contribs[5] <- sum(contribs[5], log(new_rec_prob_timing), na.rm = T)
      contrib_Lambda_rec <- 0
    } else if (state$I[j] == state$I[j-1]-1 && state$S[j] == state$S[j-1]+1) {## death in infected AND birth from infected -- S+1, I-1
      ll <- sum(ll, log(ps3/a), na.rm = T)
    } else if (state$I[j] == state$I[j-1] && state$S[j] == state$S[j-1]+1) { ## birth from recovered -- S+1, I
      ll <- sum(ll, log(ps4/a), na.rm = T)
    }
  }
  
  names(ll) <- "neg.ll"
  return(nll = -ll)
}


########################### FITTING HawkesN #############################

fitSeries <- function(history, params_init) {
  model <- list(par = c(K = NA, c = NA, theta = NA, N = NA), value = NA)
  tryCatch({
    # ## optimize with NLOPT
    model <- find.optimal.parameters(history = history, init_params = params_init)
  }, error = function(err) {
    print(paste("[fitSeries] Error in optim:  ", err))
    model <- list(par = c(K = NA, c = NA, theta = NA, N = NA), value = NA)
  })

  return(model)
}

# # function to generate random initial points
# generateRandomPoints <- function(seed=NULL){
#   init_k <- runif(5, min = .Machine$double.eps, max = 10)
#   init_c <- runif(5, min = .Machine$double.eps, max = 300)
#   init_theta <- runif(5, min = .Machine$double.eps, max = 3)
#   init_n <- floor(runif(5, min = 50, max = 5000))
#   ## 3 known good start points, all NA (which invokes the global optimizer) and all Inf (which lets IPOPT chose its own start point)
#   init_k[6:10] <- c(0.1, 0.5, 0.8, NA, Inf)
#   init_beta[6:10] <- c(0.001, 0.5, 0.8, NA, Inf)
#   init_c[6:10] <- c(10, 100, 60, NA, Inf)
#   init_theta[6:10] <- c(0.0001, 0.00001, 0.000001, NA, Inf)
#   init_n[6:10] <- c(60, 200,600, NA, Inf)
#   init <- data.frame('K' = init_k, 'beta' = init_beta,
#                      'c' = init_c, 'theta' = init_theta, 'N' = init_n)
#   return(init)
# }

find.optimal.parameters <- function(history, init_params = c(K = 0.1, c = 0.1, theta = 0.1, N = 1000), 
                                    lowerBound = c(K = 0, c = .Machine$double.eps, theta = 0, N = 1), 
                                    upperBound = c(K = Inf, c = 300, theta = Inf, N = Inf),
                                    iterations = 5000, ...) {

  if ( sum(! names(init_params) %in% names(lowerBound)) > 0 ) {
    newParamName <- names(init_params)[! names(init_params) %in% names(lowerBound)]
    warning(sprintf("You gave me the param '%s', but no boundry for it. I'll assume entire real range.", newParamName ))
    lowerBound[newParamName] <- -Inf
    upperBound[newParamName] <- Inf
  }
  lowerBound <- lowerBound[names(init_params)]
  upperBound <- upperBound[names(init_params)]
  
  ## correct initial parameters out of bounds -- should not happen, but...
  init_params[init_params > upperBound] <- upperBound[init_params > upperBound]
  init_params[init_params < lowerBound] <- lowerBound[init_params < lowerBound]
  
  ## if my initial N (population size) is below what I observe, I know it is a
  ## bad starting point. Therefore, I will correct it to start at least where I
  ## see it.
  if ("N" %in% names(init_params) && init_params["N"] < nrow(history))
    init_params["N"] <- nrow(history)
  
  res <- lbfgs(x0 = unlist(init_params),
               fn = neg.log.likelihood,
               gr = closedGradient,
               lower = lowerBound,
               upper = upperBound,
               control = list(maxeval = iterations, xtol_rel = "1e-8"),
               history = history,  ...)
  
  names(res$par) <- names(init_params)
  res$value <- neg.log.likelihood(params = res$par, history = history, ... )
  
  return(res)
}

#' Next function calculates the negative
#' log-likelihood. This is given that the optim function minimizes a function. '
#' Minimizing the -1 * log-likelihood amounts to maximizing the log-likehood. '
#' @param .cl - ugly hack to pass the cluster parameter to the gradient function
#'  in the optimization procedure. Not necessary for this function. '
neg.log.likelihood <- function(params, history) { 
  names(params) <- c("K", "c", "theta", "N")
  params <- as.list(unlist(params))
  
  ## if all good, continue
  bigT <- max(history$time)
  ## startEvent is the first event generated endogeneously (not in the initial pool of events)
  startEvent <- which.max(history$time > 0)
  endEvent <- nrow(history)
  
  ## get the formula of the log-likelihood * -1.
  ## that in the summation, we eliminated the first event, as there is no background rate in our lambda(t)
  lambdaInt <- integrateLambda(lower = 0, upper = bigT, history = history, params = params)
  lambdaSums <- sum(log(lambda(t = history$time[startEvent:endEvent], history = history, params = params)))
  retVal <- lambdaInt - lambdaSums 
  
  if ("N" %in% names(params)) {
    ## apply ridge regularizer (with alpha = 10, desired N val is the observed
    ## size of the cascade) in the HawkesN case
    retVal <- retVal
  }
  
  ## some optimizers don't like Inf, NA or NaN. Set it to the maximum real number
  if (!is.finite(retVal)) {
    warning(paste("Following params got an error value (", retVal, "): ", paste(params, collapse = ", "), sep = ""))
    retVal <- NA
  }
  
  return( retVal)
}

#' A function to calculate the value of the integral of lambda for 
#' Exponential Kernel with finite population
integrateLambda <- function(lower, upper, history, params, mmin = 1) {
  names(params) <- c("K", "c", "theta", "N")
  params <- as.list(unlist(params))
  endEvent <- nrow(history)
  
  ## closed form
  res <- sum(sapply(X = 1:nrow(history), FUN = function(i) {
    j <- i:endEvent
    ## make sure that the counter of available events never goes below zero,
    ## otherwise we get a negative integral
    avail_events <- params$N - i:(endEvent-1)
    avail_events[avail_events < 0] <- NA
    
    terms <- exp(-params$theta * (history$time[j] - history$time[i]) )
    
    val <- sum((avail_events / params$N) * (terms[-length(terms)] - terms[-1]))
    return(val)
  })) * params$K
  
  return(res)
}

#' Computes the conditional intensity of the non-homogenous Poisson process (
#' lambda(t) ) It is the equivalent of the CIF function in "simulation-model.R",
#' but updated for passing parameters as a list and calculating the intensity
#' for a vector on input (needed for numerical integration).
lambda <-  function(t, history, params = c(K = 0.024, c = 0.001, theta = 0.2, N = 1000)) {
  params <- as.list(unlist(params))
  
  res <- sapply(X = t, FUN = function(ti) {
    subst <- history[history$time <= ti,]
    return(sum(kernelFct(event = subst, t = ti, params = params)))
  })
  
  res[res == 0] <- .Machine$double.eps
  
  return(res)
}

#' calculates the influence of a given event at the givent time the event is a 2
#' elements list (mi, ti). Inclusive parameter means that if an event is present
#' at time t, then it contributes to the conditional intensity. If inclusive ==
#' F, then events at time t are removed.
kernelFct <- function(event, t, params = c(K = 0.024, theta = 0.2, N = 1000), alpha = 2.016, mmin = 1) {
  params <- unlist(params)
  # the event has 2 components: (magnitude_i, time_i)
  mat_event = matrix(unlist(event), ncol = 2, byrow = F)
  mi = mat_event[,1]
  ti = mat_event[,2]
  Nt <- min(sum(ti <= t), params["N"])

  
  # f(p_j) part - virality of a diffusion. Constant for a given diffusion. Furthermore, discount for available events.
  fun_f <- params["K"] * (1 - Nt / params["N"])
  
  fun_psi <- params["theta"] * (exp(-params["theta"] * (t - ti)))
  
  val = fun_f * fun_psi
  val[t<ti] = 0
  val[mi<mmin] = 0
  
  return(val)
}

#' A function to calculate the dervative in closed form for Finite Exponential Kernel.
closedGradient <- function(params, history){
  names(params) <- c("K", "c", "theta", "N")
  
  ## the wrapping insures that we always return as many values as parameters
  ## compose first the return vector with zero everywhere
  ret <- rep(x = 0, times = length(params))
  names(ret) <- names(params)

  # variables collection
  params <- as.list(unlist(params))
  bigT <- max(history$time)
  n <- nrow(history)
  k <- params$K
  theta <- params$theta
  N <- params$N
  
  ## number of initial events 
  i_0 <- nrow(history[history$time <= 0,])
  
  # In all the following calcuaion, res is the part of derivative coming 
  # from bigLambda part and first is derivative coming from 
  # summation of log(lambda(ti))
  
  #calculating derivative wrt k in closed form
  res <- sum(sapply(X = 1:n, FUN = function(i) {
    ## compute all the exponentials e^{-\theta (t_j - t_i)}
    j <- i:n
    terms <- exp(-params$theta * (history$time[j] - history$time[i]) )
    
    val <- sum(((N - i:(n-1)) / N) * (terms[-length(terms)] - terms[-1]))
    return(val)
  }))
  derivK <- ((n-i_0) / k) - res
  
  #calculating derivative wrt theta in closed form
  res <- sum(sapply(X = 1:n, FUN = function(i) {
    ## compute all the exponentials e^{-\theta (t_j - t_i)}
    j <- i:n
    terms <- exp(-params$theta * (history$time[j] - history$time[i]) )
    times <- history$time[j] - history$time[i]
    
    val <- sum(((N - i:(n-1)) / N) * ((times[-1] * terms[-1]) - (times[-length(times)] * terms[-length(terms)])))
    return(val)
  })) * k
  
  first <- sum(sapply(X = (i_0+1):(n), FUN = function(i) {
    
    numerator <- sum(exp(-theta * (history$time[i] - history$time[1:i-1])) * 
                       - (history$time[i] - history$time[1:i-1]))
    
    denominator <- sum(exp(-theta * (history$time[i] - history$time[1:i-1])))
    return(numerator / denominator)
  }))
  derivTheta <- ((n-i_0) / theta) + first - res
  
  #calculating derivative wrt N in closed form
  res <- sum(sapply(X = 1:n, FUN = function(i) {
    ## compute all the exponentials e^{-\theta (t_j - t_i)}
    j <- i:n
    terms <- exp(-params$theta * (history$time[j] - history$time[i]) )
    
    val <- sum( ( ( i:(n-1) ) / N^2 ) * ( terms[-length(terms)] - terms[-1] ) )
    return(val)
  })) * k
  
  ## n! = gamma(n+1), therefore deriv of log(n!) is deriv lof log(gamma(n+1)) which is digamma(n+1)
  # first <- digamma(N-i_0+1) - ((n-i_0)/N)
  first <- sum(sapply(X = (i_0+1):n, FUN = function(i) {return(1/(N-i+1))})) - ((n-i_0)/N)
  derivN <- first - res
  
  
  
  ## added a -1 because we are minimizing log-likelihood
  retval <- c(K = derivK, theta = derivTheta, N = derivN)
  
  ## put the vals we got into the return array
  ret[names(retval)] <- -1 * retval
  return(ret)
}