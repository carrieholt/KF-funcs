# Functions to run the TMB Kalman Filter  and recursive bayes Random Walk model
#===============================================
# Adapted from code from Brooke Davis and Carrie Holt
#Structure of the code copied from the sdmTMB package :
#https://github.com/pbs-assess/sdmTMB



#' Fit the analytical solution for the Kalman Filter Ricker stock recruitment model
#'
#' @param data A list or data frame containing Spawners (S) and Recruits (R) time series. No NA allowed at this point in time
#' 
#' @param silent Silent or optimization details?
#' 
#' 
#' @export
kfTMB <- function(data,  silent = FALSE, control = TMBcontrol()) {

  #===================================
  #prepare TMB input and options
  #===================================
  tmb_data <- list(
    x = data$S,
    y = log(data$R/data$S)
  )


  tmb_params <- list(
    initmeana   = lm(y~x, data=tmb_data)$coefficients[[1]],
    b           = lm(y~x, data=tmb_data)$coefficients[[2]],
    logsige     = log(1),
    logsigw     = log(1)
  )

  #to be implemented
  tmb_map <- list()
  tmb_random <- NULL

  #===================================
  # TMB fit
  #===================================

  tmb_obj <- TMB::MakeADFun(
    data = tmb_data, parameters = tmb_params, map = tmb_map,
    random = tmb_random, DLL = "Rickerkf", silent = silent)
  
  tmb_opt <- stats::nlminb(
    start = tmb_obj$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
    control = control)
  
  sd_report <- TMB::sdreport(tmb_obj)
  conv <- get_convergence_diagnostics(sd_report)

  structure(list(
    model      = tmb_opt,
    data       = data,
    tmb_data   = tmb_data,
    tmb_params = tmb_params,
    tmb_map    = tmb_map,
    tmb_random = tmb_random,
    tmb_obj    = tmb_obj,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig,
    call       = match.call(expand.dots = TRUE),
    sd_report  = sd_report),
    class      = "kfTMB")

}




#' Sample mcmc from kalman filter model uses tmbstan this function does not work well
#'
#' @param tmb_obj object from the output of the kfTMB function
#' @param iter number of iteration for mcmc sampling
#' @param chains number of mcmc chains
#' @param init Initial values for the sampler. Behaves like tmbstan
#' @param lower Vector of lower parameter bounds. Behaves like tmbstan
#' @param upper Vector of upper parameter bounds. Behaves like tmbstan
#' @param control Behaves like rstan
#' @param warmup Behaves like rstan
#' @param thin Behaves like rstan
#' 
#' @param silent Silent or optimization details?
#' 
#' @importFrom tmbstan tmbstan
#' @export
kfTMBmcmc <- function(data,  silent = FALSE, control = TMBcontrol(),
  iter=2000, chains=4, init="random", 
  lower=c(-Inf,-Inf,-Inf,-Inf), 
  upper=c(Inf,Inf,Inf,Inf),
  controlstan=list(adapt_delta = 0.95,max_treedepth = 15),
  warmup=500,
  thin = 1) {
  

  tmb_data <- list(
    x = data$S,
    y = log(data$R/data$S)
  )


  tmb_params <- list(
    initmeana   = lm(y~x, data=tmb_data)$coefficients[[1]],
    b           = lm(y~x, data=tmb_data)$coefficients[[2]],
    logsige     = log(1),
    logsigw     = log(1)
  )

  #to be implemented
  tmb_map <- list()
  tmb_random <- NULL

  
   

  #===================================
  # TMB fit
  #===================================

  tmb_obj <- TMB::MakeADFun(
    data = tmb_data, parameters = tmb_params, map = tmb_map,
    random = tmb_random, DLL = "Rickerkf", silent = silent)
  
  tmb_opt <- stats::nlminb(
    start = tmb_obj$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
    control = control)
  

  #===================================
  #prepare TMB input and options
  #===================================
   fitmcmc <- tmbstan(tmb_obj , chains = chains, iter = iter,
                init=init,
                lower=lower,
                upper=upper,
                control = controlstan, warmup = warmup,  thin = thin)
    
  mc <- extract(fitmcmc, pars=names(tmb_obj$par),
                inc_warmup=TRUE, permuted=FALSE)
  fit_summary <- summary(fitmcmc)   
  

  #===================================
  # TMB fit
  #===================================

  structure(list(
    fitmcmc    = fitmcmc,
    mcpars     = mc,
    fit_summary = fit_summary
    ))

}






#' Optimization control options stolen from sdmTMB
#'
#' Any arguments to pass to [stats::nlminb()].
#'
#' @param eval.max Maximum number of evaluations of the objective function
#'   allowed.
#' @param iter.max Maximum number of iterations allowed.
#' @param ... Anything else. See the 'Control parameters' section of
#'   [stats::nlminb()].
#'
#' @export
TMBcontrol <- function(eval.max = 1e4, iter.max = 1e4, ...) {
  list(eval.max = eval.max, iter.max = iter.max, ...)
}



#' get TMB convergence diagnostics stolen from sdmTMB
#'
#' @param sd_report A TMB sd report object
#' @export
get_convergence_diagnostics <- function(sd_report) {
  final_grads <- sd_report$gradient.fixed
  bad_eig <- FALSE
  if (!is.null(sd_report$pdHess)) {
    if (!sd_report$pdHess) {
      warning("The model may not have converged: ",
        "non-positive-definite Hessian matrix.", call. = FALSE)
    } else {
      eigval <- try(1 / eigen(sd_report$cov.fixed)$values, silent = TRUE)
      if (is(eigval, "try-error") || (min(eigval) < .Machine$double.eps * 10)) {
        warning("The model may not have converged: ",
          "extreme or very small eigen values detected.", call. = FALSE)
        bad_eig <- TRUE
      }
      if (any(final_grads > 0.01))
        warning("The model may not have converged. ",
          "Maximum final gradient: ", max(final_grads), ".", call. = FALSE)
    }
  }
  invisible(list(final_grads = final_grads, bad_eig = bad_eig))
}






#==============================================================
#recursive Bayes version


#' Fit the analytical solution for the Kalman Filter Ricker stock recruitment model
#'
#' @param data A list or data frame containing Spawners (S) and Recruits (R) time series. 
#' No NA allowed at this point in time
#' @param priorratiovar parameters for beta prior on variance allocation
#' 
#' @param silent Logical. Silent or optimization details?
#' 
#' 
#' @export
rbTMB <- function(data, priorratiovar=c(2,2), silent = FALSE, control = TMBcontrol()) {

  #===================================
  #prepare TMB input and options
  #===================================
  tmb_data <- list(
    obs_logRS = log(data$R/data$S),
    obs_S = data$S,
    prbeta1= priorratiovar[1],
    prbeta2= priorratiovar[2]
  )

  srlm <-lm(obs_logRS~obs_S, data=tmb_data)

  tmb_params <- list(
    alphao   = srlm$coefficients[[1]],
    logSmax = log(1/-srlm$coefficients[[2]]),
    rho = 0.5,
    logvarphi = log(1),
    alpha = rep(srlm$coefficients[[1]],length(tmb_data$obs_logRS))
  )

  #to be implemented
  tmb_map <- list()
  tmb_random <- "alpha"

  #===================================
  # TMB fit
  #===================================

  tmb_obj <- TMB::MakeADFun(
    data = tmb_data, parameters = tmb_params, map = tmb_map,
    random = tmb_random, DLL = "Ricker_rb_ratiovar", silent = silent)
  
  tmb_opt <- stats::nlminb(
    start = tmb_obj$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
    control = control)
  
  sd_report <- TMB::sdreport(tmb_obj)
  conv <- get_convergence_diagnostics(sd_report)

  structure(list(
    model      = tmb_opt,
    data       = data,
    tmb_data   = tmb_data,
    tmb_params = tmb_params,
    tmb_map    = tmb_map,
    tmb_random = tmb_random,
    tmb_obj    = tmb_obj,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig,
    call       = match.call(expand.dots = TRUE),
    sd_report  = sd_report),
    class      = "rbTMB")

}










#==============================================================
#Dummy function

#' Roxygen commands
#'
#' This is a dummy function who's purpose is to hold the useDynLib roxygen tag.
#' This tag will populate the namespace with compiled c++ functions upon package install.
#'
#' @useDynLib Ricker_rb_ratiovar
#' @useDynLib Rickerkf
#'
dummy <- function(){
  return(NULL)
}

# END
#***********************************************************************************

