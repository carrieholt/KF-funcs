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
kfTMB <- function(data,  silent = FALSE, control = kfTMBcontrol()) {

  #===================================
  #prepare TMB input and options
  #===================================
  tmb_data <- list(
    x = data$S,
    y = log(data$R/data$S)
  )


  tmb_params <- list(
    initmeana   = lm(y~x, data=tmb_data)$coefficients[[1]],
    loginitvara = log(1),
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
kfTMBcontrol <- function(eval.max = 1e4, iter.max = 1e4, ...) {
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


# END
#***********************************************************************************

