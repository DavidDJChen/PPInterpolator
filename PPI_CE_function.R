#***************** projection pursuit regression ***********************#
#******** Date: April 15th 2016  |  Author: Chen Daijun ****************#
#******** Project: Projection Pursuit Regression in DACE ***************#
#******** Part one: Basic Projection pursuit regression  ***************#
#******** Version : 0.0.0 **********************************************#
#'Functions build for estimating a given non-parametric model
#'@include estimate_parametric_function.R
#************** sub-function: estimate the parametric parts ************#
#*** step 7: increasing terms to fit data from Computer Experiments ***#

#---- PPR_CE: build a final PPR model for computer experiments --------#
#----------------------------------------------------------------------#
#'PPI_CE function is the main function of PPI algorithm, which input
#'the response and design matrix, max terms of PPI model, and the non-
#'parametric method(the default method is local linear regression), the
#'convergence requirement of parameter part, and the required accuracy for
#'the model.
#'
#' @export
PPI_CE <- function(response, design_matrix, max_terms, method = "ll",
                   epslon, required_accuracy){
  Y <- response
  design_matrix <- design_matrix
  epslon <- epslon
  max_terms <- max_terms
  accuracy <- required_accuracy
  res <- Y # add the 1st term
  dir_old <- matrix(runif(ncol(design_matrix)), ncol = 1)
  direction <- projection_dir_final(res, design_matrix, dir_old, epslon)
  S <- S(dir_old, res, design_matrix)
  res_pred <- predict(S, # the np_model is fixed
                      newdata = as.data.frame(design_matrix %*% direction))
  model_list <- list()
  model_term <- 1 # the number of model terms
  model_list[[1]] <- list(model = S, direction = direction)
  while(inf_norm(res_pred, res) > accuracy ){
    res <- Y - res_pred # add the 1st term
    dir_old <- matrix(runif(ncol(design_matrix)), ncol = 1)
    direction <- projection_dir_final(res, design_matrix, dir_old, epslon)
    S <- S(dir_old, res, design_matrix)
    res_pred <- predict(S, # the np_model is fixed
                        newdata = as.data.frame(design_matrix %*% direction))
    model_term <- model_term + 1
    model_list[[model_terms]] <- list(model = S, direction = direction)
  }

  return(model_list)
}
#----------------------------------------------------------------------#

#**********************************************************************#
