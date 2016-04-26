#***************** projection pursuit regression ***********************#
#******** Date: April 15th 2016  |  Author: Chen Daijun ****************#
#******** Project: Projection Pursuit Regression in DACE ***************#
#******** Part one: Basic Projection pursuit regression  ***************#
#******** Version : 0.0.0 **********************************************#

#************** sub-function: build non-parametric function ************#
#' non-parametric function
#' #'@include non_parametric_function.R
#------------------- Method_1: Local_linear model ----------------------#
#-----------------------------------------------------------------------#
#'local linear regression model
#'
#'@export
S <- function(projection_dir, residuals, design_matrix){
  library("np")
  res <- residuals
  projected <-  as.matrix(design_matrix %*% as.matrix(projection_dir))
  projected <- as.vector(projected)
  model.np_ll <- npreg(res ~ projected, regtype = "ll", bwmethod = "cv.aic",
                       gradients = TRUE)
  return(model.np_ll)
}
#-----------------------------------------------------------------------#

#---------- Method_2: fixed basis function: sigmoid function -----------#
#-----------------------------------------------------------------------#
#' signomid function
#'
#' @export
sigmoid <- function(z){
  sig_moid <- 1/(1 + exp(-z))
  return(sig_moid)
}
#-----------------------------------------------------------------------#

#------------------- Method_3: Smoothing Spline ------------------------#
#-----------------------------------------------------------------------#

#-----------------------------------------------------------------------#







