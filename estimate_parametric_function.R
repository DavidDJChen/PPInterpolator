#***************** projection pursuit regression ***********************#
#******** Date: April 15th 2016  |  Author: Chen Daijun ****************#
#******** Project: Projection Pursuit Regression in DACE ***************#
#******** Part one: Basic Projection pursuit regression  ***************#
#******** Version : 0.0.0 **********************************************#
#'Functions build for estimating a given non-parametric model
#'@include estimate_parametric_function.R
#************** sub-function: estimate the parametric parts ************#

#------------------ Using quasi-Newton method to solve -----------------#
#------------------ need responses and their derivates -----------------#

#----------- calculate the derivate of the local linear ----------------#
#--------------------- Second order Gaussian kernel---------------------#
#'2-order Gaussian kernel function
#'
#'@export
Gaussian_2 <- function(center, scatter_point, brand_width){
  x <- center
  x_i <- scatter_point
  h <- brand_width
  w <- 1 / sqrt(2 * pi) * exp(- 0.5 * ((x - as.vector(x_i)) / h) ^ 2)
  return(w)
}
#-----------------------------------------------------------------------#
#'the derivate of 2-order Gaussian kernel function
#'
#'@export
#--------------------- derivates of 2-order Gaussian kernel ------------#
derivate_Gaussian_2 <- function(center, scatter_point, brand_width){
  x <- as.vector(center)
  x_i <- scatter_point
  h <- brand_width
  der_Gaussian_2 <- 1 / sqrt(2 * pi) * exp(- 1 / 2 * (x - x_i) ^ 2 / h ^ 2) *
    (x_i - x) / h
  return(der_Gaussian_2)
}
#-----------------------------------------------------------------------#

#'S_1 function is a part of the derivate of local linear regression function
#'
#'@export
#---------------------- S_1() function: part of derivate ----------------#
S_1 <- function(center, sample_points, brand_width){
  h <- brand_width
  x <- center
  x_i <- sample_points
  s1_x <- sum(1 / h * Gaussian_2(x, x_i, h) * (x - as.vector(x_i)))
  return(s1_x)
}
#-----------------------------------------------------------------------#

#' derivate_S1 function is the derivate of S_1 function
#'
#' @export
#---------------- derivate_S1 function: derivate of S_1() --------------#
derivate_S1 <- function(center, sample_points, brand_width){
  x <- as.vector(center)
  x_i <- sample_points
  h <- brand_width
  der_S1 <- sum(1 / h * derivate_Gaussian_2(x, x_i, h) * (x - x_i) +
                  1 / h * Gaussian_2(x , x_i, h))
  return(der_S1)
}



#-----------------------------------------------------------------------#
#' S_2 funciton is another part of the derivate of the local linear regression function
#'
#'  @export
#---------------------- S_2() function: part of derivate ---------------#
S_2 <- function(center, sample_points, brand_width){
  h <- brand_width
  x <- as.vector(center)
  x_i <- as.vector(sample_points)
  s2_x <- sum(1 / h * Gaussian_2(x, x_i, h) * (x - x_i) ^ 2)
  return(s2_x)
}
#-----------------------------------------------------------------------#
#' derivate_S2 is the derivate of S_2 function
#'
#' @export
#--------------- derivate_S2 function: derivate of S_2() ---------------#
derivate_S2 <- function(center, sample_points, brand_width){
  x <- center
  x_i <- sample_points
  h <- brand_width
  der_S2 <- sum(1 / h * derivate_Gaussian_2(x, x_i, h) * (x - x_i) ^ 2 +
                  2 / h * Gaussian_2(x, x_i, h) * (x - x_i))
  return(der_S2)
}
#-----------------------------------------------------------------------#
#'function named a() is the main part of the derivate of local linear regression. g_hat'(x) = sum(a(X,h)*Y)/sum(a(X,h))
#'
#'@export
#---------------------- a: function ------------------------------------#
a <- function(center, sample_points, brand_width){
  x <- center
  x_i <- as.vector(sample_points)
  h <- brand_width
  a <- 1 / h * Gaussian_2(x, x_i, h) * (1 - S_1(x, x_i, h) / S_2(x, x_i, h) *
                                          (x - x_i))
  return(a)
}
#-----------------------------------------------------------------------#

#---------------------- derivative_a() ---------------------------------#
#' derivate_a function is the derivate of the function a
#'
#' @export
derivate_a <- function(center, sample_points, brand_width){
  h <- brand_width
  x <- as.vector(center)
  x_i <- sample_points
  der_a <- 1 / h * derivate_Gaussian_2(x, x_i, h) *
    (1 - S_1(x , x_i, h) / S_2(x, x_i, h) * (x - x_i))
  + 1 / h * Gaussian_2(x, x_i, h) * (- derivate_S1(x, x_i, h) / S_2(x, x_i, h)
                                     + S_1(x, x_i, h) / (S_2(x, x_i, h)^2)
                                     * derivate_S2(x, x_i, h) * (x - x_i)
                                     - S_1(x, x_i, h)/ S_2(x, x_i, h))
}
#-----------------------------------------------------------------------#

#---------------------- The derivate of S ------------------------------#
#' S_derivate function is the derivate of the function S()
#'
#' @export
S_derivate <- function(center, sample_points, brand_width, response){
  h <- brand_width
  x <- as.vector(center)
  x_i <- sample_points
  y <- response
  a <- a(x, x_i, h)
  der_S <- (sum(derivate_a(x, x_i, h) * y) * sum(a) - sum(a * y) *
              sum(derivate_a(x, x_i, h))) / sum(a ^ 2)
  return(der_S)
}
#-----------------------------------------------------------------------#

#***********************************************************************#

#************* Step 5: Estimating the projection direction *************#

#---------------- function : find the new projection direction ---------#
#'new_projected_dir is the function which input residuals, design_matrix,
#' old projected direction and non_parametric model to output the new
#' direction generated by the quasi-Newton algorithm.
#'
#'  @export
new_projected_dir <- function(residuals, design_matrix, old_projected_dir,
                              np_model){
  res <- residuals
  design_matrix <- design_matrix
  dir_old <- old_projected_dir
  res_pred <- predict(np_model, # the np_model is fixed
                      newdata = as.data.frame(design_matrix %*% dir_old))
  h <- np_model$bw # np_model gives the band width automatically
  res_pred_der <- matrix(rep(0, nrow(design_matrix)), ncol = 1)
  for(i in 1:nrow(design_matrix)){
    res_pred_der[i, 1] <- S_derivate(design_matrix[i, ] %*% dir_old,
                                     design_matrix %*% dir_old, h, res) #key point
  }
  weight_reg_response <- design_matrix %*% as.matrix(dir_old) +
    (res - res_pred)/res_pred_der
  weight <- (res_pred_der) ^ 2

  design_matrix <- design_matrix

  weight_lm_model <- lm(weight_reg_response ~ 0 + design_matrix, weights = weight)

  dir_new <- matrix(unname(weight_lm_model$coefficients), ncol = 1)
  return(dir_new)
}
#----------------------------------------------------------------------#

#**********************************************************************#

#************** Step 6: iterating the direction to converge ***********#

#----------- inf_norm: measure the distance between two vectors -------#
#'inf_norm function can return the infinite norm of a finite-dimensional
#'vector. That is inf_norm return the maximum value of the absolute value of
#'all vector's components
#'
#'@export
inf_norm <- function(vec_1, vec_2){
  max_elem <- max(abs(vec_1 - vec_2))
  return(max_elem)
}
#----------------------------------------------------------------------#

#---- projection_dir_final function: converge the projection_dir ------#
#'projection_dir_final function can obtain a final convergencing projection
#'direction after we input the residuals, design matrix, initial projected
#'direction and the criterion value (epslon) for a given non-parametric
#'regression model
#'
#'@export
projection_dir_final <- function(residuals, design_matrix,
                                 initial_projected_dir, epslon){
  res <- residuals
  old_dir <- initial_projected_dir
  design_matrix <- design_matrix
  epslon <- epslon
  S <- S(old_dir, res, design_matrix)
  new_dir <- new_projected_dir(res, design_matrix, old_dir, S)
  while(inf_norm(old_dir, new_dir) > epslon){
    old_dir <- new_dir
    new_dir <- new_projected_dir(res, design_matrix, old_dir, S)
  }
  return(new_dir)
}
#----------------------------------------------------------------------#

#**********************************************************************#
