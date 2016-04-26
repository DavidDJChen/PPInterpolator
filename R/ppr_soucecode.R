#***************** projection pursuit regression ***********************#
#******** Date: April 15th 2016  |  Author: Chen Daijun ****************#
#******** Project: Projection Pursuit Regression in DACE ***************#
#******** Part one: Basic Projection pursuit regression  ***************#
#******** Version : 0.0.0 **********************************************#


#**************** step 1: centerize the response ***********************#

centerize <- function(response, design_matrix){
  Y <- response   # initial response
  X <- design_matrix #
  all_one_mat <- matrix(rep(1, (nrow(X))^2), nrow = nrow(X))
  X_centered <- X - 1/nrow(X)* all_one_mat %*% X
  cent_outcome <- list(cent_resp = Y, cent_design = X_centered)
  return(cent_outcome)  # centerizing formula: Y - mean(Y)
} # return a vector

#***********************************************************************#

#**************** step 2: initializaion ********************************#

initial <- function(centered_resp){ # input centered response
  res <- centered_resp # initialize residuals
  M <- 0               # initialize counter
  initial_outcome <- list(residuals = res, counter = M)# output as a list
  return(initial_outcome)
}

#***********************************************************************#

#**************** step 3: calculate the figure_of_merit ****************#
#************************** criterion of fit ***************************#

figure_of_merit <- function(residuals, fit_values_alpha){
  fig_of_mer <- 1 - sum((r - fit_values_alpha)^2) / (sum( r^2 ))
  return(fig_of_mer) # inputs: residuals and fitted values
}

#***********************************************************************#

#**************** step 4: univariate smoothing *************************#

#****************** Method 1: fixed S(x) *******************************#
#******************** Using sigmodi function ***************************#

sigmoid <- function(z){
  sig_moid <- 1/(1 + exp(-z))
  return(sig_moid)
}

#***********************************************************************#

#***************** Method 2: local linear regression *******************#
#********************* Using packages: np + cubature *******************#

# calculate the local linear response
#-----------------------------------------------------------------------#
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

#----------- calculate the derivate of the local linear ----------------#
#--------------------- Second order Gaussian kernel---------------------#
Gaussian_2 <- function(center, scatter_point, brand_width){
  x <- center
  x_i <- scatter_point
  h <- brand_width
  w <- 1 / sqrt(2 * pi) * exp(- 0.5 * ((x - as.vector(x_i)) / h) ^ 2)
  return(w)
}
#-----------------------------------------------------------------------#

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

#---------------------- S_1() function: part of derivate ----------------#
S_1 <- function(center, sample_points, brand_width){
  h <- brand_width
  x <- center
  x_i <- sample_points
  s1_x <- sum(1 / h * Gaussian_2(x, x_i, h) * (x - as.vector(x_i)))
  return(s1_x)
}
#-----------------------------------------------------------------------#

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

#---------------------- S_2() function: part of derivate ---------------#
S_2 <- function(center, sample_points, brand_width){
  h <- brand_width
  x <- as.vector(center)
  x_i <- as.vector(sample_points)
  s2_x <- sum(1 / h * Gaussian_2(x, x_i, h) * (x - x_i) ^ 2)
  return(s2_x)
}
#-----------------------------------------------------------------------#

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
inf_norm <- function(vec_1, vec_2){
  max_elem <- max(abs(vec_1 - vec_2))
  return(max_elem)
}
#----------------------------------------------------------------------#

#---- projection_dir_final function: converge the projection_dir ------#
projection_dir_final <- function(residuals, design_matrix,
                                 initial_projected_dir, epslon){
  res <- residuals
  old_dir <- initial_projected_dir
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

#*** step 7: increasing terms to fit data from Computer Experiments ***#

#---- PPR_CE: build a final PPR model for computer experiments --------#
#----------------------------------------------------------------------#
PPR_CE <- function(response, design_matrix, max_terms, method = "ll",
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





