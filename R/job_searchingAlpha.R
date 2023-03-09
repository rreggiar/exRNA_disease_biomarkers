library(tidyverse)
library(glmnet)
library(pROC)

doMC::registerDoMC(cores =12)


lapply(seq(0,1,0.1), function(alpha) {
  set.seed(2022-03-07)
  
  print(alpha)

  mult_model <- 
    cv.glmnet(x = elife_x_in_train,
              y = elife_y_in_train,
              family = 'multinomial',
              type.measure = 'class',
              nfolds = 20, 
              alpha = alpha, 
              parallel = TRUE,
              keep = TRUE)
  
  
  coef_tmp <- Reduce(cbind, coef(mult_model, s = mult_model$lambda.1se))
  
  beta <- coef_tmp[apply(coef_tmp != 0, 1, any), ]
  
  colnames(beta) <- names(coef(mult_model, s = mult_model$lambda.1se))
  
  beta
  
  inverse.logit.calc <- function(y_hat) { exp(y_hat)/(1+exp(y_hat)) }
  softmax.logit.calc <- function(y_hat) { exp(y_hat)/sum(exp(y_hat)) }
  
  y_hat_matrix <- mult_model$fit.preval[,,which(mult_model$lambda == mult_model$lambda.min)]
  
  # multi-nomial so need to apply softmax to each row (sample)
  mult_response_matrix <- 
    t(apply(y_hat_matrix, 1, softmax.logit.calc))
  
  mult_response_matrix
  
  # which column is the max -- so need to pull from column names
  y_hat_predictions <- 
    colnames(mult_response_matrix)[apply(mult_response_matrix, 1, which.max)]
  
  multi_roc <- 
    pROC::multiclass.roc(response = elife_y_in_train, 
                         predictor = mult_response_matrix, 
                         levels = c('0','1','2', '3', '4', '5'))
  
  
  caret::confusionMatrix(data = as.factor(y_hat_predictions),
                         reference = as.factor(elife_y_in_train)) -> cm_out
  
  ret_lst <- list(cm_out)
  names(ret_lst) <- alpha
  
  ret_lst
  
}) %>% purrr::flatten() -> alpha_list

return(alpha_list)
