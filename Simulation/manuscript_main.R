library(tidyverse)
library(MASS)
library(BB)
library(numDeriv) 
library(conflicted)
library(dplyr)


# sim_data
sim_data <- function(N, para){
  beta_Vw <- para$beta_Vw
  beta_Z1 <- para$beta_Z1
  beta_Z2 <- para$beta_Z2
  beta_Z3 <- para$beta_Z3
  beta_Z4 <- para$beta_Z4
  beta_Y1 <- para$beta_Y1
  beta_Y0 <- para$beta_Y0
  Vc1 <- rnorm(N, mean=0, 2)
  Vc <- rnorm(N, mean=0, 1)
  Vw <- Vc + rnorm(N,mean=0,1) 
  Z <- rbinom(N, size = 1, prob = plogis(cbind(1, Vc1, Vc) %*% beta_Z1)) 
  Y1 <- cbind(1, Vc1, Vc^3, Vw) %*% unlist(beta_Y1) + rnorm(N)  # Y1 = beta1 + beta2 * Vc1 + beta3 * Vc^3 + beta4 * Vw + epsilion
  Y0 <- cbind(1, Vc1, Vc^3, Vw) %*% unlist(beta_Y0) + rnorm(N)  # Y1 = alpha1 + alpha2 * Vc + alpha3 * Vc^2 + alpha4 * Vw + epsilion
  Y <- (Z * Y1 + (1 - Z) * Y0) |> as.vector()
  return(list(Vc = Vc, Vw = Vw, Z = Z, Y = Y, Y1 = Y1, Y0 = Y0, Vc1 = Vc1))
}

# parameters
N = 1000
beta_Vw <- 1
beta_Z1 <- c(-0.1, 0.4, -0.2)
beta_Z2 <- c(-0.1, 0.4, -0.2, 0)
beta_Z3 <- c(-0.1, -0.2)
beta_Z4 <- c(-0.1, -0.2, 0)
beta_Y1 <- list(alpha0 = 2, alpha_Vc1 = rep(-3, 1), alpha_Vc3 = rep(2, 1), alpha_Vw = rep(-3, 1))
beta_Y0 <- list(eta0 = 2, eta_Vc1 = rep(-3, 1), eta_Vc3 = rep(1, 1), eta_Vw = rep(-3, 1))
paras <- list(beta_Vw = beta_Vw, beta_Z1 = beta_Z1, beta_Z2 = beta_Z2, beta_Z3 = beta_Z3, beta_Z4 = beta_Z4, beta_Y0 = beta_Y0, beta_Y1 = beta_Y1)


# Q_beta  initial value
initial_para <- list(Q_true_Vcw = c(1.5, -2, 1, -2, -2, 1), 
                     Q_true_Vc = c(1.7, -2, 1, -2, 0.5),
                     Q_wrong_Vcw = rep(1, 3),
                     Q_wrong_Vc = rep(1, 2))

# e_alpha  initial value
e_par <- list(
  e_true_Vc = c(-0.1,0.4,-0.2),
  e_true_Vcw = c(-0.1,0.4,-0.2, 0),
  e_wrong_Vc =c(-0.1,0.4),
  e_wrong_Vcw =c(-0.1,0.4,0) 
)



source("～/manuscript_base_function.R")

# output
param_settings <- list(
  list(covarite_set = 'c',  outcome_model = TRUE,  type = 'ipw', e_model = TRUE,  pscore = 'known',     e_par = e_par),
  list(covarite_set = 'c',  outcome_model = TRUE,  type = 'ipw', e_model = FALSE, pscore = 'known',     e_par = e_par),
  list(covarite_set = 'c',  outcome_model = FALSE, type = 'ipw', e_model = TRUE,  pscore = 'known',     e_par = e_par),
  list(covarite_set = 'c',  outcome_model = FALSE, type = 'ipw', e_model = FALSE, pscore = 'known',     e_par = e_par),
  list(covarite_set = 'cw', outcome_model = TRUE,  type = 'ipw', e_model = TRUE,  pscore = 'known',     e_par = e_par),
  list(covarite_set = 'cw', outcome_model = TRUE,  type = 'ipw', e_model = FALSE, pscore = 'known',     e_par = e_par),
  list(covarite_set = 'cw', outcome_model = FALSE, type = 'ipw', e_model = TRUE,  pscore = 'known',     e_par = e_par),
  list(covarite_set = 'cw', outcome_model = FALSE, type = 'ipw', e_model = FALSE, pscore = 'known',     e_par = e_par),
  list(covarite_set = 'c',  outcome_model = TRUE,  type = 'ipw', e_model = TRUE,  pscore = 'not known', e_par = e_par),
  list(covarite_set = 'c',  outcome_model = TRUE,  type = 'ipw', e_model = FALSE, pscore = 'not known', e_par = e_par),
  list(covarite_set = 'c',  outcome_model = FALSE, type = 'ipw', e_model = TRUE,  pscore = 'not known', e_par = e_par),
  list(covarite_set = 'c',  outcome_model = FALSE, type = 'ipw', e_model = FALSE, pscore = 'not known', e_par = e_par),
  list(covarite_set = 'cw', outcome_model = TRUE,  type = 'ipw', e_model = TRUE,  pscore = 'not known', e_par = e_par),
  list(covarite_set = 'cw', outcome_model = TRUE,  type = 'ipw', e_model = FALSE, pscore = 'not known', e_par = e_par),
  list(covarite_set = 'cw', outcome_model = FALSE, type = 'ipw', e_model = TRUE,  pscore = 'not known', e_par = e_par),
  list(covarite_set = 'cw', outcome_model = FALSE, type = 'ipw', e_model = FALSE, pscore = 'not known', e_par = e_par)
)

# Parallel experiments
folder_path <- "~/Desktop/result_both/"
dir.create(folder_path, showWarnings = FALSE)  

library(foreach)
library(doParallel)
num_cores <- detectCores() - 1
cl <- parallel::makeCluster(num_cores)
registerDoParallel(cl)
n_sim <- 200

results <- foreach(seed_id = 1:n_sim,
                   .packages = c('tidyverse', 'MASS', 'BB', 'numDeriv', 'conflicted', 'dplyr')) %dopar% {
                     
                     set.seed(seed_id)   
                     data <- sim_data(N, paras)
                     Vc1 <- data[['Vc1']]
                     Vc  <- data[['Vc']]
                     Vw  <- data[['Vw']]
                     Z   <- data[['Z']]
                     Y   <- data[['Y']]
                     
                     
                     temp_usual <- vector('list', length = 16)
                     temp_imp   <- vector('list', length = 16)
                     
                     for (j in seq_along(param_settings)) {
                       args <- param_settings[[j]]
                       
                       # usual estimator
                       temp_usual[[j]] <- usual_est(
                         para          = initial_para,
                         covarite_set  = args$covarite_set,
                         outcome_model = args$outcome_model,
                         type          = args$type,
                         e_model       = args$e_model,
                         pscore        = args$pscore,
                         e_par         = args$e_par
                       )
                       
                       # improved estimator
                       temp_imp[[j]] <- imp_est(
                         para          = initial_para,
                         covarite_set  = args$covarite_set,
                         outcome_model = args$outcome_model,
                         type          = args$type,
                         e_model       = args$e_model,
                         pscore        = args$pscore,
                         e_par         = args$e_par
                       )
                     }
                     
                     save(temp_usual, temp_imp,
                          file = paste0(folder_path, seed_id, ".rdata"))
                     return(list(usual = temp_usual, imp = temp_imp))
                   }
parallel::stopCluster(cl)

# save
setwd(folder_path)
file_list  <- list.files(pattern = "\\.rdata$", ignore.case = TRUE)
num_files  <- length(file_list)


tau_usual <- matrix(NA, nrow = num_files, ncol = 16)
var_usual <- matrix(NA, nrow = num_files, ncol = 16)
tau_imp   <- matrix(NA, nrow = num_files, ncol = 16)
var_imp   <- matrix(NA, nrow = num_files, ncol = 16)

for (i in seq_along(file_list)) {
  load(file_list[i])
  tau_usual[i, ] <- sapply(temp_usual, function(x) x$tau)
  var_usual[i, ] <- sapply(temp_usual, function(x) x$Var)
  tau_imp[i, ]   <- sapply(temp_imp,   function(x) x$tau)
  var_imp[i, ]   <- sapply(temp_imp,   function(x) x$Var)
}

write.csv(tau_usual, "tau_matrix_usual_1000_200.csv", row.names = FALSE)
write.csv(var_usual, "var_matrix_usual_1000_200.csv", row.names = FALSE)
write.csv(tau_imp,   "tau_matrix_imp_1000_200.csv",   row.names = FALSE)
write.csv(var_imp,   "var_matrix_imp_1000_200.csv",   row.names = FALSE)

# calculate CP value 
true_tau <- 0

calc_cp <- function(tau_mat, var_mat) {
  se_mat <- sqrt(var_mat) / sqrt(N)
  lower  <- tau_mat - 1.96 * se_mat
  upper  <- tau_mat + 1.96 * se_mat
  colSums((true_tau >= lower) & (true_tau <= upper), na.rm = TRUE) / num_files
}

cp_usual <- calc_cp(tau_usual, var_usual)
cp_imp   <- calc_cp(tau_imp,   var_imp)

print(round(cp_usual, 3))
print(round(cp_imp,   3))