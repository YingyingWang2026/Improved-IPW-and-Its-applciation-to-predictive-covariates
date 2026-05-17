
library(tidyverse)
library(MASS)
library(BB)
library(numDeriv) 
library(conflicted)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)

N_vec       <- c(200,500,1000)
n_sim       <- 1000
n_cores     <- max(1L, detectCores() - 1L)
source("~/sup_sensitivity_function.R")
output_base <- "~/sup_results_eta_grid"


# data generation-----------------------------------------------------


sim_data <- function(N, eta) {
  
  Vc1 <- rnorm(N, mean = 0, sd = sqrt(2))   
  Vc  <- rnorm(N, mean = 0, sd = 1)         
  Vw  <- runif(N, min = -1, max = 1)        
  
  pz <- plogis(-0.1 + 0.4 * Vc1 - 0.2 * Vc)
  Z  <- rbinom(N, size = 1, prob = pz)
  
  eps0 <- rnorm(N, mean = 0, sd = 1)
  eps1 <- rnorm(N, mean = 0, sd = 1)
  
  Y0 <- 2 - 3 * Vc1 + Vc^3 - 3 * Vw + eps0
  Y1 <- 2 - 3 * Vc1 + 2 * Vc^3 - 3 * Vw + eta + eps1
  
  Y <- Z * Y1 + (1 - Z) * Y0
  
  return(list(
    Vc1 = Vc1,
    Vc  = Vc,
    Vw  = Vw,
    Z   = Z,
    Y   = as.vector(Y),
    Y1  = as.vector(Y1),
    Y0  = as.vector(Y0)
  ))
}




# propensity score coefficients
beta_Z1 <- c(-0.1, 0.4, -0.2)     
beta_Z2 <- c(-0.1, 0.4, -0.2, 0)   
beta_Z3 <- c(-0.1, -0.2)           
beta_Z4 <- c(-0.1, -0.2, 0)        


# Sensitivity parameters
eta_values  <- seq(-0.5, 0.5, length.out = 11)

e_par <- list(
  e_true_Vc   = c(-0.1,  0.4, -0.2),
  e_true_Vcw  = c(-0.1,  0.4, -0.2, 0),
  e_wrong_Vc  = c(-0.1,  0.4),
  e_wrong_Vcw = c(-0.1,  0.4,  0)
)

initial_para <- list(
  Q_true_Vcw  = c(1.5, -2,  1, -2, -2,  1),
  Q_true_Vc   = c(1.7, -2,  1, -2,  0.5),
  Q_wrong_Vcw = rep(1, 3),
  Q_wrong_Vc  = rep(1, 2)
)


#output-----------------------------------------------------
calculate_summary <- function(estimator, covarite_set, outcome_model, e_model, pscore, e_par) {
  if (estimator == "usual") {
    return(
      usual_est(
        para = initial_para,
        covarite_set = covarite_set,
        outcome_model = outcome_model,
        e_model = e_model,
        pscore = pscore,
        e_par = e_par
      )
    )
  } else if (estimator == "imp") {
    return(
      imp_est(
        para = initial_para,
        covarite_set = covarite_set,
        outcome_model = outcome_model,
        e_model = e_model,
        pscore = pscore,
        e_par = e_par
      )
    )
  } else if (estimator == "imp_std") {
    return(
      imp_est_std(
        para = initial_para,
        covarite_set = covarite_set,
        outcome_model = outcome_model,
        e_model = e_model,
        pscore = pscore,
        e_par = e_par
      )
    )
  } else if (estimator == "aipw_std") {
    return(
      aipw_est_std(
        para = initial_para,
        covarite_set = covarite_set,
        outcome_model = outcome_model,
        e_model = e_model,
        pscore = pscore,
        e_par = e_par
      )
    )
  }
  
  else {
    stop("Unknown estimator")
  }
}

# prepare for many settings

param_grid <- expand.grid(
  estimator = c("usual", "imp", "imp_std","aipw_std"),
  covarite_set = c("c", "cw"),
  outcome_model = c(TRUE, FALSE),
  e_model = c(TRUE, FALSE),
  pscore = c("known", "not known"),
  stringsAsFactors = FALSE
)

param_settings <- lapply(seq_len(nrow(param_grid)), function(i) {
  list(
    estimator = param_grid$estimator[i],
    covarite_set = param_grid$covarite_set[i],
    outcome_model = param_grid$outcome_model[i],
    e_model = param_grid$e_model[i],
    pscore = param_grid$pscore[i],
    e_par = e_par
  )
})



collect_summary_results <- function(param_settings) {
  
  
  res_df <- do.call(
    rbind,
    lapply(param_settings, function(s) {
      
      out <- tryCatch(
        {
          do.call(calculate_summary, s)
        },
        error = function(e) {
          warning(
            paste0(
              "Failed at estimator=", s$estimator,
              ", covarite_set=", s$covarite_set,
              ", outcome_model=", s$outcome_model,
              ", e_model=", s$e_model,
              ", pscore=", s$pscore,
              ". Error: ", e$message
            )
          )
          list(tau = NA, Var = NA)
        }
      )
      
      data.frame(
        estimator = s$estimator,
        covarite_set = s$covarite_set,
        outcome_model = s$outcome_model,
        e_model = s$e_model,
        pscore = s$pscore,
        tau = out$tau,
        AVar = out$Var,#####
        stringsAsFactors = FALSE
      )
    })
  )
  
  res_list <- split(res_df, res_df$estimator)
  
  return(list(
    result_list = res_list, 
    result_table = res_df   ))  
}




fns_to_export <- c(
  "usual_est", "imp_est", "imp_est_std", "aipw_est_std",
  "calculate_summary", "collect_summary_results",
  "sim_data",
  "VAR","b_fun", "B_fun","build_Q_design","score_fun", "fit_beta_std",
  "m_fun_baseline_ipw_known", "m_fun_baseline_ipw_specified",
  "m_fun_ipw_known", "m_fun_ipw_specified",
  "m_fun_imp_std_known", "m_fun_imp_std_specified",
  "m_fun_aipw_known", "m_fun_aipw_specified",
  "Q_Vcw_true", "Q_Vc_true", "Q_Vcw_wrong", "Q_Vc_wrong",
  "Q_Vcw_std_true", "Q_Vc_std_true", "Q_Vcw_std_wrong", "Q_Vc_std_wrong",
  "pscore_Vcw_true", "pscore_Vc_true", "pscore_Vcw_wrong", "pscore_Vc_wrong",
  "e_Vcw_true", "e_Vc_true", "e_Vcw_wrong", "e_Vc_wrong",
  "Gamma_ipw", "fai_ipw", "log_fun"
)


run_one_sim <- function(sim_id, N, eta, eta_id, param_settings, initial_para, e_par) {
  
  seed_value <- as.integer(paste0(eta_id, N, sim_id))
  set.seed(seed_value)
  
  data_now <- sim_data(N, eta)
  assign("Vc1", data_now[["Vc1"]], envir = globalenv())
  assign("Vc",  data_now[["Vc"]],  envir = globalenv())
  assign("Vw",  data_now[["Vw"]],  envir = globalenv())
  assign("Z",   data_now[["Z"]],   envir = globalenv())
  assign("Y",   data_now[["Y"]],   envir = globalenv())
  assign("N",   N,                  envir = globalenv())
  
  out <- tryCatch(
    collect_summary_results(param_settings),
    error = function(e) {
      message(sprintf("[eta=%.4f | N=%d | sim=%04d] ERROR: %s",
                      eta, N, sim_id, e$message))
      return(NULL)
    }
  )
  
  if (is.null(out)) return(NULL)
  
  df        <- out$result_table
  df$Var    <- df$AVar / N 
  df$bias   <- df$tau - eta 
  df$sim_id <- sim_id
  df$N      <- N
  df$eta    <- eta
  df
}

# Main Loop
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

for (eta_id in seq_along(eta_values)) {
  
  eta <- eta_values[eta_id]
  
  eta_label <- sprintf("eta_%+.4f", eta)
  eta_label <- gsub("\\+", "pos", eta_label)
  eta_label <- gsub("-", "neg", eta_label)
  eta_label <- gsub("\\.", "p", eta_label)
  
  output_root <- file.path(output_base, eta_label)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  
  for (N in N_vec) {
    
    out_dir <- file.path(output_root, paste0("N", N))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    cl <- makeCluster(n_cores, type = "PSOCK")
    registerDoParallel(cl)
    
    clusterExport(
      cl,
      varlist = c(fns_to_export, "run_one_sim",
                  "param_settings", "initial_para", "e_par"),
      envir = globalenv()
    )
    
    clusterEvalQ(cl, {
      library(tidyverse)
      library(MASS)
      library(BB)
      library(numDeriv)
      library(conflicted)
      library(dplyr)
    })
    
    foreach(
      sim_id = seq_len(n_sim),
      .packages  = c("tidyverse", "MASS", "BB", "numDeriv", "conflicted", "dplyr"),
      .errorhandling = "pass",
      .verbose = FALSE
    ) %dopar% {
      
      df_sim <- run_one_sim(
        sim_id        = sim_id,
        N             = N,
        eta           = eta,
        eta_id        = eta_id,
        param_settings = param_settings,
        initial_para  = initial_para,
        e_par         = e_par
      )
      
      if (!is.null(df_sim)) {
        fname <- file.path(out_dir, sprintf("sim_%04d.csv", sim_id))
        write.csv(df_sim, fname, row.names = FALSE)
      }
      
      invisible(NULL)
    }
    
    stopCluster(cl)
  }
  
  #output the summary results
  
  all_results <- lapply(N_vec, function(N) {
    
    out_dir <- file.path(output_root, paste0("N", N))
    files <- list.files(out_dir, pattern = "\\.csv$", full.names = TRUE)
    
    if (length(files) == 0) {
      warning(sprintf("No csv files found for eta = %.6f, N = %d", eta, N))
      return(NULL)
    }
    
    data_list <- lapply(files, read.csv)
    other_cols <- setdiff(colnames(data_list[[1]]), c("tau", "bias","AVar","Var"))
    
    mean_tau_var <- Reduce(`+`, lapply(data_list, function(df) {
      df[, c("tau", "bias","AVar","Var"), drop = FALSE]
    })) / length(files)
    
    cover_list <- lapply(data_list, function(df) {
      se    <- sqrt(df$AVar) / sqrt(N)
      lower <- df$tau - 1.96 * se
      upper <- df$tau + 1.96 * se
      as.numeric(eta >= lower & eta <= upper)
    })
    cp <- Reduce(`+`, cover_list) / length(files)
    
    other_data <- data_list[[1]][, other_cols, drop = FALSE]
    
    result <- cbind(other_data,mean_tau_var, cp = cp)
    result$N   <- N
    result$eta <- eta
    return(result)
  })
  
  all_results <- do.call(rbind, all_results)
  
  write.csv(
    all_results,
    file.path(output_root, "row_means_by_N.csv"),
    row.names = FALSE
  )
  
}
