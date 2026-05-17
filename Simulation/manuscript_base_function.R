# base_function-----------------------------------------------------

# propensity score model

pscore_Vc_true <- function(Z) {
  (Z * plogis(cbind(1, Vc1, Vc) %*% beta_Z1) + (1 - Z) * (1 - plogis(cbind(1, Vc1, Vc) %*% beta_Z1))) |> as.vector()
}

pscore_Vcw_true <- function(Z) {
  (Z * plogis(cbind(1, Vc1, Vc, Vw) %*% beta_Z2) + (1 - Z) * (1 - plogis(cbind(1, Vc1, Vc, Vw) %*% beta_Z2))) |> as.vector()
}

pscore_Vc_wrong <- function(Z) { 
  (Z * plogis(cbind(1, Vc) %*% beta_Z3) + (1 - Z) * (1 - plogis(cbind(1,Vc) %*% beta_Z3))) |> as.vector()
}
pscore_Vcw_wrong <-function(Z){
  (Z * plogis(cbind(1, Vc, Vw) %*% beta_Z4) + (1 - Z) * (1 - plogis(cbind(1, Vc , Vw) %*% beta_Z4))) |> as.vector()
}


# estimated e(x,alpha)

e_Vc_true <- function(par, Z) {
  (Z * plogis(cbind(1, Vc1, Vc) %*% par) + (1 - Z) * (1 - plogis(cbind(1, Vc1, Vc) %*% par))) |> as.vector()
}

e_Vcw_true <- function(par, Z) {
  (Z * plogis(cbind(1, Vc1, Vc,Vw) %*% par) + (1 - Z) * (1 - plogis(cbind(1, Vc1, Vc, Vw) %*% par))) |> as.vector()
}


##### wrong e(V,Z) model where V %in% {Vc, Vcw}

e_Vc_wrong <- function(par, Z) {
  (Z * pnorm(cbind(1, Vc) %*% par) + (1 - Z) * (1 - pnorm(cbind(1,  Vc) %*% par))) |> as.vector()
}

e_Vcw_wrong <- function(par, Z) {
  (Z * pnorm(cbind(1, Vc, Vw) %*% par) + (1 - Z) * (1 - pnorm(cbind(1, Vc, Vw) %*% par))) |> as.vector()
}



log_fun <- function(par, Z, e_V) {
  Z * log(e_V(par, 1)) + (1 - Z) * log(e_V(par, 0))
}

score_fun <- function(par, Z, e_V) {
  jacobian(log_fun, x = par, Z = Z, e_V = e_V) |> colSums()
}


Gamma_ipw <- function(para, Z, Q_summary, e_V = e_V, e_alpha) {
  Q1 <- Q_summary(para, Z = 1)[['Q']]; Q0 <- Q_summary(para, Z = 0)[['Q']]
  fai <- function(e_alpha, Z) {
    Z * Y / e_V(e_alpha, Z = 1) - (1 - Z) * Y / e_V(e_alpha, Z = 0) - Z / e_V(e_alpha, Z = 1) * Q1 + 
      (1 - Z) / e_V(e_alpha, Z = 0) * Q0
  }
  -jacobian(fai, x = e_alpha, Z = Z) |> colMeans()
} 

fai_ipw <- function(para, Z, Q_summary, e_V = e_V, e_alpha) {
  partial_e <- jacobian(e_V, x = e_alpha, Z = 1) 
  partial_Qz <- Q_summary(para, Z = Z)[['partial_Q']] 
  M <- partial_e / (e_V(e_alpha, Z = Z) ^ 2)
  L <- vector('list', N)   
  for (i in 1:N) {
    L[[i]] <- partial_Qz[i, ] %o% M[i, ]
  }
  -(L |> reduce(`+`)) / N 
}


# outcome model

Q_Vcw_true <- function(para, Z) {
  
  list(partial_Q = cbind(1, Vc1, Vc^3, Vc, Vw, Z * Vc^3), Q = as.vector(cbind(1, Vc1, Vc^3, Vc, Vw, Z * Vc^3) %*% para))
}
Q_Vc_true <- function(para, Z) {
  
  list(partial_Q = cbind(1, Vc1, Vc^3, Vc, Z * Vc^3), Q = as.vector(cbind(1, Vc1, Vc^3, Vc, Z * Vc^3) %*% para)) 
}

#####  wrong outcome model Q(V,Z) where %in% {Vc, Vcw}

Q_Vcw_wrong <- function(para, Z) {
  list(partial_Q = cbind(1, Z*Vc, Z*Vw), Q = as.vector(cbind(1, Z*Vc, Z*Vw) %*% para))
}
Q_Vc_wrong <- function(para, Z) {
  list(partial_Q = cbind(1, Z*Vc), Q = as.vector(cbind(1,  Z*Vc) %*% para)) 
}




### m_function

m_fun_baseline_ipw_known <- function(parameter, Q_summary, e_V) {
  
  tau <- parameter
  
  m1 <- ((-1) ^ (1 - Z) * Y / e_V(Z) - tau) |> as.matrix()
  
  m1
} 




m_fun_baseline_ipw_specified <- function(parameter, Q_summary, e_V) {
  q <- parameter |> names() |> startsWith('e') |> sum()
  
  e_alpha <- parameter[1 : q] 
  tau <- parameter[length(parameter)]
  
  m1 <- jacobian(log_fun, x = e_alpha, Z = Z, e_V = e_V) 
  m2 <- (-1) ^ (1 - Z) * Y / e_V(e_alpha, Z) - tau
  m <- cbind(m1, m2)
  
  m
} 


m_fun_ipw_known <- function(parameter, Q_summary, e_V) {
  para <-  parameter[1 : (length(parameter) - 1)]
  tau <- parameter[length(parameter)]
  Q <- Q_summary(para, Z = Z)[['Q']]; partial_Qz <- Q_summary(para, Z = Z)[['partial_Q']]
  
  m1 <- (Y - Q) * partial_Qz / e_V(Z) ^ 2 
  m2 <- (-1) ^ (1 - Z) * (Y - Q) / e_V(Z) - tau
  m <- cbind(m1, m2)
  m
}

m_fun_ipw_specified_origin <- function(parameter, Q_summary, e_V) {
  q <- parameter |> names() |> startsWith('e') |> sum(); n_Q_beta <- (length(parameter) - (q + 1) ^ 2) / (q + 1)
  e_alpha <- parameter[1 : q]; 
  psi <- parameter[(q + 1) : (q + q ^ 2)]; 
  fai <- parameter[(q + q ^ 2 + 1) : (q + q ^ 2 + q * n_Q_beta)]
  ksi <- parameter[(q + q ^ 2 + q * n_Q_beta + 1) : (2 * q + q ^ 2 + q * n_Q_beta)]
  para <- parameter[(2 * q + q ^ 2 + q * n_Q_beta + 1) : (2 * q + q ^ 2 + q * n_Q_beta + n_Q_beta)]
  tau <- parameter[2 * q + q ^ 2 + q * n_Q_beta + n_Q_beta + 1]
  
  m1 <- jacobian(log_fun, x = e_alpha, Z = Z, e_V = e_V) 
  m2 <- t(psi - (m1 |> apply(1,as_tibble) |> map(~ .$value) |> map(~ . %*% t(.)) |> map(c) |> reduce(cbind))) 
  
  partial_e <- jacobian(e_V, x = e_alpha, Z = 1) 
  partial_Qz <- Q_summary(para, Z = Z)[['partial_Q']] 
  M <- partial_e / (e_V(e_alpha, Z = Z) ^ 2)
  L <- vector('list', N)
  for (i in 1:N) {
    L[[i]] <- partial_Qz[i, ] %o% M[i, ]
  }
  m3 <- t(fai + (L |> map(c) |> reduce(cbind)))
  
  Q1 <- Q_summary(para, Z = 1)[['Q']]; Q0 <- Q_summary(para, Z = 0)[['Q']]; Q <- Q_summary(para, Z = Z)[['Q']]
  Fai <- function(e_alpha, Z) {
    Z * Y / e_V(e_alpha, Z = 1) - (1 - Z) * Y / e_V(e_alpha, Z = 0) - Z / e_V(e_alpha, Z = 1) * Q1 + 
      (1 - Z) / e_V(e_alpha, Z = 0) * Q0
  }
  m4 <- t(ksi + t(jacobian(Fai, x = e_alpha, Z = Z)))
  
  m5 <- (1 / (e_V(e_alpha, Z) ^ 2) * (partial_Qz + t(matrix(fai, nrow = n_Q_beta) %*% solve(matrix(psi, nrow = q)) %*% t(partial_e))) * 
           as.vector(Y - Q - t(ksi %*% solve(matrix(psi, nrow = q)) %*% t(partial_e))))
  
  m6 <- (-1) ^ (1 - Z) * (Y - Q) / e_V(e_alpha, Z) - tau
  
  m <- cbind(m1, m2, m3, m4, m5, m6)
  
  m
}



m_fun_ipw_specified <- function(parameter, Q_summary, e_V) {
  q <- parameter |> names() |> startsWith('e') |> sum(); n_Q_beta <- (length(parameter) - q - 1)
  
  e_alpha <- parameter[1 : q]; 
  Q_beta <- parameter[(q + 1) : (q + n_Q_beta)]
  tau <- parameter[q + n_Q_beta + 1]
  
  S <- jacobian(log_fun, x = e_alpha, Z = Z, e_V = e_V)
  H <- t(S) %*% S / N
  fai <- fai_ipw(Q_beta, Z = Z, Q_summary, e_V = e_V, e_alpha = e_alpha)
  ksi <- Gamma_ipw(Q_beta, Z = Z, Q_summary, e_V = e_V, e_alpha = e_alpha)
  
  
  m1 <- jacobian(log_fun, x = e_alpha, Z = Z, e_V = e_V) 
  
  partial_e <- jacobian(e_V, x = e_alpha, Z = 1) 
  partial_Qz <- Q_summary(Q_beta, Z = Z)[['partial_Q']]
  
  Q1 <- Q_summary(Q_beta, Z = 1)[['Q']]; Q0 <- Q_summary(Q_beta, Z = 0)[['Q']]; Q <- Q_summary(Q_beta, Z = Z)[['Q']]
  
  
  m5 <- (1 / (e_V(e_alpha, Z) ^ 2) * (partial_Qz + t(fai %*% solve(H) %*% t(partial_e))) * 
           as.vector(Y - Q - t(ksi %*% solve(H) %*% t(partial_e))))
  
  m6 <- (-1) ^ (1 - Z) * (Y - Q) / e_V(e_alpha, Z) - tau
  
  
  m <- cbind(m1, m5, m6)
  
  m
}




### variance function
b_fun <- function(m_fun, parameter = parameter, Q_summary, e_V) {
  b <- apply(m_fun(parameter, Q_summary, e_V), 2, mean)
  b
}
B_fun <- function(m_fun, parameter = parameter, Q_summary, e_V) {
  B <- jacobian(b_fun, x = parameter, m_fun = m_fun, Q_summary = Q_summary, e_V = e_V)
  B
}
VAR <- function(m_fun, parameter = parameter, Q_summary, e_V) {
  B <- B_fun(m_fun, parameter = parameter, Q_summary = Q_summary, e_V = e_V)
  m <- m_fun(parameter = parameter, Q_summary = Q_summary, e_V = e_V)
  D <- t(m) %*% m / N
  VAR <- ginv(B) %*% D %*% t(ginv(B))
  VAR[ncol(VAR), ncol(VAR)]
}



# estimating function 

# usual IPW estimator
usual_est <- function(para, data, covarite_set = 'cw', outcome_model = FALSE, type = 'ipw', pscore = 'known', e_model = TRUE, e_par = e_par) {
  para <- if(outcome_model) {
    switch(
      covarite_set,
      'cw' = para[['Q_true_Vcw']],
      'c' = para[['Q_true_Vc']]
    ) 
  } else {
    switch(
      covarite_set,
      'cw' = para[['Q_wrong_Vcw']],
      'c' = para[['Q_wrong_Vc']]
    )
  }
  
  Q_summary <- if(outcome_model) {
    switch(
      covarite_set,
      'cw' = Q_Vcw_true,
      'c' = Q_Vc_true
    )
  } else {
    switch(
      covarite_set, 
      'cw' = Q_Vcw_wrong,
      'c' = Q_Vc_wrong
    )
  }
  
  
  eq_function <- function(para, Q_summary = Q_summary) {
    Q <- Q_summary(para, Z = Z)[['Q']]; partial_Qz <- Q_summary(para, Z = Z)[['partial_Q']]
    ((Y - Q) * partial_Qz) |> colMeans()
  }
  
  if(pscore == 'known') {
    if(e_model){
      e_V <-  switch(
        covarite_set,
        'cw' = pscore_Vcw_true,
        'c' = pscore_Vc_true)
    } else {
      e_V <-  switch(
        covarite_set,
        'cw' = pscore_Vcw_wrong,
        'c' = pscore_Vc_wrong)
    }
    
    if(type == 'ipw') {
      Q_beta <- NA
      tau <- ((-1) ^ (1 - Z) * Y / e_V(Z)) |> mean()
      
      parameter <- tau
      Var <- VAR(m_fun_baseline_ipw_known, parameter = parameter, Q_summary = Q_summary, e_V = e_V)
      
    } 
    
    
  } else {
    if(e_model) {
      e_V <- switch(
        covarite_set,
        'cw' =  e_Vcw_true,
        'c' = e_Vc_true
      ) 
      e_par <- switch(
        covarite_set,
        'cw' = e_par[['e_true_Vcw']],
        'c' = e_par[['e_true_Vc']]
      ) 
    } else {
      e_V <- switch(
        covarite_set,
        'cw' =  e_Vcw_wrong,
        'c' = e_Vc_wrong
      ) 
      e_par <- switch(
        covarite_set,
        'cw' = e_par[['e_wrong_Vcw']],
        'c' = e_par[['e_wrong_Vc']]
      )
    }
    #11
    temp <- tryCatch({
      BBsolve(par = e_par, fn = score_fun, Z = Z, e_V = e_V)
    }, error = function(e) {
      warning("BBsolve failed: ", conditionMessage(e))
      return(NULL)
    })
    
    if (is.null(temp) || is.null(temp$convergence) || temp$convergence != 0) {
      warning("BBsolve did not converge. Returning NA.")
      return(list(Q_beta = NA, tau = NA, Var = NA))
    }
    
    
    
    e_alpha <- temp$par
    S <- jacobian(log_fun, x = e_alpha, Z = Z, e_V = e_V)
    H <- t(S) %*% S / N
    
    if(type == 'ipw') {
      Q_beta <- NA
      if(temp$convergence == 0) {
        tau <- ((-1) ^ (1 - Z) * Y / e_V(e_alpha, Z)) |> mean()
        
        parameter <- c(e_alpha = e_alpha, tau)
        Var <- VAR(m_fun_baseline_ipw_specified, parameter = parameter, Q_summary = Q_summary, e_V = e_V)
        
      } else {
        tau <- Var <- NA
      }
    }
    
  }
  return(list(Q_beta = Q_beta, tau = tau, Var = Var))
}



# improved IPW estimator 
imp_est <- function(para, data, covarite_set = 'cw', outcome_model = FALSE, type = 'ipw', pscore = 'known', e_model = TRUE, e_par = e_par) {
  para <- if(outcome_model) {
    switch(
      covarite_set,
      'cw' = para[['Q_true_Vcw']],
      'c' = para[['Q_true_Vc']]
    ) 
  } else {
    switch(
      covarite_set,
      'cw' = para[['Q_wrong_Vcw']],
      'c' = para[['Q_wrong_Vc']]
    )
  }
  
  Q_summary <- if(outcome_model) {
    switch(
      covarite_set,
      'cw' = Q_Vcw_true,
      'c' = Q_Vc_true
    )
  } else {
    switch(
      covarite_set, 
      'cw' = Q_Vcw_wrong,
      'c' = Q_Vc_wrong
    )
  }
  
  
  if(pscore == 'known') {
    if(e_model){
      e_V <-  switch(
        covarite_set,
        'cw' = pscore_Vcw_true,
        'c' = pscore_Vc_true)
    } else {
      e_V <-  switch(
        covarite_set,
        'cw' = pscore_Vcw_wrong,
        'c' = pscore_Vc_wrong)
    }
    
    eq_function <- function(para, Q_summary = Q_summary, type, e_V = e_V) {
      Q <- Q_summary(para, Z = Z)[['Q']]; partial_Qz <- Q_summary(para, Z = Z)[['partial_Q']]
      partial_Q1 <- Q_summary(para, Z = 1)[['partial_Q']]
      partial_Q0 <- Q_summary(para, Z = 0)[['partial_Q']]
      
      if(type == 'ipw') {
        ((Y - Q) * partial_Qz / e_V(Z) ^ 2) |> colMeans() 
      } 
    }
    
    temp <- BBsolve(par = para, fn = eq_function, type = type, Q_summary = Q_summary, e_V = e_V, control=list(maxit = 100, trace = TRUE), method = 2)
    if (outcome_model | e_model) {
      if(temp$convergence == 0){
        Q_beta <- temp$par
        Q1 <- Q_summary(Q_beta, Z = 1)[['Q']]; Q0 <- Q_summary(Q_beta, Z = 0)[['Q']]
        if(type == 'ipw') {
          tau <- (Z * (Y - Q1) / e_V(Z = 1) - (1 - Z) * (Y - Q0) / e_V(Z = 0)) |> mean() 
          parameter <- c(Q_beta, tau)
          Var <- VAR(m_fun_ipw_known, parameter = parameter, Q_summary = Q_summary, e_V = e_V)
        } 
      } 
      else{Q_beta <- rep(NA, length(para))  
      tau <- NA
      Var <- NA}
    } else {
      Q_beta <- temp$par
      Q1 <- Q_summary(Q_beta, Z = 1)[['Q']]; Q0 <- Q_summary(Q_beta, Z = 0)[['Q']]
      if(type == 'ipw') {
        tau <- (Z * (Y - Q1) / e_V(Z = 1) - (1 - Z) * (Y - Q0) / e_V(Z = 0)) |> mean() 
        parameter <- c(Q_beta, tau)
        Var <- VAR(m_fun_ipw_known, parameter = parameter, Q_summary = Q_summary, e_V = e_V)
      } 
    }
    
  } else {
    if(e_model) {
      e_V <- switch(
        covarite_set,
        'cw' =  e_Vcw_true,
        'c' = e_Vc_true
      ) 
      e_par <- switch(
        covarite_set,
        'cw' = e_par[['e_true_Vcw']],
        'c' = e_par[['e_true_Vc']]
      ) 
    } else {
      e_V <- switch(
        covarite_set,
        'cw' =  e_Vcw_wrong,
        'c' = e_Vc_wrong
      ) 
      e_par <- switch(
        covarite_set,
        'cw' = e_par[['e_wrong_Vcw']],
        'c' = e_par[['e_wrong_Vc']]
      )
    }
    
    temp1 <- BBsolve(par = e_par, fn = score_fun, Z = Z, e_V = e_V)
    e_alpha <- temp1$par
    partial_e <- jacobian(e_V, x = e_alpha, Z = 1) |> t()
    S <- jacobian(log_fun, x = e_alpha, Z = Z, e_V = e_V)
    H <- t(S) %*% S / N
    
    eq_function <- function(para, Q_summary = Q_summary, e_V = e_V, type) {
      Q <- Q_summary(para, Z = Z)[['Q']]; Q1 <- Q_summary(para, Z = 1)[['Q']]; Q0 <- Q_summary(para, Z = 0)[['Q']]
      partial_Qz <- Q_summary(para, Z = Z)[['partial_Q']]
      partial_Q1 <- Q_summary(para, Z = 1)[['partial_Q']]
      partial_Q0 <- Q_summary(para, Z = 0)[['partial_Q']]
      if(type == 'ipw') {
        gamma <- Gamma_ipw(para, Z = Z, Q_summary, e_V = e_V, e_alpha = e_alpha)
        fai <- fai_ipw(para, Z = Z, Q_summary, e_V = e_V, e_alpha = e_alpha)
        bg <- (fai %*% solve(H) %*% partial_e) |> t()
        g <- (gamma %*% solve(H) %*% partial_e) |> t() |> as.vector()
        # estimating functions
        (1 / (e_V(e_alpha, Z) ^ 2) * (partial_Qz +  bg) * (Y - Q - g)) |> colMeans()
        
      } 
    }
    
    temp2 <- tryCatch({
      BBsolve(par = para, fn = eq_function, type = type, Q_summary = Q_summary, e_V = e_V,
              control = list(maxit = 1000, tol = 1e-02, trace = TRUE), method = 2)
    }, error = function(e) {
      warning("BBsolve temp2 (Q_beta) failed: ", conditionMessage(e))
      return(NULL)
    })
    
    if (is.null(temp2) || is.null(temp2$convergence) || temp2$convergence != 0) {
      warning("BBsolve temp2 (Q_beta) did not converge. Returning NA.")
      return(list(Q_beta = NA, tau = NA, Var = NA))
    }
    
    print(temp2$par)
    if (outcome_model | e_model) {
      if(temp1$convergence == 0 & temp2$convergence == 0) {
        Q_beta <- temp2$par
        Q1 <- Q_summary(Q_beta, Z = 1)[['Q']]; Q0 <- Q_summary(Q_beta, Z = 0)[['Q']]
        if(type == 'ipw') {
          tau <- (Z * (Y - Q1) / e_V(e_alpha, Z = 1) - (1 - Z) * (Y - Q0) / e_V(e_alpha, Z = 0)) |> mean()
          
          parameter <- c(e_alpha = e_alpha, Q_beta, tau)
          Var <- VAR(m_fun_ipw_specified, parameter = parameter, Q_summary = Q_summary, e_V = e_V) 
          
        }  
      } else {
        Q_beta <- temp2$par
        tau <- Var <- NA
      }
    } else {
      Q_beta <- temp2$par
      Q1 <- Q_summary(Q_beta, Z = 1)[['Q']]; Q0 <- Q_summary(Q_beta, Z = 0)[['Q']]
      if(type == 'ipw') {
        tau <- (Z * (Y - Q1) / e_V(e_alpha, Z = 1) - (1 - Z) * (Y - Q0) / e_V(e_alpha, Z = 0)) |> mean()
        
        parameter <- c(e_alpha = e_alpha, Q_beta, tau)
        Var <- VAR(m_fun_ipw_specified, parameter = parameter, Q_summary = Q_summary, e_V = e_V) 
        
      } 
    }
    
  }
  return(list(Q_beta = Q_beta, tau = tau, Var = Var))
}
