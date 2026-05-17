
# Outcome model 

# Q_Vcw model
Q_Vcw <- function(para, Z) {
  list(partial_Q = cbind(1, Xc, Xw), Q = as.vector(cbind(1, Xc, Xw) %*% para))
}

# Q_Vc model 
Q_Vc <- function(para, Z) {
  list(partial_Q = cbind(1, Xc), Q = as.vector(cbind(1, Xc) %*% para))
}


# Propensity score model

# e_Vcw model
e_Vcw <- function(par, Z) {
  (Z * plogis(cbind(1, Xc,Xw) %*% par) + (1 - Z) * (1 - plogis(cbind(1, Xc,Xw) %*% par))) |> as.vector()
}

# e_Vc model
e_Vc <- function(par, Z) {
  (Z * plogis(cbind(1, Xc) %*% par) + (1 - Z) * (1 - plogis(cbind(1, Xc) %*% par))) |> as.vector()
}



# Log likelihood 
log_fun <- function(par, Z, e_V) {
  Z * log(e_V(par, 1)) + (1 - Z) * log(e_V(par, 0))
} 

# score functions 
score_fun <- function(par, Z, e_V = e_V) {
  jacobian(log_fun, x = par, Z = Z, e_V = e_V) |> colSums()
}

# Gamma, fai 
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


# m_fun 
m_fun_ipw_specified <- function(parameter, Q_summary, e_V) {
  q <- parameter |> names() |> startsWith('e_alpha') |> sum(); n_Q_beta <- (length(parameter) - (q + 1) ^ 2) / (q + 1)
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
  
  m5 <- (1 / (e_V(e_alpha, Z) ^ 2) * (partial_Qz + t(matrix(fai, nrow = n_Q_beta) %*% ginv(matrix(psi, nrow = q)) %*% t(partial_e))) * 
           as.vector(Y - Q - t(ksi %*% ginv(matrix(psi, nrow = q)) %*% t(partial_e))))
  
  m6 <- (-1) ^ (1 - Z) * (Y - Q) / e_V(e_alpha, Z) - tau
  
  m <- cbind(m1, m2, m3, m4, m5, m6)
  
  m
}



m_fun_baseline_ipw_specified <- function(parameter, Q_summary, e_V) {
  q <- parameter |> names() |> startsWith('e_alpha') |> sum()
  
  e_alpha <- parameter[1 : q] 
  tau <- parameter[length(parameter)]
  
  m1 <- jacobian(log_fun, x = e_alpha, Z = Z, e_V = e_V) 
  m2 <- (-1) ^ (1 - Z) * Y / e_V(e_alpha, Z) - tau
  m <- cbind(m1, m2)
  
  m
}





# variance
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
# improved IPW estimator
est_improved_summary <- function(para, covarite_set = 'cw', type = 'ipw', e_par = e_par) {
  
  para <- switch(
    covarite_set,
    'cw' = para[['Q_Vcw']],
    'c' = para[['Q_Vc']]
  ) 
  
  Q_summary <- switch(
    covarite_set,
    'cw' = Q_Vcw,
    'c' = Q_Vc
  )
  
  e_V <- switch(
    covarite_set,
    'cw' =  e_Vcw,
    'c' = e_Vc
  ) 
  
  e_par <- switch(
    covarite_set,
    'cw' = e_par[['e_Vcw']],
    'c' = e_par[['e_Vc']]
  ) 
  
  if(covarite_set == 'cw'){
    temp1 <- glm(Z ~ Xc+Xw , family = binomial('logit'))
  }else{
    temp1 <- glm(Z ~ Xc, family = binomial('logit'))
  }
  e_alpha <- temp1$coef
  
  eps<-0.05
  e_V_raw <- e_V
  e_V <- function(par,Z){
    ps <- plogis(cbind(1,if(covarite_set=='cw')cbind(Xc,Xw)else Xc)%*%par)
    ps <- pmin(pmax(ps,eps),1-eps)
    (Z*ps+(1-Z)*(1-ps))|>as.vector()}
  
  partial_e <- jacobian(e_V, x = e_alpha, Z = 1) |> t()
  S <- jacobian(log_fun, x = e_alpha, Z = Z, e_V = e_V)
  H <- t(S) %*% S / N
  
  eq_function <- function(para, Q_summary = Q_summary, e_V = e_V, type) {
    Q <- Q_summary(para, Z = Z)[['Q']]; Q1 <- Q_summary(para, Z = 1)[['Q']]; Q0 <- Q_summary(para, Z = 0)[['Q']]
    if(type == 'ipw') {
      gamma <- Gamma_ipw(para, Z = Z, Q_summary, e_V = e_V, e_alpha = e_alpha)
      g <- (gamma %*% ginv(H) %*% t(S)) |> t() |> as.vector()
      mean(((-1) ^ (1 - Z) * (Y - Q) / e_V(e_alpha, Z) - g)^2)
      
    } 
  }
  if(covarite_set == 'cw'){
    para <- lm(Y ~ Xc + Xw)$coefficients
  }else{ 
    para <- lm(Y ~ Xc)$coefficients
  }
  
  temp2 <- BBoptim(par = para, fn = eq_function, type = type, Q_summary = Q_summary,
                   e_V = e_V, control=list(maxit = 100, trace = TRUE, ftol = 1e-02), method = 2)
  
  Q_beta <- temp2$par
  Q1 <- Q_summary(Q_beta, Z = 1)[['Q']]; Q0 <- Q_summary(Q_beta, Z = 0)[['Q']]
  if(type == 'ipw') {
    tau <- (Z * (Y - Q1) / e_V(e_alpha, Z = 1) - (1 - Z) * (Y - Q0) / e_V(e_alpha, Z = 0)) |> mean()
    fai <- fai_ipw(Q_beta, Z = Z, Q_summary, e_V = e_V, e_alpha = e_alpha)
    gamma <- Gamma_ipw(Q_beta, Z = Z, Q_summary, e_V = e_V, e_alpha = e_alpha) 
    g <- (gamma %*% ginv(H) %*% t(S)) |> t() |> as.vector() 
    Var <- mean(((-1) ^ (1 - Z) * (Y - Q1) / e_V(e_alpha, Z) - g)^2)  
  } 
  return(list(Q_beta = Q_beta, tau = tau, Var = Var))
}


# usual IPW estimator
est_baseline_summary <- function(para, covarite_set = 'cw', type = 'ipw', e_par = e_par) {
  
  para <- switch(
    covarite_set,
    'cw' = para[['Q_Vcw']],
    'c' = para[['Q_Vc']]
  ) 
  
  Q_summary <- switch(
    covarite_set,
    'cw' = Q_Vcw,
    'c' = Q_Vc
  )
  
  e_V <- switch(
    covarite_set,
    'cw' =  e_Vcw,
    'c' = e_Vc
  ) 
  
  e_par <- switch(
    covarite_set,
    'cw' = e_par[['e_Vcw']],
    'c' = e_par[['e_Vc']]
  ) 
  
  
  eq_function <- function(para, Q_summary = Q_summary) {
    Q <- Q_summary(para, Z = Z)[['Q']]; partial_Qz <- Q_summary(para, Z = Z)[['partial_Q']]
    ((Y - Q) * partial_Qz) |> colMeans()
  }
  if(covarite_set == 'cw'){
    temp1 <- glm(Z ~ Xc+Xw , family = binomial('logit'))
    Q_beta <- rep(0, ncol( cbind(1, Xc, Xw)))
  }else{
    Q_beta <- rep(0, ncol( cbind(1, Xc)))
    temp1 <- glm(Z ~ Xc, family = binomial('logit'))
  }
  e_alpha <- temp1$coef

  eps<-0.05
  e_V_raw <- e_V
  e_V <- function(par,Z){
    ps <- plogis(cbind(1,if(covarite_set=='cw')cbind(Xc,Xw)else Xc)%*%par)
    ps <- pmin(pmax(ps,eps),1-eps)
    (Z*ps+(1-Z)*(1-ps))|>as.vector()}
  S <- jacobian(log_fun, x = e_alpha, Z = Z, e_V = e_V)
  H <- t(S) %*% S / N
  
  
  
  
  if(type == 'ipw') {
    tau <- ((-1) ^ (1 - Z) * Y / e_V(e_alpha, Z)) |> mean()
    gamma <- Gamma_ipw(Q_beta, Z = Z, Q_summary, e_V = e_V, e_alpha = e_alpha) 
    g <- (gamma %*% ginv(H) %*% t(S)) |> t() |> as.vector() 
    Var <- mean(((-1) ^ (1 - Z) * (Y - 0) / e_V(e_alpha, Z) - g)^2) 
  } 
  return(list(Q_beta = Q_beta, tau = tau, Var = Var))
}


# Q_beta  initial value
initial_para <- list(
  Q_Vcw = c(rep(1, ncol(cbind(1, Xc, Xw)))), 
  Q_Vc = c(rep(1, ncol(cbind(1, Xc))))
)


e_par <- list(
  e_Vc = rep(0.01, 1 + ncol(Xc)),        
  e_Vcw = rep(0.01, 1 + ncol(Xc)+ncol(Xw)) 
)
