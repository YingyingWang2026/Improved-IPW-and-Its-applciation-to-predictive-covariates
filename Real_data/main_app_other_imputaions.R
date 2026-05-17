library(mice)
library(missForest)
library(tidyverse)
library(MASS)
library(BB)
library(numDeriv)

set.seed(123)
df <- read.csv("~/data.csv")
missing_rates <- sapply(df, function(x) round(mean(is.na(x)), 3))
print(data.frame(column = names(df), missing_rate = missing_rates))

cov_raw <- c(
  "AGE", "SEX","NATION",  "MARITAL_STATUS", "HEIGHT", "WEIGHT", "BP_HIGH", "BP_LOW",
  "BMI","HYPERTENTION", "HYPERLIPIDEMIA", "A_S",
  "CEREBRAL_APOPLEXTY", "CAROTID_ARTERY_STENOSIS", "FLD",
  "CIRRHOSIS", "CLD", "PANCREATIC_DISEASE", "BILIARY_TRACT_DISEASE",
  "NEPHROPATHY", "RENAL_FALIURE", "NERVOUS_SYSTEM_DISEASE",
  "CHD", "MI", "CHF", "ARRHYTHMIAS", "RESPIRATORY_SYSTEM_DISEASE",
  "LEADDP", "HEMATONOSIS", "RHEUMATIC_IMMUNITY", "PREGNANT",
  "ENDOCRINE_DISEASE", "MEN", "PCOS", "DIGESTIVE_CARCINOMA",
  "UROLOGIC_NEOPLASMS", "GYNECOLGICAL_TUMOR", "BREAST_TUMOR",
  "LUNG_TUMOR", "INTRACRANIAL_TUMOR", "OTHER_TUMOR", "GLU",
  "GLU_2H", "GSP", "TG", "TC", "HDL_C",
  "LDL_C", "FBG", "UPR_24", "BU", "SCR", "UCR", "SUA", "HB",
  "CP", "INS", "PCV", "PLT", "ESR", "TBILI", "DBILI", "TP",
  "ALB", "LDH_L", "ALT", "AST", "GGT", "ALP", 
  "PT", "PTA", "APTT", "ALB_CR", "LPS", "CA199",
  "CRP", "IBILI", "GLO"
)

is_binary <- function(x) {
  ux <- unique(x[!is.na(x)])
  length(ux) <= 2 && all(ux %in% c(0, 1))
}

dat0 <- df[, c("label", "HBA1C", cov_raw)]
binary_vars    <- cov_raw[sapply(dat0[, cov_raw], is_binary)]
continuous_vars <- setdiff(cov_raw, binary_vars)

alpha          <- 0.05
VOTE_THRESHOLD <- 2

classify <- function(df_imp) {
  Z <- df_imp$label
  Y <- df_imp$HBA1C
  dat <- set_names(
    data.frame(Z = Z, Y = Y, df_imp[, cov_raw]),
    c("Z", "Y", paste0("X", 1:length(cov_raw)))
  )
  covs <- paste0("X", 1:length(cov_raw))
  
  sep <- covs[sapply(covs, function(x) {
    m <- suppressWarnings(glm(as.formula(paste("Z ~", x)), family = binomial(), data = dat))
    any(m$fitted.values < 1e-6 | m$fitted.values > 1 - 1e-6)
  })]
  covs <- setdiff(covs, sep)
  
  p_Z <- sapply(covs, function(x) {
    model <- glm(as.formula(paste("Z ~", x)), family = binomial(), data = dat)
    summary(model)$coefficients[x, "Pr(>|z|)"]
  })
  p_Y <- sapply(covs, function(x) {
    model <- lm(as.formula(paste("Y ~ Z +", x)), data = dat)
    summary(model)$coefficients[x, "Pr(>|t|)"]
  })
  
  idx <- as.integer(sub("X", "", covs))
  set_names(
    ifelse(p_Z < alpha & p_Y < alpha, "Confounder",
           ifelse(p_Z >= alpha & p_Y < alpha, "Predictor of Y", "others")),
    cov_raw[idx]
  )
}

# method 1: multiple imputation (mice) 
meth <- make.method(dat0)
meth["label"] <- meth["HBA1C"] <- ""
meth[binary_vars]    <- "logreg"
meth[continuous_vars] <- "pmm"
pred <- make.predictorMatrix(dat0)
pred["label", ] <- pred["HBA1C", ] <- 0
pred[cov_raw, c("label", "HBA1C")] <- 1

imp1 <- mice(dat0, m = 5, maxit = 20, method = meth,
             predictorMatrix = pred, seed = 2026, printFlag = FALSE)

mice_res  <- lapply(1:5, function(i) classify(complete(imp1, i)))
mice_vars <- unique(unlist(lapply(mice_res, names)))
res1 <- sapply(mice_vars, function(v) {
  votes <- sapply(mice_res, function(r) r[v])
  tb1 <- table(votes[!is.na(votes)])
  names(tb1)[which.max(tb1)]
})

#  method 2: missForest 
dat_rf <- dat0
dat_rf[binary_vars] <- lapply(dat_rf[binary_vars], factor)

imp2 <- missForest(dat_rf, maxiter = 10, ntree = 100, verbose = TRUE)$ximp
imp2[binary_vars] <- lapply(imp2[binary_vars], function(x) as.numeric(as.character(x)))

res2 <- classify(imp2)


all_vars <- unique(c(names(res1), names(res2)))
vote <- data.frame(
  var          = all_vars,
  n_confounder = sapply(all_vars, function(v) sum(c(res1[v], res2[v]) == "Confounder",   na.rm = TRUE)),
  n_predictor_Y = sapply(all_vars, function(v) sum(c(res1[v], res2[v]) == "Predictor of Y", na.rm = TRUE))
) %>% mutate(
  Final = case_when(
    n_confounder  >= VOTE_THRESHOLD ~ "Confounder",
    n_predictor_Y >= VOTE_THRESHOLD ~ "Predictor of Y",
    TRUE                            ~ "not_selected"
  )
)
print(vote)

stable_confounders <- vote$var[vote$Final == "Confounder"]
stable_predictor_Y <- vote$var[vote$Final == "Predictor of Y"]
print(stable_confounders)
print(stable_predictor_Y)


cov_avg <- Reduce("+", lapply(1:5, function(i) complete(imp1, i)[, cov_raw])) / 5

source("～/base_funtion.R")


make_dataset <- function(Z_vec, Y_vec, cov_mat) {
  dat <- set_names(
    data.frame(Z_vec, Y_vec, cov_mat),
    c("Z", "Y", paste0("X", seq_along(cov_raw)))
  )
  Xc <<- scale(as.matrix(dat[, paste0("X", which(cov_raw %in% stable_confounders)), drop = FALSE]))
  Xw <<- scale(as.matrix(dat[, paste0("X", which(cov_raw %in% stable_predictor_Y)),  drop = FALSE]))
  Z  <<- dat$Z
  Y  <<- dat$Y
  N  <<- length(Y)
}

run_estimators <- function(sfx) {
  set.seed(2026)
  z95 <- qnorm(0.975)
  
  est <- list(
    bc  = est_baseline_summary(para = initial_para, covarite_set = 'c',  type = 'ipw', e_par = e_par),
    bcw = est_baseline_summary(para = initial_para, covarite_set = 'cw', type = 'ipw', e_par = e_par),
    ic  = est_improved_summary (para = initial_para, covarite_set = 'c',  type = 'ipw', e_par = e_par),
    icw = est_improved_summary (para = initial_para, covarite_set = 'cw', type = 'ipw', e_par = e_par)
  )
  
  for (nm in names(est))
    cat(sfx, nm, "tau:", est[[nm]][["tau"]], "  Var:", est[[nm]][["Var"]], "\n")
  
  ci <- lapply(est, function(e) {
    se <- sqrt(e[["Var"]] / N)
    c(lower = e[["tau"]] - z95 * se, upper = e[["tau"]] + z95 * se)
  })
  for (nm in names(ci))
    cat(sfx, nm, "95% CI: [", ci[[nm]]["lower"], ",", ci[[nm]]["upper"], "]\n")
  
  list(
    tau    = c(bc = est$bc[["tau"]],  bcw = est$bcw[["tau"]],
               ic = est$ic[["tau"]],  icw = est$icw[["tau"]]),
    var    = c(bc = est$bc[["Var"]],  bcw = est$bcw[["Var"]],
               ic = est$ic[["Var"]],  icw = est$icw[["Var"]]),
    ci_bc  = ci$bc,  ci_bcw = ci$bcw,
    ci_ic  = ci$ic,  ci_icw = ci$icw
  )
}

flatten_results <- function(r, i) {
  s <- as.character(i)
  c(
    setNames(r$tau, paste0(c("baseline_ipw_tau_c", "baseline_ipw_tau_cw",
                             "improved_ipw_Vc_tau_c", "improved_ipw_Vcw_tau_cw"), s)),
    setNames(r$var, paste0(c("baseline_ipw_var_c", "baseline_ipw_var_cw",
                             "improved_ipw_Vcw_var_c", "improved_ipw_Vcw_var_cw"), s)),
    setNames(c(r$ci_bc, r$ci_bcw, r$ci_ic, r$ci_icw),
             paste0(c("ci_lower1","ci_upper1","ci_lower2","ci_upper2",
                      "ci_lower3","ci_upper3","ci_lower4","ci_upper4"), s))
  )
}


make_dataset(
  Z_vec   = complete(imp1, 1)$label,
  Y_vec   = complete(imp1, 1)$HBA1C,
  cov_mat = cov_avg
)
res_imp1 <- run_estimators("imp1")


make_dataset(
  Z_vec   = imp2$label,
  Y_vec   = imp2$HBA1C,
  cov_mat = imp2[, cov_raw]
)
res_imp2 <- run_estimators("imp2")


all_vals <- c(flatten_results(res_imp1, 1), flatten_results(res_imp2, 2))
write.csv(
  data.frame(variable = names(all_vals), value = all_vals),
  "imp_app_result.csv",
  row.names = FALSE
)
