
library(tidyverse)
library(MASS)
library(BB)
library(numDeriv)

set.seed(123)
df <- read.csv("～/data.csv")
missing_rates <- sapply(df, function(x) round(sum(is.na(x)) / length(x), 3))
print(data.frame(column = names(df), missing_rate = missing_rates))

cov_raw <- c(
  "AGE", "SEX", "NATION", "MARITAL_STATUS", "HEIGHT", "WEIGHT", "BP_HIGH", "BP_LOW",
  "BMI", "HYPERTENTION", "HYPERLIPIDEMIA", "A_S",
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

dat0 <- df[, c("label", "HBA1C", cov_raw)]
binary_vars     <- cov_raw[sapply(dat0[, cov_raw], function(x) {
  ux <- unique(x[!is.na(x)]); length(ux) <= 2 && all(ux %in% c(0, 1))
})]
continuous_vars <- setdiff(cov_raw, binary_vars)

#imputation 
imp_vars <- function(dat, binary_vars, continuous_vars) {
  dat_imp <- dat
  for (v in binary_vars) {
    miss_idx <- which(is.na(dat_imp[[v]]))
    if (length(miss_idx) == 0) next
    dat_imp[[v]][miss_idx] <- rbinom(length(miss_idx), 1, mean(dat_imp[[v]], na.rm = TRUE))
  }
  for (v in continuous_vars) {
    miss_idx <- which(is.na(dat_imp[[v]]))
    if (length(miss_idx) == 0) next
    dat_imp[[v]][miss_idx] <- rnorm(length(miss_idx),
                                    mean(dat_imp[[v]], na.rm = TRUE),
                                    sd(dat_imp[[v]],   na.rm = TRUE))
  }
  dat_imp
}

dat_imp <- imp_vars(dat0, binary_vars, continuous_vars)
Z <- dat_imp$label
Y <- dat_imp$HBA1C
names(dat_imp) <- c("Z", "Y", paste0("X", seq_along(cov_raw)))

#covariate classification 
alpha      <- 0.05
covariates <- paste0("X", seq_along(cov_raw))

p_Z <- sapply(covariates, function(x) {
  glm(as.formula(paste("Z ~", x)), family = binomial(), data = dat_imp) |>
    summary() |> (\(s) s$coefficients[x, "Pr(>|z|)"])()
})
print(p_Z)

p_Y <- sapply(covariates, function(x) {
  lm(as.formula(paste("Y ~ Z +", x)), data = dat_imp) |>
    summary() |> (\(s) s$coefficients[x, "Pr(>|t|)"])()
})
print(p_Y)

classification <- data.frame(
  Covariate = covariates,
  Type = ifelse(p_Z < alpha & p_Y < alpha, "Confounder",
                ifelse(p_Y < alpha, "Predictor of Y",
                       ifelse(p_Z < alpha, "Predictor of Z", "Irrelevant")))
)
print(classification)

Xc <- scale(as.matrix(dat_imp[, classification$Covariate[classification$Type == "Confounder"],     drop = FALSE]))
Xw <- scale(as.matrix(dat_imp[, classification$Covariate[classification$Type == "Predictor of Y"], drop = FALSE]))
N  <- length(Y)

# estimation 
source("～/base_funtion.R")
set.seed(1234)
z95 <- qnorm(0.975)

est <- list(
  bc  = est_baseline_summary(para = initial_para, covarite_set = 'c',  type = 'ipw', e_par = e_par),
  bcw = est_baseline_summary(para = initial_para, covarite_set = 'cw', type = 'ipw', e_par = e_par),
  ic  = est_improved_summary (para = initial_para, covarite_set = 'c',  type = 'ipw', e_par = e_par),
  icw = est_improved_summary (para = initial_para, covarite_set = 'cw', type = 'ipw', e_par = e_par)
)

for (nm in names(est))
  cat(nm, "tau:", est[[nm]][["tau"]], "  Var:", est[[nm]][["Var"]], "\n")

ci <- lapply(est, function(e) {
  se <- sqrt(e[["Var"]] / N)
  c(lower = e[["tau"]] - z95 * se, upper = e[["tau"]] + z95 * se)
})
for (nm in names(ci))
  cat(nm, "95% CI: [", ci[[nm]]["lower"], ",", ci[[nm]]["upper"], "]\n")