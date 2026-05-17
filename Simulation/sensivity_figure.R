
library(tidyverse)

output_base <- "~/sup_results_eta_grid"
N_vec       <- c(200, 500, 1000)  
eta_dirs <- list.dirs(output_base, recursive = FALSE, full.names = TRUE)

# prepare for figures-----------------------------------------------------
for (eta_dir in eta_dirs) {
  eta_label <- basename(eta_dir)
  
  results_by_N <- lapply(N_vec, function(N) {
    
    out_dir <- file.path(eta_dir, paste0("N", N))
    if (!dir.exists(out_dir)) return(NULL)
    
    files <- list.files(out_dir, pattern = "\\.csv$", full.names = TRUE)
    if (length(files) == 0) return(NULL)
    
    data_list <- lapply(files, read.csv)
    
    has_na  <- sapply(data_list, function(df) any(is.na(df)))
    data_list <- data_list[!has_na]
    if (length(data_list) == 0) return(NULL)
    
    other_cols <- setdiff(colnames(data_list[[1]]), c("tau", "bias","AVar","Var"))
    tau_mat    <- do.call(cbind, lapply(data_list, `[[`, "tau"))
    var_mat    <- do.call(cbind, lapply(data_list, `[[`, "Var"))
    Avar_mat    <- do.call(cbind, lapply(data_list, `[[`, "AVar"))
    bias_mat    <- do.call(cbind, lapply(data_list, `[[`, "bias"))
    
    cbind(
      data_list[[1]][, other_cols, drop = FALSE],
      data.frame(tau = rowMeans(tau_mat), Var = rowMeans(var_mat), 
                 Avar = rowMeans(Avar_mat),  bias = rowMeans(bias_mat))
      
    )
  })
  
  results_by_N <- Filter(Negate(is.null), results_by_N)
  if (length(results_by_N) == 0) next
  
  out_csv <- file.path(eta_dir, "row_means_by_N_fig.csv")
  write.csv(do.call(rbind, results_by_N), out_csv, row.names = FALSE)
}

files <- list.files(
  output_base,
  pattern = "row_means_by_N_fig\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

res_all <- files |>
  lapply(read.csv) |>
  bind_rows()


# mse figure-----------------------------------------------------
mse_df <- res_all |>
  dplyr::filter(
    N == 200,# replace 500,1000
    estimator %in% c("usual", "imp", "aipw_std")
  ) |>
  mutate(
    estimator = recode(
      estimator,
      "usual" = "IPW",
      "imp" = "Improved IPW",
      "aipw_std" = "AIPW"
    ),
    
    bias = bias,
    var_tau = Var,
    
    mse = bias^2 + var_tau,
    
    qe_label = paste0(
      "Q=", ifelse(outcome_model, "correct", "wrong"),
      ", e=", ifelse(e_model, "correct", "wrong")
    ),
    
    row_label = paste0(
      ifelse(covarite_set == "c", "C", "CW"),
      ", PS=", pscore
    ),
    
    case_label = paste0(row_label, "\n", qe_label),
    
    qe_order = case_when(
      outcome_model == TRUE  & e_model == TRUE  ~ 1,
      outcome_model == FALSE & e_model == TRUE  ~ 2,
      outcome_model == TRUE  & e_model == FALSE ~ 3,
      outcome_model == FALSE & e_model == FALSE ~ 4
    ),
    
    row_order = case_when(
      covarite_set == "c"  & pscore == "known"     ~ 1,
      covarite_set == "c"  & pscore == "not known" ~ 2,
      covarite_set == "cw" & pscore == "known"     ~ 3,
      covarite_set == "cw" & pscore == "not known" ~ 4
    ),
    
    facet_order = row_order * 10 + qe_order
  ) |>
  arrange(facet_order) |>
  mutate(
    case_label = factor(case_label, levels = unique(case_label))
  )

p_mse <- ggplot(
  mse_df,
  aes(
    x = eta,
    y = mse,
    color = estimator,
    group = estimator,
    linetype = estimator,
    shape = estimator
  )
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_wrap(~ case_label, scales = "free_y", ncol = 4) +
  labs(
    x = expression(eta),
    y = "MSE",
    color = "Estimator",
    linetype = "Estimator",
    shape = "Estimator",
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

print(p_mse)

ggsave(
  file.path(output_base, "mse_comparison_N200_16cases.pdf"),
  p_mse,
  width = 14,
  height = 10,
  dpi = 300
)

# variance figure-----------------------------------------------------
var_df <- res_all |>
    dplyr::filter(
      N == 200,# replace 500,1000
      estimator %in% c("usual", "imp_std", "aipw_std")
    ) |>
    mutate(
      estimator = recode(
        estimator,
        "usual" = "IPW",
        "imp_std" = "Improved IPW(standard Q)",
        "aipw_std" = "AIPW"
      ),
      

      qe_label = paste0(
        "Q=", ifelse(outcome_model, "correct", "wrong"),
        ", e=", ifelse(e_model, "correct", "wrong")
      ),
      
  
      row_label = paste0(
        ifelse(covarite_set == "c", "C", "CW"),
        ", PS=", pscore
      ),

      case_label = paste0(row_label, "\n", qe_label),
      

      qe_order = case_when(
        outcome_model == TRUE  & e_model == TRUE  ~ 1,
        outcome_model == FALSE & e_model == TRUE  ~ 2,
        outcome_model == TRUE  & e_model == FALSE ~ 3,
        outcome_model == FALSE & e_model == FALSE ~ 4
      ),
      

      row_order = case_when(
        covarite_set == "c"  & pscore == "known"     ~ 1,
        covarite_set == "c"  & pscore == "not known" ~ 2,
        covarite_set == "cw" & pscore == "known"     ~ 3,
        covarite_set == "cw" & pscore == "not known" ~ 4
      ),
      

      facet_order = row_order * 10 + qe_order
    ) |>
    arrange(facet_order) |>
    mutate(
      case_label = factor(case_label, levels = unique(case_label))
    )
  

  p_var <- ggplot(
    var_df,
    aes(x = eta, y = Var,
        color = estimator,
        group = estimator,
        linetype = estimator,
        shape = estimator)
  ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    facet_wrap(~ case_label, scales = "free_y", ncol = 4) +
    labs(
      x = expression(eta),
      y = "Estimated variance",
      color = "Estimator",
      linetype = "Estimator",
      shape = "Estimator",

    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      strip.text = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    )
  
  print(p_var)
  
  ggsave(
    filename = file.path(output_base, "variance_comparison_N200_16cases.pdf"),
    plot = p_var,
    width = 14,
    height = 10,
    dpi = 300
  )


