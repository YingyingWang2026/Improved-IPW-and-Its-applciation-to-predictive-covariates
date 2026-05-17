# 加载必要的包
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
####覆盖率不可能是1
N=1000
true_tau<--0.5
folder_path <- "/Users/yiyayiyayo/🧑‍🎓/IPW_adj/4.17-simpleseed/sup_results_-0.5/N1000"
csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

# 如果文件数量不是1000，给出提示
if(length(csv_files) != 1000) {
  warning(paste("找到", length(csv_files), "个CSV文件，不是1000个。"))
}

# 初始化一个空列表，用于存储每种情况的覆盖结果（每个实验一个逻辑向量）
all_results <- list()

# 遍历每个CSV文件
for (i in seq_along(csv_files)) {
  # 读取CSV文件
  df <- read_csv(csv_files[i], show_col_types = FALSE)
  
  # 检查必要的列是否存在
  if(!all(c("tau", "Var") %in% colnames(df))) {
    stop(paste("文件", csv_files[i], "中缺少 tau 或 Var 列"))
  }
  

  df <- df %>%
    mutate(
      se = sqrt(Var)/sqrt(N),                     # 标准误
      lower = tau - 1.96 * se,            # 置信区间下限
      upper = tau + 1.96 * se,            # 置信区间上限
      cover = (true_tau >= lower & true_tau <= upper)  # 是否落入区间
    )
  
  # 保存本次实验的覆盖结果（一个逻辑向量，长度=64）
  all_results[[i]] <- df$cover
}

# 将列表转换为数据框（1000行 × 64列，每行是一次实验，每列是一种情况）
results_df <- as.data.frame(do.call(rbind, all_results))
colnames(results_df) <- paste0("case_", 1:ncol(results_df))

# 计算每种情况的覆盖率（在1000次实验中的比例）
coverage_rates <- results_df %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(everything(), names_to = "Case", values_to = "CoverageRate")

# 输出覆盖率结果
print(coverage_rates)

# 可选：保存覆盖率结果到CSV
write_csv(coverage_rates, "/Users/yiyayiyayo/🧑‍🎓/IPW_adj/4.17-simpleseed/sup_results_-0.5/N1000_coverage_rates_64cases.csv")

# 可选：保存每个实验的详细覆盖结果（1000×64）
write_csv(results_df, "detailed_coverage_all_experiments.csv")

# 可视化覆盖率（直方图或点图）
library(ggplot2)
ggplot(coverage_rates, aes(x = Case, y = CoverageRate)) +
  geom_point() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "95% Confidence Interval Coverage Rate for 64 Cases (1000 Replicates)",
       y = "Coverage Rate", x = "Case")