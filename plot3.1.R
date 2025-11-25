suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(rstatix)
})

compare_methods <- function(all_datasets,
                            metrics = c("accuracy", "weighted_f1"),
                            save_prefix = NULL) {
  stopifnot(is.list(all_datasets), length(all_datasets) >= 1)
  
  # --- 固定方法顺序与显示名 ---
  method_order <- c("scANMF", "scCATCH", "ScType", "scPred", "SingleR")
  name_map <- c(
    "result_scanmf" = "scANMF",
    "result_scCATCH" = "scCATCH",
    "result_sctype" = "ScType",
    "result_scpred" = "scPred",
    "result_sigleR" = "SingleR"
  )
  
  extract_one_dataset <- function(ds_list, dataset_name) {
    method_names <- names(ds_list)
    
    cat("\n🧩 Dataset:", dataset_name, "\n")
    
    map_dfr(method_names, function(method) {
      runs <- ds_list[[method]]
      if (is.null(runs)) return(NULL)
      cat("  ├─", method, ": ")
      
      # 提取每个方法的多次实验结果
      df <- map_dfr(seq_along(runs), function(i) {
        res <- runs[[i]]
        
        # --- 安全检查部分 ---
        if (is.null(res)) return(NULL)
        if (is.atomic(res) && all(is.na(res))) return(NULL)
        if (is.list(res) && length(res) == 0) return(NULL)
        if (!is.list(res)) return(NULL)
        
        # 提取 metric_cal 输出
        vals <- lapply(metrics, function(m) {
          if (!is.null(res[[m]]) && length(res[[m]]) == 1) {
            as.numeric(res[[m]])
          } else {
            NA_real_
          }
        })
        names(vals) <- metrics
        
        if (all(is.na(unlist(vals)))) return(NULL)
        
        data.frame(
          dataset = dataset_name,
          method = name_map[[method]] %||% method,
          run = i,
          as.data.frame(vals, optional = TRUE),
          check.names = FALSE
        )
      })
      
      cat(nrow(df), "valid runs\n")
      return(df)
    })
  }
  
  # --- 汇总所有数据集 ---
  metric_df <- bind_rows(map2(all_datasets, names(all_datasets), extract_one_dataset))
  mapping <- c(
    "accuracy" = "Accuracy",
    "weighted_f1" = "Weighted-F1"
  )
  colnames(metric_df) <-  ifelse(colnames(metric_df)  %in% names(mapping), mapping[colnames(metric_df) ], colnames(metric_df) )
  
  if (nrow(metric_df) == 0) {
    stop("❌ 没有可用的结果（全部 NA 或结构不匹配）")
  }
  
  # --- 转为长表 ---
  metric_long <- metric_df %>%
    mutate(method = factor(method, levels = method_order)) %>%
    pivot_longer(cols = all_of(ifelse(metrics %in% names(mapping), mapping[metrics], metrics)), names_to = "metric", values_to = "value")
  
  # --- 计算统计量 ---
  summary_tbl <- metric_long %>%
    group_by(dataset, method, metric) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd   = sd(value, na.rm = TRUE),
      n    = sum(!is.na(value)),
      se   = sd / sqrt(pmax(n, 1)),
      .groups = "drop"
    ) %>%
    arrange(metric, dataset, factor(method, levels = method_order))
  
  # --- 显著性检验 ---
  tests <- metric_long %>%
    group_by(dataset, metric) %>%
    rstatix::pairwise_wilcox_test(value ~ method, p.adjust.method = "BH") %>%
    ungroup()
  
  # --- 画图：箱线图 ---
  p_box <- ggplot(metric_long, aes(x = method, y = value, fill = method)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.18, size = 1.3, alpha = 0.6) +
    facet_grid(metric ~ dataset, scales = "free_y") +
    scale_x_discrete(limits = method_order) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none") +
    labs(title = "",
         x = NULL, y = "Score")
  
  # --- 画图：平均±标准误 ---
  p_bar <- ggplot(summary_tbl, aes(x = method, y = mean, fill = method)) +
    geom_col(alpha = 0.9) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    facet_grid(metric ~ dataset, scales = "free_y") +
    scale_x_discrete(limits = method_order) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none") +
    labs(title = "Average Performance (Mean ± SE)",
         x = NULL, y = "Mean ± SE")
  
  # --- 保存图像 ---
  if (!is.null(save_prefix)) {
    ggsave(paste0(save_prefix, "_boxplot.png"), p_box, width = 12, height = 6, dpi = 300)
    ggsave(paste0(save_prefix, "_barplot.png"), p_bar, width = 12, height = 6, dpi = 300)
  }
  
  # --- 输出汇总 ---
  cat("\n✅ Summary of valid runs per dataset-method:\n")
  print(summary_tbl %>%
          group_by(dataset, method) %>%
          summarise(n = sum(n, na.rm = TRUE), .groups = "drop"))
  
  invisible(list(
    metric_df   = metric_df,
    metric_long = metric_long,
    summary     = summary_tbl,
    tests       = tests,
    plots       = list(box = p_box, bar = p_bar)
  ))
}


  


# 假设你已经有每个数据集的 5 个方法结果：
# dataset_X 是一个 list，包含 result_scanmf/result_scCATCH/...，每个是长度10的list
out <- compare_methods(
  all_datasets = list(
    "Lawlor" = Lawlor,
    "Segerstolpe" = Segerstolpe,
    "Muraro" = Muraro,
    "Romanov"= Romanov,
    "Zeisel" =Zeisel
  ),
  metrics = c("accuracy", "weighted_f1"),
  save_prefix = "cross_pancreas"
)


# 看表
out$summary            # 各方法均值/方差/SE/有效重复数
out$tests              # 显著性检验结果
# 画图对象（可直接 print）
out$plots$box
out$plots$bar

all_datasets = list(
  "Lawlor" = Lawlor,
  "Segerstolpe" = Segerstolpe,
  "Muraro" = Muraro,
  "Romanov"= Romanov
)



save(all_datasets,file = 'in_multiple.RData')

load( 'in_multiple.RData')


out <- compare_methods(
  all_datasets,
  metrics = c("accuracy", "weighted_f1"),
  save_prefix = "cross_pancreas"
)


# 看表
out$summary            # 各方法均值/方差/SE/有效重复数
out$tests              # 显著性检验结果
# 画图对象（可直接 print）
out$plots$box
out$plots$bar


