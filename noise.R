library(tidyverse)

# ---- 1. 原始数据（你已有） ----
df <- tribble(
  ~label_ret, ~type,         ~noise, ~mean,   ~sd,
  
  # 10% label retained
  "10% label retained", "Label noise",   0,  0.9972, 0.0036,
  "10% label retained", "Label noise",   20, 0.9875, 0.0124,
  "10% label retained", "Label noise",   40, 0.9732, 0.0324,
  "10% label retained", "Label noise",   60, 0.9062, 0.0950,
  "10% label retained", "Marker noise",  20, 0.9882, 0.0237,
  "10% label retained", "Marker noise",  40, 0.9847, 0.0211,
  "10% label retained", "Marker noise",  60, 0.9740, 0.0321,
  
  # 20% label retained
  "20% label retained", "Label noise",   0,  0.9981, 0.0024,
  "20% label retained", "Label noise",   20, 0.9862, 0.0116,
  "20% label retained", "Label noise",   40, 0.9680, 0.0443,
  "20% label retained", "Label noise",   60, 0.9592, 0.0639,
  "20% label retained", "Marker noise",  20, 0.9894, 0.0274,
  "20% label retained", "Marker noise",  40, 0.9793, 0.0234,
  "20% label retained", "Marker noise",  60, 0.9617, 0.0358,
  
  # 30% label retained
  "30% label retained", "Label noise",   0,  0.9965, 0.0057,
  "30% label retained", "Label noise",   20, 0.9905, 0.0062,
  "30% label retained", "Label noise",   40, 0.9577, 0.0476,
  "30% label retained", "Label noise",   60, 0.8917, 0.1141,
  "30% label retained", "Marker noise",  20, 0.9935, 0.0120,
  "30% label retained", "Marker noise",  40, 0.9787, 0.0397,
  "30% label retained", "Marker noise",  60, 0.9451, 0.0284
)

# ---- 2. 使用 SE 替代 SD ----
df <- df %>% mutate(se = sd / sqrt(20))

# ---- 3. 绘图 ----
p <- ggplot(df, aes(x = noise, y = mean, color = type)) +
  geom_line(aes(linetype = type), linewidth = 0.7) +   # PDF 中线条更细更好看
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = pmax(mean - se, 0),
      ymax = pmin(mean + se, 1)
    ),
    width = 4,
    linewidth = 0.6
  ) +
  facet_wrap(~ label_ret, nrow = 1) +
  scale_color_manual(values = c(
    "Label noise"  = "#1f78b4",
    "Marker noise" = "#ff7f00"
  )) +
  scale_linetype_manual(values = c(
    "Label noise"  = "solid",
    "Marker noise" = "dashed"
  )) +
  scale_x_continuous(breaks = c(0,20,40,60)) +
  scale_y_continuous(limits = c(0.85, 1.0)) +
  labs(
    x = "Noise level (%)",
    y = "Accuracy",
    color = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 13) +     # 全局字体全部设为 13
  theme(
    strip.background = element_rect(fill = "white", color = "black"), # 保留白框标题
    strip.text = element_text(size = 13, face = "bold"),
    
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    
    legend.position = "bottom",
    legend.text = element_text(size = 13),
    legend.title = element_blank(),
    legend.key.width = unit(1.2, "cm"),
    
    panel.grid.minor = element_blank()
  )


p
# 输出到 PDF（矢量图）
pdf("noise_plot.pdf", width = 7, height = 4, useDingbats = FALSE)

print(p)

dev.off()
