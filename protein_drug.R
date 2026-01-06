library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)

# 1. 读取原始Excel数据（仅改路径，其他不动）
data <- read_excel(
  path = "/Users/",  # 替换为你的文件路径
  sheet = 1,
  col_names = TRUE
)

# 2. 直接构建“蛋白-药物”关联矩阵（不筛选，保留所有数据）
# 去重：避免同一“蛋白-药物”组合重复显示
heatmap_raw <- data %>%
  # 1. 仅保留“有效药物”记录（过滤无关联药物）
  filter(Drug_Type == "有效药物") %>%  # 直接用原始列名Drug_Type
  # 2. 按“疾病-基因”分组去重（避免同一组合的多个药物重复标记）
  group_by(GeneSymbol, Drug) %>%    # 直接用原始列名disease和GeneSymbol
  summarise(
    has_assoc = 1,  # 有有效药物关联→标记为1
    .groups = "drop"  # 取消分组，避免警告
  )


# 3. 生成完整网格（包含所有“蛋白-药物”组合，无关联为NA）
all_proteins <- unique(data$GeneSymbol)  # Y轴：所有蛋白
all_drugs <- unique(data$Drug)          # X轴：所有药物
heatmap_grid <- expand.grid(
  GeneSymbol = all_proteins,
  Drug = all_drugs,
  stringsAsFactors = FALSE
)

# 4. 合并关联标记（有=1，无=NA）
heatmap_data <- heatmap_grid %>%
  left_join(heatmap_raw, by = c("GeneSymbol", "Drug")) %>%
  mutate(
    assoc_status = ifelse(is.na(has_assoc), "No Link", "Linked")
  )

# 5. 绘制热图（Y=蛋白，X=药物）
ggplot(heatmap_data, aes(x = GeneSymbol, y = Drug)) +
  # 热图主体
  geom_tile(aes(fill = assoc_status), color = "gray80", size = 0.1) +
  # 颜色配置
  scale_fill_manual(values = c("No Link" = "white", "Linked" = "#2E86AB")) +
  # 标题与轴标签
  labs(
    title = "Protein-Drug Link Heatmap",
    x = "Drug",
    y = "Protein (GeneSymbol)",
    fill = "Link Status"
  ) +
  # 主题优化（避免标签重叠）
  theme_minimal() +
  theme(
    # X轴药物名旋转，缩小字号
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    # Y轴蛋白名缩小字号
    axis.text.y = element_text(size = 7),
    # 标题居中
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    # 移除网格线
    panel.grid = element_blank(),
    # 图例位置
    legend.position = "right"
  ) +
  # 热图无间隙
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

# 6. 保存热图（替换为你的保存路径）
ggsave(
  filename = "/Users/",
  width = ifelse(length(all_drugs) > 20, 16, 12),
  height = ifelse(length(all_proteins) > 20, 15, 10),
  dpi = 300,
  device = "pdf"
)

cat("热图已保存！")