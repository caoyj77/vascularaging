library(readr)
library(readxl)
library(dplyr)
library(tidyverse)
# 读取输入蛋白列表
protein_list <- read_xlsx("/Users/")
protein_list <- list("RBP5")
your_proteins<-protein_list$angina
your_proteins_df <- data.frame(GeneSymbol = your_proteins)
your_proteins_df <- data.frame(GeneSymbol = protein_list)
# ==== 3. 读取DGIdb提供的 interactions.tsv 文件 ====
dgidb_data <- read_tsv("/Users/", comment = "#")  # 跳过注释行
#gidb_interactions <- read_tsv("/Users/")
library(data.table)
dgidb_data <- fread("/Users/")
glimpse(dgidb_data)

# ==== 4. 提取所有出现在DGIdb中的基因名 ====
dgidb_genes <- unique(dgidb_data$gene_name)

# ==== 5. 标注你的蛋白是否为已知成药靶点 ====
annotated <- your_proteins_df %>%
  mutate(DGIdb_Druggable = GeneSymbol %in% dgidb_genes)

# 筛选 gene_name 为 "RBP5" 的所有药物交互
rbp5_targets <- dgidb_data %>%
  filter(gene_name == "RBP5")
# ==== 6. 可选：添加具体药物信息（合并药物名） ====
# 整理 drug-target 信息（每个基因对应的药物名合并为一列）
dgidb_grouped <- dgidb_data %>%
  select(gene_name, drug_name) %>%
  group_by(gene_name) %>%
  summarise(Drugs = paste(unique(drug_name), collapse = "; "), .groups = "drop")

# 合并药物信息到你的数据
final_result <- annotated %>%
  left_join(dgidb_grouped, by = c("GeneSymbol" = "gene_name"))

# ==== 7. 导出结果 ====
write_csv(final_result, "/Users/")

# ==== 8. 简单统计：有多少是成药靶点 ====
table(final_result$DGIdb_Druggable)

