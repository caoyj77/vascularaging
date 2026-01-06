rm(list=ls())  
library(data.table)  
library(dplyr) 
cov_df <- read.csv('/Users/')
nmr_df_filtered  <- read.csv('/Users/')
label_df <- read.csv('/Users/')
column_names <- colnames(nmr_df)
metabolite_list <- column_names[3:length(column_names)]
metabolite_list <- nmr_df$Label
nmr_df_filtered <- nmr_df[, c('eid', metabolite_list)]

metabolite_data <- merge(cov_df, nmr_df_filtered, by = "eid")
print(colnames(metabolite_data))
colnames(metabolite_data)[3] <- "Age"

#消除协变量影响
# 确保使用 lm() 时协变量没有 NA 值

for (metabolite in metabolite_list) {
  formula <- as.formula(paste0(metabolite, " ~ Smoke + SBP + FastingTime + BL_HYPT"))
  metabolite_data[!is.na(metabolite_data[, metabolite]) & !is.na(metabolite_data$Smoke) &
                    !is.na(metabolite_data$SBP) & !is.na(metabolite_data$FastingTime) &
                    !is.na(metabolite_data$BL_HYPT), metabolite] <- resid(lm(formula, data = metabolite_data))
}

# for (metabolite in metabolite_list) {
#   formula <- as.formula(paste0(metabolite, " ~ Smoke + SBP + FastingTime"))
#   metabolite_data[!is.na(metabolite_data[, metabolite]) & !is.na(metabolite_data$Smoke) &
#                     !is.na(metabolite_data$SBP) & !is.na(metabolite_data$FastingTime) , metabolite] <- resid(lm(formula, data = metabolite_data))
# }

metabolite_data <- metabolite_data[ , -c(2:12)]

#处理每个疾病的匹配结果并计算z score
result <- data.frame(Metabolite = metabolite_list)
print(result)
bins <- seq(0,15, by = 2)  # 设置年份分组：从0到15，步长为1年
print(bins)
bins_level <- paste0(bins[-length(bins)], "-", bins[-1])  # 定义每个区间的标签，如 "0-1", "1-2", "2-3", ...
# for (bin in bins_level) {
#   result[[bin]] <- NA
# }

disease_columns <- list(
  c('SAH', 'SAH_targe', 'SAH_yrs'),
  c('ICH', 'ICH_targe', 'ICH_yrs'),
  c('CI', 'CI_trage', 'CI_yrs'),
  c('angina', 'angina_targe', 'angina_yrs'),
  c('AS', 'AS_targe', 'AS_yrs'),
  c('PAD', 'PAD_targe', 'PAD_yrs'),
  c('stroke', 'stroke_targe', 'stroke_yrs'),
  c('hf', 'hf_trage', 'hf_yrs'),
  c('CAD', 'CAD_trage', 'CAD_yrs'),
  c('MI', 'MI_targe', 'MI_yrs')
  # c('stroke','target_y','BL2Target_yrs')
)

# 循环遍历每个疾病列
for (disease_info in disease_columns) {
  disease <- disease_info[1]
  target_col <- disease_info[2]
  years_col <- disease_info[3]
  
  disease_result <- result
  matched_data <- fread(paste0("/Users/", disease, "_matched_result.txt"), sep = "\t", header = T)
  # 检查是否有缺失的列名，若有则用默认名称替换
  if(any(is.na(names(matched_data)))){
   na_idx <- which(is.na(names(matched_data)))
   names(matched_data)[na_idx] <- paste0("V", na_idx)
 }
  # 检查是否有缺失的列名或空字符串列名，若有则用默认名称替换
  commented_code <- function() {col_names <- names(matched_data)
  for (i in seq_along(col_names)) {
    if (is.na(col_names[i]) || col_names[i] == "") {
      new_name <- paste0("V", i)
      # 检查新列名是否重复
      while (new_name %in% col_names) {
        i <- i + 1
        new_name <- paste0("V", i)
      }
      col_names[i] <- new_name
    }
  }
  names(matched_data) <- col_names}
  
  print(names(matched_data))

  
  matched_data <- matched_data %>%
    group_by(subclass) %>%
    mutate(BL2Target_yrs_for_1 = ifelse(target_y == 1, get(years_col), 0)) %>% #如果target_y == 1，BL2Target_yrs_for_1=BL2Target_yrs，else=0
    mutate(BL2Target_yrs = sum(BL2Target_yrs_for_1)) %>%
    ungroup() %>%
    select(-BL2Target_yrs_for_1) %>%
    bind_rows()
  matched_data <- as.data.frame(matched_data)
  # print(matched_data)
  
  matched_data <- matched_data[matched_data$BL2Target_yrs <=15 & matched_data$BL2Target_yrs >= 0, ]
  matched_data$BL2Target_yrs <- matched_data$BL2Target_yrs
  matched_data$bin <- cut(matched_data$BL2Target_yrs, bins)
  matched_data <- matched_data[!is.na(matched_data$bin), ]


  for (current_bin in unique(matched_data$bin)) {
    data <- matched_data[matched_data$bin == current_bin, ]
    if (nrow(data) == 0) {
      print(paste("Skipping empty bin:", current_bin))
      next
    }
    # print(data)
    case_ids <- data[data$target_y == 1,]$eid
    control_ids <- data[data$target_y == 0,]$eid
    if (length(case_ids) == 0 || length(control_ids) == 0) {
      print(paste("Skipping bin due to lack of cases or controls:", current_bin))
      next
    }
    case_metabolite_data <- metabolite_data[metabolite_data$eid %in% case_ids,]
   # case_metabolite_means <- colMeans(case_metabolite_data, na.rm = TRUE)
    control_metabolite_data <- metabolite_data[metabolite_data$eid %in% control_ids,][,-1]
    control_metabolite_means <- colMeans(control_metabolite_data, na.rm = TRUE)
    control_metabolite_sd <- apply(control_metabolite_data, 2, sd, na.rm = TRUE)
    control_metabolite_sd[is.na(control_metabolite_sd)] <- 1e-10
    if (any(control_metabolite_sd == 0, na.rm=TRUE)) {
      print(paste("Warning: Zero standard deviation for some metabolites in bin:", current_bin))
      control_metabolite_sd[control_metabolite_sd == 0] <- 1e-10  # 或者选择跳过这些代谢物
      
    }
    case_control_meta_zscore <- (case_metabolite_means - control_metabolite_means) / control_metabolite_sd
    disease_result[[current_bin]] <- case_control_meta_zscore
    # disease_result[[current_bin]] <- control_metabolite_means
    # print(disease_result)
    
  }

  # writepath <- paste0("/Users/", disease, "_metabolite_matrix_result.txt")
  # fwrite(matched_data, writepath, sep = " ", row.names = F, col.names = T, quote = F, na = "NA")
  writepath <- paste0("/Users/", disease, "_metabolite_zscore_result.txt")
  fwrite(disease_result, writepath, sep = " ", row.names = F, col.names = T, quote = F, na = "NA")

}

print(bins_level)
