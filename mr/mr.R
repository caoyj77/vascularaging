# 安装所需的 R 包（如果尚未安装）
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")

# 加载 R 包
library(TwoSampleMR)
library(dplyr)
library(data.table)
library(vroom)
library(MRPRESSO)
library(ieugwasr)

# 读取 UKB PPP 蛋白 GWAS 数据 1000kb、0.01、5e-8
clump_data<-vroom('/Users/')
clump_other<-vroom('/Users/')
exposure_dat <- vroom('/Users/') 
cluster <- read.csv('/Users/')
decode <- read.csv('/Users/')
head(ukb_prot)
colnames(clump_data)
exposure_dat <- exposure_dat[exposure_dat$exposure %in% cluster$pn, ]

# 读取 FinnGen 中风 GWAS 数据
finn_stroke <- vroom('/Users/') 
head(finn_stroke)

exposure_dat <- format_data(
  dat = clump_data,
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  phenotype_col = "gene",
  type = "exposure"
)

snp_clump <- ld_clump_local(
  dplyr::tibble(rsid = decode$SNP, pval = decode$pval.exposure),  # 暴露数据的 SNP 和 p 值
  clump_kb = 1000,   # clumping 窗口大小（kb）
  clump_r2 = 0.01,   # clumping 的 r² 阈值
  clump_p = 1e-6,    # exposure已经筛选过了
  #bfile = "/Users//data_maf0.01_rs_ref",
  bfile = "/Users/mr/1kg.v3/EUR",  # LD 参考面板路径
  plink_bin = "/Users/plink2"  # PLINK 可执行文件路径
)

output_base <- "/var/"
# 尝试读取 .clumps 文件
clumps_file <- paste0(output_base, ".clumps")
if (file.exists(clumps_file)) {
  # 读取 .clumps 文件内容
  decode_data <- read.table(clumps_file, header = FALSE)
  # 使用 clumps_data 进行后续分析
  #print(clumps_data)
} else {
  stop("PLINK 未生成 .clumps 文件")
}
colnames(clumps_data) <- c("1", "2", "rsid","4","5","6","7","8","9","10","11")
# 假设 snp_clump 中有一个 'rsid' 列存储 SNP ID
exposure_dat <- exposure_dat[exposure_dat$SNP %in% clumps_data$rsid, ]
cluster <- read.csv('/Users/')
exposure_dat <- exposure_dat[exposure_dat$exposure %in% cluster$pn, ]
write.csv(exposure_dat, file = "/Users/", row.names = FALSE)

exposure_dat <- format_data(
  dat = clump_data,
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  phenotype_col = "gene",
  type = "exposure"
)

outcome_dat <- format_data(
  dat = finn_stroke,
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "P",
  type = "outcome"
)
rm(column_names)
harmonized_data <- harmonise_data(decode, outcome_dat)
harmonized_data <- harmonise_data(decode, outcome_dat)
write.csv(harmonized_data, file = "/Users/", row.names = FALSE)

mr_results <- mr(harmonized_data)
write.csv(mr_results, file = "/Users/", row.names = FALSE)

dat.mr_heterogeneity<-mr_heterogeneity(harmonized_data)   
write.csv(dat.mr_heterogeneity, file = "/Users/", row.names = FALSE)
dat.mr_pleiotropy_test<-mr_pleiotropy_test(harmonized_data) 
write.csv(dat.mr_pleiotropy_test, file = "/Users/", row.names = FALSE)


loo <- mr_leaveoneout(harmonized_data)
print(mr_results)

# MR-PRESSO 检测离群值
presso <- mr_presso(
  BetaOutcome = "beta.outcome", 
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome", 
  SdExposure = "se.exposure",
  data = harmonized_data,
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  NbDistribution = 10000000
)

# 离群值校正后的结果
print(presso$`MR-PRESSO results`$`Outlier-corrected`)
# Step 6: 可视化结果
# 1. 散点图
mr_scatter_plot(mr_results, harmonized_data )

# 2. 森林图
mr_forest_plot(mr_singlesnp(harmonized_data ))

# 3. 漏斗图
mr_funnel_plot(mr_singlesnp(harmonized_data ))

###
input = merge(clump_data,outcome_dat,by="SNP")
input = input[input['gene'] =='MMP12', ]
input <- input[!duplicated(input$SNP), ]
maf_threshold <- 0.05
input <- input[input$eaf > maf_threshold, ]
input$maf <- pmin(input$eaf, 1 - input$eaf)
input <- input[input$eaf.exposure > maf_threshold, ]
input$maf <- pmin(input$eaf.exposure, 1 - input$eaf.exposure)
library(coloc)
exposure_coloc <- list(snp = input$SNP,pvalues = input$pval.exposure, N = input$samplesize.exposure, type = "quant", beta = input$beta.exposure, varbeta = input$se.exposure^2,MAF = input$maf)
exposure_coloc <- list(snp = input$SNP,pvalues = input$pval, N = input$samplesize, type = "quant", beta = input$beta, varbeta = input$se^2,MAF = input$maf)
outcome_coloc <- list(snp = input$SNP,pvalues = input$pval.outcome, N = 1832, type = "cc", beta = input$beta.outcome, varbeta = input$se.outcome^2, MAF = input$maf)
exposure_coloc <- list(snp = input$SNP,pvalues = input$pval.exposure, N = input$samplesize.exposure, type = "quant", beta = input$beta.exposure, varbeta = input$se.exposure^2)
outcome_coloc <- list(snp = input$SNP,pvalues = input$pval.outcome, N = 1832, type = "cc", beta = input$beta.outcome, varbeta = input$se.outcome^2)
coloc_results <- coloc.abf(dataset1 = exposure_coloc, dataset2 = outcome_coloc)
print(coloc_results$summary)
library(locuscomparer)
gwaslist <- list()
result <- list()
if(coloc_results$summary[6]>0.01){
  input$H4 = coloc_results$summary[6]
  result[[i]] = input
  gwas = cbind(input$SNP,input$pval)
  colnames(gwas) = c("rsid","pval")
  eqtl = cbind(input$SNP,input$pval.outcome)
  colnames(eqtl) = c("rsid","pval")
  write.table(gwas,"gwas.tsv",col.names = T,row.names = F,sep="\t",quote = F)
  write.table(eqtl,"eqtl.tsv",col.names = T,row.names = F,sep="\t",quote = F)
  gwaslist[[i]] = input$gene
  locuscompare("eqtl.tsv","gwas.tsv",legend =FALSE,snp = "rs1007350")  
  ggsave(paste(i,"locus.png"),width=10,height = 5)
}
