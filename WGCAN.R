library(WGCNA)

# 假设protein_data是一个包含蛋白表达数据的矩阵（行是样本，列是蛋白）
external_df <- read.csv("/Users/")
significant_markers <- external_df [external_df$pval_bfi < 0.05, "protein_code"] %>% unlist()
protein_expr <- read.csv("/Users/", header = TRUE, row.names = 1)
valid_markers <- significant_markers[significant_markers %in% colnames(protein_expr)]
protein_exp <- protein_expr[, valid_markers]
print(significant_markers)
# 预处理数据
protein_data <- t(scale(t(protein_exp))) # 转置后标准化

# 计算相关性矩阵
cor_matrix <- cor(protein_data, method = "pearson")

# 选择软阈值
power <- pickSoftThreshold(protein_data, powerVector = c(1:20), verbose = 5)$powerEstimate

# 构建加权邻接矩阵
adjacency_matrix <- adjacency(cor_matrix, power = power)

# 计算拓扑重叠矩阵
TOM_matrix <- TOMsimilarity(adjacency_matrix)

# 绘制拓扑重叠矩阵热图
TOMplot(dissTOM, geneTree, as.character(dynamicColors), main = "Topological Overlap Matrix")

# 聚类分析
dissTOM <- 1 - TOM_matrix
geneTree <- hclust(as.dist(dissTOM), method = "average")

# 进行层次聚类
geneTree = hclust(as.dist(dissTOM), method = "average")

# 设定模块检测的参数
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

# 将模块编号转换为颜色
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors) 
# 绘制拓扑重叠矩阵热图
TOMplot(dissTOM, geneTree, as.character(dynamicColors), main = "Topological Overlap Matrix")

# 划分模块
modules <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE)

# 查看每个蛋白所属的模块
table(modules)

threshold <- 0.1
edges <- which(TOM_matrix > threshold, arr.ind = TRUE)
library(readxl)
protein_names <- read_excel("/Users/")
protein_names <- as.vector(protein_names)
protein_name <- protein_names[["pn"]]

# 获取相关的蛋白对
edge_list <- data.frame(from = protein_name[edges[, 1]], to = protein_name[edges[, 2]], weight = TOM_matrix[edges])

# 保存边数据为csv格式
write.csv(edge_list, "/Users/", row.names = FALSE)


# 假设modules是模块分配结果，protein_names是蛋白名称
node_info <- data.frame(protein = protein_names, module = modules)

# 检查数据框内容
head(node_info)

# 保存节点信息为csv格式
write.csv(node_info, "/Users/", row.names = FALSE)


#####
# 加载WGCNA包
library(WGCNA)
options(stringsAsFactors = FALSE)

# 检查数据是否合格（行是样本，列是基因）
datExpr <- protein_data  # 替换成你实际的数据变量

# 1. 自动选择软阈值 power（或你已选定）
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 2. 构建网络并识别模块
net <- blockwiseModules(
  datExpr,
  power = sft$powerEstimate, # 或自定义数值
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.1,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "TOM",
  verbose = 3
)

# 3. 提取模块颜色标签
moduleColors <- labels2colors(net$colors)
table(moduleColors)  # 看看是不是大部分基因都被归到一个模块了

# 4. 计算模块特征基因
MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes

# 5. 计算kME（每个基因与模块特征基因的相关性）
kMEtable <- signedKME(datExpr, MEs)

# 6. 选择目标模块（例如：blue模块）
targetModule <- "blue"
moduleGenes <- moduleColors == targetModule
kME_target <- kMEtable[moduleGenes, paste0("kME", targetModule)]

# 7. 找出模块内kME最大的基因（hub gene）
hubGene <- rownames(kMEtable[moduleGenes, ])[which.max(kME_target)]
cat("Hub gene in module", targetModule, "is:", hubGene, "\n")
