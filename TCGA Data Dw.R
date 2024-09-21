options(warn=-1)

library(GDCRNATools)
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("MultiAssayExperiment")
library("maftools")
library("dplyr")
library("ComplexHeatmap")
library(GenomicRanges)
library(DESeq2)
library(ggplot2)

# get a list of projects
#gdcprojects <- getGDCprojects()
#getProjectSummary('TCGA-BRCA')

# build a query to retrieve gene expression data without barcode ------------
query_TCGA <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts', 
                       data.type = 'Gene Expression Quantification',
                       access = 'open')
getResults(query_TCGA)

# download data - GDCdownload

GDCdownload(query=query_TCGA)

#tcga_brca_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)

# prepare data
tcga_brca_data <- GDCprepare(query_TCGA)

colnames(colData(tcga_brca_data))
table(tcga_brca_data@colData$ajcc_pathologic_stage)

dim(assay(tcga_brca_data)) #Gene Expression Data
brca_gene_exp_matrix <- assay(tcga_brca_data)


# building a query for clinical data
query_TCGA <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Clinical',
                       data.type = 'Clinical Supplement')
clinical <- GDCquery_clinic("TCGA-BRCA")
output_clinical <- getResults(query_TCGA)


#tpm matrix
brca_matrix <- assay(tcga_brca_data, 'tpm_unstrand')


#DEG Analysis

clinical_data = colData(tcga_brca_data)
group = factor(clinical_data$definition)


group = relevel(group, ref="Solid Tissue Normal")
design = model.matrix(~group)
head(design)

dge = DGEList( # creating a DGEList object
  counts=assay(tcga_brca_data),
  samples=colData(tcga_brca_data),
  genes=as.data.frame(rowData(tcga_brca_data)))

# filtering
keep = filterByExpr(dge,design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) #  use rm() to remove objects from memory if you don't need them anymore

dge = calcNormFactors(dge,method="TMM")
v = voom(dge,design,plot=TRUE)

fit = lmFit(v, design)
fit = eBayes(fit)

topGenes = topTable(fit, coef=1, sort.by="p", number = 'INF')
print(topGenes)

write.csv(topGenes, file="topgenes.csv")



#DEG Analysis (Method 2)


# normalization of genes
dataNorm <- TCGAanalyze_Normalization(
  tabDF = tcga_brca_data, 
  geneInfo =  geneInfoHT
)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt),
  typesample = c("NT")
)

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt), 
  typesample = c("TP")
)

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,samplesNT],
  mat2 = dataFilt[,samplesTP],
  Cond1type = "Solid Tissue Normal",
  Cond2type = " Primary Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT",
  
)

dataDEGs$P.Value.adjusted <- p.adjust(dataDEGs$PValue, method = "BH")

# Select differentially expressed genes based on adjusted p-values
differentially_expressed_genes <- dataDEGs[dataDEGs$PValue.adjusted < 0.05, ]

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = dataDEGs,
  typeCond1 = "Primary Tumor",
  typeCond2 = "Solid Tissue Normal",
  TableCond1 = dataFilt[,samplesTP],
  TableCond2 = dataFilt[,samplesNT]
)

write.csv(dataDEGsFiltLevel, file="DEG_Samples.csv")



library(gplots)
library(readxl)
library(plyr)
library(pheatmap)
c <- read.csv("Heatmap_BRCA.csv")

library(ggplot2)
library(reshape2)
c <- melt(c, id.vars = "X")

ggplot(c, aes(x = LogFC, y = X)) +
  geom_tile(aes(fill = LogFC), colour = "white") +
  scale_fill_gradient(low = "cyan", high = "deeppink3") +
  theme_minimal()


Heatmap <- read.csv("Heatmap_BRCA.csv")
df <- data.matrix(Heatmap)

# Add classic arguments like main title and axis title
heatmap(df, Colv = NA, Rowv = NA, scale="column", col = coul, xlab="LogFC", ylab="X", main="heatmap")

heatmap.2(df, scale = "none", col = bluered(200), key.xlab = "logFC",
          trace = "none", density.info = "none",cexCol=1 #,cexRow = 1
)


df <- data.frame(
  GeneName = c("UCN3",
           "MUC2",
           "MYL2",
           "CGA",
           "CSAG1",
           "MAGEA12",
           "ACTL8",
           "MAGEA1",
           "IBSP",
           "MAGEA3",
           "KLHL1",
           "MYH2",
           "CKM",
           "MIR1-1HG",
           "MYH7",
           "PPP1R3A",
           "MYL1",
           "STRIT1",
           "C10orf71",
           "NRAP"),
  LogFC = c(8.90947631,
           8.670556009,
           -6.752772361,
           8.580617516,
           8.382966999,
           8.170402832,
           7.763106993,
           7.759037053,
           7.73295213,
           7.620344371,
           7.678597833,
          -7.263850008,           
          -7.081435584,
           -6.887413735,
           -6.611442989,
           -6.603212387,
           -6.490826284,
           -6.431899113,
           -6.086854318,
           -6.043186995
  )
)

heatmap.2(cbind(df$LogFC, df$LogFC), trace="n", Colv = NA, 
          dendrogram = "row", labCol = "", key.xlab = "LogFC", labRow = df$GeneName, cexRow = 0.75, col = col<- colorRampPalette(c("purple", "blue", "purple"))(256)
)



# Load the necessary libraries (if not already loaded)
library(WGCNA)

# Load your TPM data
# Replace 'brca_gene_exp_matrix' with your actual TPM matrix
# If you haven't already loaded it, you can uncomment the line below:
# brca_gene_exp_matrix <- assay(tcga_brca_data)

# Perform WGCNA

# 1. Data preprocessing (optional but recommended)
# Log2 transform and normalize your TPM data
tpm_data <- log2(brca_matrix + 1)

# 2. Create a sample trait data frame (containing sample IDs and potentially other clinical metadata)
# Replace 'sample_metadata' with your actual sample metadata if available
# You might need to adjust this step to match your metadata format
sample_metadata <- data.frame(sampleID = colnames(tpm_data))

# 3. Build a signed weighted adjacency matrix using a soft thresholding power
power <- 6  # You can adjust this value; it's often chosen through analysis
adjacency_matrix <- adjacency(datExpr = tpm_data, power = power)


# 4. Calculate the Topological Overlap Matrix (TOM)
TOM <- TOMsimilarity(adjacency_matrix)

# 5. Hierarchical clustering and module identification
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
modules <- cutreeStatic(geneTree, cutHeight = 0.25, minSize = 30)

# 6. Module-trait relationships (correlation between module eigengenes and clinical traits)
module_colors <- labels2colors(modules)

# Count the number of genes in each module
module_gene_counts <- table(module_colors)

# Display the number of genes in each module
print(module_gene_counts)


MEs <- moduleEigengenes(tpm_data, modules)


nGenes = ncol(tpm_data)
nSamples = nrow(tpm_data)
MEs0 = moduleEigengenes(tpm_data, module_colors)$eigengenes
MEs1 = orderMEs(MEs0)
moduleTraitCor = cor(MEs1, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = signif(moduleTraitPvalue, 1),
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# 7. Visualization
# Plot the dendrogram and module colors
pdf("Dendrogram_ModuleColors.pdf")
plotDendroAndColors(geneTree, module_colors)
dev.off()

moduleTraitHeatmap(modules)
# Plot module-trait relationships
pdf("ModuleTraitRelationships.pdf")
plotMEs(MEs, module_colors = module_colors, shape = "rect")
dev.off()



#method 2


# Load the necessary libraries (if not already loaded)
library(WGCNA)

# Load your TPM data
# Replace 'brca_gene_exp_matrix' with your actual TPM matrix
# If you haven't already loaded it, you can uncomment the line below:
# brca_gene_exp_matrix <- assay(tcga_brca_data)

# Perform WGCNA

# 1. Data preprocessing (optional but recommended)
# Log2 transform and normalize your TPM data
tpm_data <- log2(brca_gene_exp_matrix + 1)

# 2. Create a sample trait data frame (containing sample IDs and potentially other clinical metadata)
# Replace 'sample_metadata' with your actual sample metadata if available
# You might need to adjust this step to match your metadata format
sample_metadata <- data.frame(sampleID = colnames(tpm_data))

# 3. Build the network and modules using blockwiseModules

power <- 6  
network_data <- blockwiseModules(
  datExpr = tpm_data,
  traitData = sample_metadata,
  power = power,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  verbose = 0
)

module_colors <- network_data$colors

# Assuming wgcna_network is your WGCNA network object obtained from blockwiseModules








# Visualize module dendrogram and module colors
plotDendroAndColors(network_data$dendrograms[[1]], network_data$colors)




plotModuleTrait(network_data, 
                useColor = TRUE, 
                showCor = TRUE, 
                moduleOrder = "size", 
                traitOrder = "alphabetical", 
                geneGroupLabels = FALSE)

# Visualize module-trait relationships (heatmap)
moduleTraitHeatmap(network_data)
moduleTraitHeatmap(tpm_data, moduleColors = module_colors)


#method 3


# 4. Transpose module_colors and calculate correlations
correlations <- cor(t(module_colors), tpm_data, method = "spearman")

# 5. Create a heatmap to visualize the module-trait relationships
heatmap(correlations, 
        Colv = NA,  # Disable clustering of columns (traits)
        scale = "none",  # Do not scale the data
        col = colorRampPalette(c("blue", "white", "red"))(50),  # Color scheme
        main = "Module-Trait Relationships Heatmap")


#method  4

library(impute)


# Example: Log2 transformation
tpm_data_log2 <- log2(brca_matrix + 1)
# Create a matrix of gene correlations
gene_correlation_matrix <- cor(tpm_data_log2)

# Manually determine the soft threshold using the scale-free topology criterion
sft <- pickSoftThreshold(datExpr = gene_correlation_matrix, powerVector = c(1:20), verbose = 5)

# Select the power for which the scale-free topology fit curve starts to plateau
soft_threshold <- sft$powerEstimate

# Check if the soft threshold is within the valid range
if (soft_threshold < 1 | soft_threshold > 30) {
  stop("Soft threshold must be between 1 and 30.")
}

# Create a WGCNA network using the selected soft threshold
modules <- blockwiseModules(gene_correlation_matrix, power = soft_threshold,
                            TOMType = "unsigned", minModuleSize = 30,
                            mergeCutHeight = 0.25, reassignThreshold = 0.04)

# Access the adjacency matrix from the network
adjacency_matrix <- modules$TOM

# Plot the dendrogram and module colors
plotDendroAndColors(modules$dendrograms[[1]], modules$colors)

# Visualize the module-trait relationships
moduleTraitHeatmap(modules)

# Perform functional enrichment analysis
# Example: using the `enrichR` package
# Install 'enrichR' if not already installed: install.packages("enrichR")
library(enrichR)
enrich_results <- enrichR(modules$MEs$moduleColors, organism = "human")


#Method 5

# Assuming brca_gene_exp_matrix is your gene expression data
# Assuming samples are in columns and genes are in rows

# Log-transform the data (if it's not in log scale)
logged_data <- log1p(brca_gene_exp_matrix)


# Create a numeric matrix
numeric_data <- as.data.frame(t(logged_data))

# Check if your data is correctly structured
str(logged_data)

# Soft Thresholding
soft_threshold <- pickSoftThreshold(numeric_data, powerVector = 1:10, networkType = "unsigned")

# Build the network using the selected soft threshold
adjacency_matrix <- adjacency(datExpr = numeric_data, power = soft_threshold$powerEstimate, type = "unsigned")

# Construct the Topological Overlap Matrix (TOM)
TOM <- TOMsimilarity(adjacency_matrix)

# Hierarchical clustering of the TOM
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Module detection using the Dynamic Tree Cut algorithm
modules <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)


# Visualization of the network and modules
plotDendroAndColors(dendro = geneTree, colors = modules)

# Calculate module eigengenes
MEs <- moduleEigengenes(numeric_data, modules)

# Calculate module-trait relationships
module_trait_cor <- moduleTraitRelationship(moduleColors = modules, traits = trait_data)

# Plot the heatmap of module-trait relationships
heatmap(module_trait_cor$correlationMatrix,
        Colv = NA, Rowv = NA,
        col = colorRampPalette(c("blue", "white", "red"))(50),
        main = "Module-Trait Relationships")

