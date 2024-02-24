# Install packages for Seurat Analysis
install.packages("Seurat")
install.packages("tidyverse")
reticulate::py_install(packages = 'umap-learn')
#install.packages("leidenAlg")
#install.packages('installr')
install.packages('Matrix')
install.packages("reticulate")
#install.Rtools()
install.packages('Matrix')

#Load libraries 
library(Seurat)
library(tidyverse)
library(leidenAlg)
library(Matrix)
library(installr)


#Load Data 
#Sometimes have to change the path from \ --> /
data_dir <- "C:/Users/orsmi/OneDrive/Desktop/CSB_185 Project_(R)/DropSeq.Cortex.7days.24hrs.TBI.Sham/Barcodes, features, matrix"
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
#Need to remove all parts of the name of the file that does not include just the barcode part 
data <- Read10X(data.dir = data_dir)
metadata <- read.table(file ="C:/Users/orsmi/OneDrive/Desktop/CSB_185 Project_(R)/DropSeq.Cortex.7days.24hrs.TBI.Sham/GSE180862_DropSeq.Cortex.7days.24hrs.TBI.Sham.metaData.tsv",sep='\t',header=TRUE)
yang.seurat.obj = CreateSeuratObject(counts = data, project="Frontal Cortex", min.cells = 3, min.features = 200)
str(yang.seurat.obj)

#Initializing Condition Column (Only upper case letters = 7 days (S or T) and Sham or TBI = 24 hours)
yang.seurat.obj$Condition <- NA # or any other initialization value

#Sham
yang.seurat.obj$Condition[yang.seurat.obj$orig.ident == "S1Cortex"] <- "Sham"
yang.seurat.obj$Condition[yang.seurat.obj$orig.ident == "S2Cortex"] <- "Sham"
yang.seurat.obj$Condition[yang.seurat.obj$orig.ident == "S3Cortex"] <- "Sham"
yang.seurat.obj$Condition[yang.seurat.obj$orig.ident == "Sham4Cortex"] <- "Sham"
yang.seurat.obj$Condition[yang.seurat.obj$orig.ident == "Sham6Cortex"] <- "Sham"
yang.seurat.obj$Condition[yang.seurat.obj$orig.ident == "Sham2Cortex"] <- "Sham"
#TBI
yang.seurat.obj$Condition[yang.seurat.obj$orig.ident == "T1Cortex"] <- "TBI"
yang.seurat.obj$Condition[yang.seurat.obj$orig.ident == "T2Cortex"] <- "TBI"
yang.seurat.obj$Condition[yang.seurat.obj$orig.ident == "T3Cortex"] <- "TBI"
yang.seurat.obj$Condition[yang.seurat.obj$orig.ident == "TBI1Cortex"] <- "TBI"
yang.seurat.obj$Condition[yang.seurat.obj$orig.ident == "TBI3Cortex"] <- "TBI"
yang.seurat.obj$Condition[yang.seurat.obj$orig.ident == "TBI5Cortex"] <- "TBI"

#Initializing Timepoint Column
yang.seurat.obj$Timepoint <- NA # or any other initialization value

#Sham
yang.seurat.obj$Timepoint[yang.seurat.obj$orig.ident == "S1Cortex"] <- "7days"
yang.seurat.obj$Timepoint[yang.seurat.obj$orig.ident == "S2Cortex"] <- "7days"
yang.seurat.obj$Timepoint[yang.seurat.obj$orig.ident == "S3Cortex"] <- "7days"
yang.seurat.obj$Timepoint[yang.seurat.obj$orig.ident == "T1Cortex"] <- "7days"
yang.seurat.obj$Timepoint[yang.seurat.obj$orig.ident == "T2Cortex"] <- "7days"
yang.seurat.obj$Timepoint[yang.seurat.obj$orig.ident == "T3Cortex"] <- "7days"

#TBI
yang.seurat.obj$Timepoint[yang.seurat.obj$orig.ident == "Sham4Cortex"] <- "24hrs"
yang.seurat.obj$Timepoint[yang.seurat.obj$orig.ident == "Sham6Cortex"] <- "24hrs"
yang.seurat.obj$Timepoint[yang.seurat.obj$orig.ident == "Sham2Cortex"] <- "24hrs"
yang.seurat.obj$Timepoint[yang.seurat.obj$orig.ident == "TBI1Cortex"] <- "24hrs"
yang.seurat.obj$Timepoint[yang.seurat.obj$orig.ident == "TBI5Cortex"] <- "24hrs"
yang.seurat.obj$Timepoint[yang.seurat.obj$orig.ident == "TBI3Cortex"] <- "24hrs"

# 1. QC -------
View(yang.seurat.obj@meta.data)
# % MT reads
yang.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(yang.seurat.obj, pattern = "^mt-")
View(yang.seurat.obj@meta.data)

VlnPlot(yang.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(yang.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering -----------------
yang.seurat.obj <- subset(yang.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                            percent.mt < 15)

# 3. Normalize data ----------
yang.seurat.obj <- NormalizeData(yang.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 4. Identify highly variable features --------------
yang.seurat.obj <- FindVariableFeatures(yang.seurat.obj, selection.method = "vst", nfeatures = 3000)

# 5. Scaling -------------
all.genes <- rownames(yang.seurat.obj)
yang.seurat.obj <- ScaleData(yang.seurat.obj, features = all.genes)

# 6. Perform Linear dimensionality reduction --------------
yang.seurat.obj <- RunPCA(yang.seurat.obj, features = VariableFeatures(object = yang.seurat.obj))
# determine dimensionality of the data
ElbowPlot(yang.seurat.obj)

# 7. Clustering ------------ #Might have to install leidenalg using pip install leidenalg
yang.seurat.obj <- FindNeighbors(yang.seurat.obj, dims = 1:20)
yang.seurat.obj <- FindClusters(yang.seurat.obj, resolution = 1, algorithm  = 4)

distance, total,displacement, tortuosity,vertical = [0]