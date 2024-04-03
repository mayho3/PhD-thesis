# source("/data/mayerlab/cmayer/Project_CA/GitLab/cloneseq/lib5.R")
library(readr)
library(purrr)
library(UpSetR)
library(monocle3)
library(patchwork)
library(VennDiagram)

options(echo=TRUE)

dir.create("/home/mpg08/ho.may/crispr", showWarnings = FALSE)
setwd("/home/mpg08/ho.may/crispr") #set working directory 
####### Load LGE Data 3


load("/datastore_share/bioinformatics/datasets/seuratobjects/MUC28072_GEpb4b.Rdata")
DefaultAssay(E16) <- "RNA"
Idents(E16) <- "CCA_Assignment"
E16_subset_inhib <- subset(E16, idents = unique(E16@meta.data$CCA_Assignment)[-c(3,16,18)])
E16_subset_inhib <- subset(x = E16_subset_inhib, cells = WhichCells(E16_subset_inhib, expression = Neurod2 > 2 | Neurod6 > 2),invert = T)
Neurod2 >2 | Neurod6 >2), invert = T)
Idents(E16_subset_inhib) <- droplevels(Idents(E16_subset_inhib))
# MUC28072 <- subset(x = MUC28072, 

#        idents = c("Col4a1","e_Neurod2/6/Satb2","e_Neurod2/Epha5","e_Neurod6/Nfix","e_Nfix/Tcf4/Lhx2","e_Ppp2r2b/Nfib"), 
#        invert = TRUE)
# MUC28072 <- subset(x = MUC28072, cells = WhichCells(MUC28072, expression = Neurod2 > 2 | Neurod6 > 2),invert = T)
# Idents(MUC28072) <- droplevels(Idents(MUC28072))
# 
# load("/datastore_share/bioinformatics/datasets/seuratobjects/CA303.Rdata")
# load("/datastore_share/bioinformatics/datasets/seuratobjects/CA300.Rdata")
# load("/datastore_share/bioinformatics/datasets/seuratobjects/CA302.Rdata")
# load("/datastore_share/bioinformatics/datasets/seuratobjects/CA299.Rdata")
# load("/datastore_share/bioinformatics/datasets/seuratobjects/CA301.Rdata")
# load("/datastore_share/bioinformatics/datasets/seuratobjects/CA298.Rdata")
# 
CA303@meta.data[,"dataset"] <- "CA303"
DefaultAssay(CA303) <- "RNA"
CA303 <- subset(x = CA303, cells = WhichCells(CA303, expression = Neurod2 > 2 | Neurod6 > 2),invert = T)
# 
CA300@meta.data[,"dataset"] <- "CA300"
DefaultAssay(CA300) <- "RNA"
CA300 <- subset(x = CA300, cells = WhichCells(CA300, expression = Neurod2 > 2 | Neurod6 > 2),invert = T)
# 
CA302@meta.data[,"dataset"] <- "CA302"
DefaultAssay(CA302) <- "RNA"
CA302 <- subset(x = CA302, cells = WhichCells(CA302, expression = Neurod2 > 2 | Neurod6 > 2),invert = T)
# 
CA299@meta.data[,"dataset"] <- "CA299"
DefaultAssay(CA299) <- "RNA"
CA299 <- subset(x = CA299, cells = WhichCells(CA299, expression = Neurod2 > 2 | Neurod6 > 2),invert = T)
# 
CA301@meta.data[,"dataset"] <- "CA301"
DefaultAssay(CA301) <- "RNA"
CA301 <- subset(x = CA301, cells = WhichCells(CA301, expression = Neurod2 > 2 | Neurod6 > 2),invert = T)
# 
CA298@meta.data[,"dataset"] <- "CA298"
DefaultAssay(CA298) <- "RNA"
CA298 <- subset(x = CA298, cells = WhichCells(CA298, expression = Neurod2 > 2 | Neurod6 > 2),invert = T)
# 
# 
# #Integration
data.list <- c(E16_subset_inhib,CA303,CA300,CA302,CA299,CA301,CA298)
# for (i in seq_along(data.list)) {data.list[[i]] <- SCTransform(data.list[[i]],vars.to.regress = c("nCount_RNA","nFeature_RNA"))
# }
features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = features)
Eminence.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT",
                                           anchor.features = features)
Eminence.combined.sct <- IntegrateData(anchorset = Eminence.anchors, normalization.method = "SCT")
Eminence.combined.sct <- FindVariableFeatures(Eminence.combined.sct, assay = "integrated")
Eminence.combined.sct <- RunPCA(Eminence.combined.sct, verbose = FALSE, assay = "integrated")
Eminence.combined.sct <- RunUMAP(Eminence.combined.sct, dims = 1:15, verbose = FALSE, assay = "integrated")
Eminence.combined.sct <- FindNeighbors(Eminence.combined.sct, dims = 1:15)
Eminence.combined.sct <- FindClusters(Eminence.combined.sct, resolution = 0.8)
save(Eminence.combined.sct, file = "/home/mpg08/ho.may/crispr/E16.GE.combined.sct.Rdata")

load("/datastore_share/Users/cmayer/temp/Eminence.combined.sct.GEpb6.Rdata")

Eminence.combined.sct@meta.data$neuron_class_all <- plyr::mapvalues(E16_GE_assign_240222$CCA_Assignment_2, from = E16_GE_assign_240222$Cluster, to = E16_GE_assign_240222$class)

Eminence.combined.sct@meta.data[, "GE"] <- plyr::mapvalues(Eminence.combined.sct@meta.data$dataset,
                                                           from = c("CA303","CA300","CA302","CA299","CA301","CA298"),
                                                           to = c("LGE","LGE","CGE","CGE","MGE","MGE"
                                                           ))

DimPlot(Eminence.combined.sct, group.by = "GE")

FeaturePlot(Eminence.combined.sct, c("Adora2a","Tac1","Drd2"))

FeaturePlot(Eminence.combined.sct, c("Nes"))
ggsave("FeaturePlot_sct.pdf", width = 300, height = 300, units = "mm",scale = 1,dpi=600)
DimPlot(Eminence.combined.sct, group.by = "seurat_clusters", label = T)
ggsave("DimPlot_ident.pdf", width = 300, height = 300, units = "mm",scale = 1,dpi=600)
DimPlot(Eminence.combined.sct, group.by = "dataset", label = T, )
ggsave("DimPlot_dataset.pdf", width = 300, height = 300, units = "mm",scale = 1,dpi=600)

Idents(object = Eminence.combined.sct) <- "seurat_clusters"
marker.Eminence.combined.sct <- FindAllMarkers(Eminence.combined.sct, only.pos = T,logfc.threshold = 0.98,min.pct =0.7)
write_csv(marker.Eminence.combined.sct, path = "/home/mpg08/ho.may/crispr/marker.Eminence.combin.csv")
Eminence.combined.sct@meta.data$CCA_Assignment_2 <- plyr::mapvalues(Eminence.combined.sct@meta.data$seurat_clusters, 
                                                                    from = E16_GE_assign_240222$Cluster, to = E16_GE_assign_240222$Assignment)
