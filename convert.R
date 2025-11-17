library(Seurat)
library(SeuratDisk)
library(MuDataSeurat)


setwd("~/mtr2/tfm/code")
filepath <- "daty/032224_L4_all_cells_Seurat5.rds"
# seur <- readRDS("daty/032224_L4_all_cells_Seurat5.rds")
seur <- readRDS(filepath)
seurat_object_new <- seur
seurat_object_new[["RNA"]] <- as(object = seurat_object_new[["RNA"]], Class = "Assay") 

seur_h5s <- SeuratDisk::as.h5Seurat(seurat_object_new, "test.h5Seurat", overwrite = T)
# seur_h5s <- SeuratDisk::LoadH5Seurat("test.h5Seurat")
seur_h5ad <- SeuratDisk::Convert(seur_h5s,dest = "h5ad", overwrite = T)



# MuDataSeurat::WriteH5MU(seurat_object_new, "test.h5mu", )
# seur_mu <- MuDataSeurat::ReadH5MU("test.h5mu")


# Trabajo con seurat directamente

seurat_object_new


unique(seurat_object_new$orig.ident)
seurat_object_new$Cell.type <- as.factor(seurat_object_new$Cell.type)

sort(table(seurat_object_new$Cell.type), decreasing = T)


