library(Seurat)

stopifnot(packageVersion("SeuratObject") >= "5.0.0")

genes <- c("med-2", "med-1", "end-3", "end-1", "elt-7", "elt-2", "tbx-35", "ceh-51")
merge <- c(MSa = "MSx", MSp = "MSx", MSxa = "MSxx", MSxp = "MSxx")
stages <- c("earlyEMS", "EMS", "earlyEandMS", "E", "earlyEx", "Ex", "earlyExx", "Exx", "Exxx",
            "MS", "earlyMSx", "MSx", "MSa", "MSp", "earlyMSxx", "MSxa", "MSxp")

# Ahringer J, van der Burght SN, Carelli F, Appert A, Dong Y, Buttress T, Hill M (2024).
# C. elegans early embryo single-nucleus multiome dataset [Data set]. 
# Zenodo. doi:10.5281/zenodo.14261732
url <- "https://zenodo.org/api/records/14261732/files/C_elegans_wt_early_embryo_GEX_ATAC_multiome.rds/content"
path <- "data/raw/C_elegans_wt_early_embryo_GEX_ATAC_multiome.rds"
if (!file.exists(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    options(timeout = 3600)
    download.file(url, path, mode = "wb")
}
obj <- readRDS(path)
obj <- subset(obj, subset=cell_type %in% stages)

counts <- t(as.matrix(GetAssayData(obj, assay = "RNA", layer = "counts")[genes, ]))
stage <- ifelse(obj$cell_type %in% names(merge), merge[obj$cell_type], obj$cell_type)
write.table(
    data.frame(cell_id = colnames(obj), cell_type = stage, counts, check.names = FALSE), 
    "data/ems_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
