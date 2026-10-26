# Cluster + profile the rebuilt balanced CyTOF reference (per species).
# De novo annotation by marker profiles; label names harmonized to manuscript set.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
source(file.path(script_dir, "portable.R"))
opt <- parse_cli("Optional reference reclustering; not a reconstruction of saved source labels.",
                 "species")
input_dir <- option_path(opt, "input-dir", "data")
out_dir <- option_path(opt, "output-dir", "data", input = FALSE)
species <- species_options(opt)
library(Seurat)
library(data.table)

cluster_species <- function(species) {
  df <- fread(file.path(input_dir, paste0("cytof_", species, "_balanced.csv")), data.table = FALSE)
  df$donor <- neutral_donors(df$donor, species, "cytof")
  meta <- df[, c("donor", "stimulation")]
  expr <- df[, setdiff(colnames(df), c("donor", "stimulation"))]
  so <- CreateSeuratObject(counts = t(as.matrix(expr)), assay = "CyTOF",
                           meta.data = meta, project = species)
  so <- NormalizeData(so, normalization.method = "CLR", margin = 2, verbose = FALSE)
  so <- ScaleData(so, features = rownames(so), verbose = FALSE)
  so <- RunPCA(so, features = rownames(so), npcs = 20, verbose = FALSE)
  so <- FindNeighbors(so, dims = 1:20, verbose = FALSE)
  so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
  so <- RunUMAP(so, dims = 1:20, verbose = FALSE)
  saveRDS(so, file.path(out_dir, paste0("cytof_", species, "_balanced_clustered.rds")))

  # per-cluster mean on arcsinh scale (interpretable)
  cl <- so$seurat_clusters
  prof <- aggregate(as.matrix(expr), by = list(cluster = cl), FUN = mean)
  cat("\n=== ", species, " clusters:", length(unique(cl)), " ===\n")
  print(table(cl))
  list(obj = so, prof = prof, expr = expr, clusters = cl)
}


# compact top-marker report per cluster
report <- function(res, species) {
  prof <- res$prof
  rownames(prof) <- prof$cluster
  prof <- prof[, -1]
  cat("\n########## ", species, " top markers per cluster (arcsinh mean, z across clusters) ##########\n")
  z <- scale(as.matrix(prof))
  rownames(z) <- rownames(prof)
  for (cl in rownames(z)) {
    top <- sort(z[cl, ], decreasing = TRUE)[1:6]
    cat(sprintf("cluster %2s (n=%4d): %s\n", cl, sum(res$clusters == cl),
                paste(sprintf("%s=%.1f", names(top), top), collapse = ", ")))
  }
}
for (sp in species) report(cluster_species(sp), sp)
