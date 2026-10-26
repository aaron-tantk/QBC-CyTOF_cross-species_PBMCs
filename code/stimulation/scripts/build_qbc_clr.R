# Rebuild QBC MARIO inputs: per-donor CLR normalization (Seurat, margin=2),
# matching Seurat-QBC-clustering.Rmd (raw counts -> CLR -> 10k cells/donor).
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
source(file.path(script_dir, "portable.R"))
opt <- parse_cli("Per-file QBC CLR normalization. Manifest columns: species,donor,path; file order is retained.",
                 c("input", "species"))
manifest_path <- option_path(opt, "input")
out_dir <- option_path(opt, "output-dir", "data", input = FALSE)
species <- species_options(opt)
library(Seurat)
library(data.table)

build_species <- function(files, donors, out_path) {
  mats <- list()
  for (i in seq_along(files)) {
    df <- fread(files[i], data.table = FALSE)
    if (!grepl("^CD|^HLA|^Ig", colnames(df)[1])) df <- df[, -1]  # drop index/cell-ID col
    stopifnot(ncol(df) == 134)
    so <- CreateSeuratObject(counts = t(as.matrix(df)), assay = "Protein", project = donors[i])
    so <- NormalizeData(so, normalization.method = "CLR", margin = 2, assay = "Protein", verbose = FALSE)
    mat <- t(as.matrix(GetAssayData(so, layer = "data")))
    mats[[i]] <- data.frame(mat, check.names = FALSE, donor = donors[i])
    cat(donors[i], ":", nrow(mat), "cells x", ncol(mat), "features; CD19 range:",
        round(range(mat[, "CD19"]), 3), "\n")
  }
  combined <- do.call(rbind, mats)
  fwrite(combined, out_path)
  cat("wrote", out_path, ":", nrow(combined), "cells\n\n")
  combined
}

manifest <- fread(manifest_path, data.table = FALSE, colClasses = "character")
stopifnot(all(c("species", "donor", "path") %in% names(manifest)),
          all(manifest$species %in% c("human", "rhesus")))
for (sp in species) {
  m <- manifest[manifest$species == sp, , drop = FALSE]
  if (nrow(m) == 0 || anyDuplicated(m$donor)) stop("Need unique per-file neutral donors for ", sp)
  if (!all(grepl(paste0("^", sp, "_qbc_[0-9]{2,}$"), m$donor)))
    stop("Manifest donor labels must use neutral species_qbc_01 indices")
  files <- ifelse(grepl("^(/|[A-Za-z]:)", m$path), m$path,
                  file.path(dirname(manifest_path), m$path))
  if (!all(file.exists(files))) stop("A manifest input file is missing")
  build_species(files, m$donor, file.path(out_dir, paste0("qbc_", sp, "_30k.csv")))
}
