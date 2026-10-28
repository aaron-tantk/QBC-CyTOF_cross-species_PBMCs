# Cluster corrected fused datasets on the 59 cross-species QBC markers,
# marker-verify clusters, and regenerate the Table S5 purity audit.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
source(file.path(script_dir, "portable.R"))
opt <- parse_cli("Optional fused clustering and Table S5 composition audit; source labels are NOT overwritten.",
                 c("species", "tag", "markers", "rds-output-dir"))
species <- species_options(opt, allow_both = FALSE)
tag <- if (is.null(opt$tag) || opt$tag == "") "" else paste0("_", opt$tag)
if (!tag %in% c("", "_pbmc")) stop("--tag must be pbmc or omitted")
FUSED <- option_path(opt, "input-dir", "fused")
OUT <- option_path(opt, "output-dir", "results", input = FALSE)
RDS_OUT <- if (is.null(opt[["rds-output-dir"]])) OUT else
  option_path(opt, "rds-output-dir", input = FALSE)
marker_path <- if (is.null(opt$markers))
  file.path(script_dir, "..", "markers", "Table_S4_cross_species_59_markers.csv") else opt$markers
if (!file.exists(marker_path)) stop("Marker panel not found; use --markers")
library(Seurat)
library(data.table)
library(dplyr)

# --- 59 cross-species markers (Table S4), mapped to QBC column names ---
s4 <- fread(marker_path)
s4_names <- s4$Antibody
name_map <- c("FcεRIα" = "FcRI", "HLA-A,B,C" = "HLA-ABC", "Integrin β7" = "integrinB7")
s4_qbc <- ifelse(s4_names %in% names(name_map), name_map[s4_names], s4_names)

# --- load fused pairs ---
fused <- fread(file.path(FUSED, paste0("fused_", species, tag, ".csv")), data.table = FALSE)
fused$QBC_donor <- neutral_donors(fused$QBC_donor, species, "qbc")
fused$CyTOF_donor <- neutral_donors(fused$CyTOF_donor, species, "cytof")
meta_cols <- c("QBC_donor", "CyTOF_donor", "stimulation", "CyTOF_annotation",
               "CyTOF_broad_lineage", "qbc_row_id", "cytof_row_id", "within_pair_anchor_r")
meta <- fused[, meta_cols]
expr <- fused[, setdiff(colnames(fused), meta_cols)]
cat(species, "fused:", nrow(fused), "pairs x", ncol(expr), "features\n")

so <- CreateSeuratObject(counts = t(as.matrix(expr)), assay = "Protein", meta.data = meta)
# CLR normalize (as in the original stimulations Rmd), then scale the 59 QBC markers
so <- NormalizeData(so, normalization.method = "CLR", margin = 2, verbose = FALSE)
xs_markers <- paste0(s4_qbc, "-QBC")
xs_markers <- intersect(xs_markers, rownames(so))
cat("cross-species markers used:", length(xs_markers), "\n")
so <- ScaleData(so, features = xs_markers, verbose = FALSE)
so <- RunPCA(so, features = xs_markers, npcs = 20, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:20, verbose = FALSE)
so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
so <- FindClusters(so, resolution = 1.2, verbose = FALSE)  # finer, for purity sensitivity
so <- RunUMAP(so, dims = 1:20, verbose = FALSE)
saveRDS(so, file.path(RDS_OUT, paste0("fused_", species, tag, "_clustered.rds")))

# --- marker verification per cluster (QBC-side decisive markers, arcsinh-free CLR scale) ---
decisive <- paste0(c("CD19","CD20","CD3","CD4","CD8","CD14","CD16","CD56","HLA-DR","CD11c",
                     "CD123","CD1c","CD41","CD61","CD66b","CD45RA","CCR7","CD161","CD33","CD11b","IgM"), "-QBC")
decisive <- intersect(decisive, rownames(so))
raw <- GetAssayData(so, layer = "counts")[decisive, ]  # fused CLR input values
cl <- so$Protein_snn_res.0.5
prof <- sapply(decisive, function(m) tapply(raw[m, ], cl, mean))
colnames(prof) <- gsub("-QBC", "", colnames(prof))
cat("\n=== ", species, " res0.5 cluster marker profiles (mean) ===\n")
print(round(prof, 2))
write.csv(round(prof, 3), file.path(OUT, paste0(species, tag, "_cluster_marker_profiles.csv")))

# --- Table S5 regeneration: per cluster, CyTOF partner composition + concordance ---
make_s5 <- function(so, cluster_col) {
  md <- so@meta.data
  md$cl <- md[[cluster_col]]
  rows <- lapply(split(md, md$cl), function(d) {
    tab <- sort(table(d$CyTOF_annotation), decreasing = TRUE)
    top <- names(tab)[1]
    pct <- round(100 * as.numeric(tab[1]) / nrow(d), 1)
    comp <- paste(sprintf("%s %d%%", names(tab)[1:min(3, length(tab))],
                          round(100 * as.numeric(tab[1:min(3, length(tab))]) / nrow(d))), collapse = "; ")
    data.frame(cluster = d$cl[1], n_cells = nrow(d),
               transferred_CyTOF_label = top, purity_pct = pct,
               top_CyTOF_composition = comp,
               mean_within_pair_anchor_r = round(mean(d$within_pair_anchor_r), 3),
               frac_pair_r_gt_0.5 = round(mean(d$within_pair_anchor_r > 0.5), 3))
  })
  do.call(rbind, rows)
}
s5_new <- make_s5(so, "Protein_snn_res.0.5")
write.csv(s5_new, file.path(OUT, paste0("Table_S5_regenerated_", species, tag, "_res0.5.csv")), row.names = FALSE)
s5_fine <- make_s5(so, "Protein_snn_res.1.2")
write.csv(s5_fine, file.path(OUT, paste0("Table_S5_regenerated_", species, tag, "_res1.2.csv")), row.names = FALSE)
cat("\n=== regenerated Table S5 (", species, ", res0.5) ===\n")
print(s5_new)
cat("\nres1.2 purity summary: median", median(s5_fine$purity_pct), "% across", nrow(s5_fine), "clusters\n")
