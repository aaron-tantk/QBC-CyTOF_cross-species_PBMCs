# Cross-method concordance (Bokai comment 1): Harmony (label-free) integration of
# QBC + CyTOF on the 18 shared anchors, 1-NN label/partner transfer, compared to MARIO.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
source(file.path(script_dir, "portable.R"))
opt <- parse_cli("Optional Harmony-1NN comparator (full reference). Run MARIO --prepare-only first.",
                 c("species", "matches-dir"))
species <- species_options(opt, allow_both = FALSE)
DATA <- option_path(opt, "input-dir", "mario_out")
OUT <- option_path(opt, "output-dir", "results", input = FALSE)
MATCHES <- option_path(opt, "matches-dir", "mario_out")
library(Seurat)
library(harmony)
library(data.table)
library(dplyr)
library(aricode)

anchors <- if (species == "human") {
  c("CD11b","CD11c","CD123","CD14","CD16","CD161","CD19","CD1c","CD3","CD33",
    "CD4","CD45","CD45RA","CD56","CD7","CD8","HLA-DR","IgM")
} else {
  c("CD11b","CD11c","CD123","CD14","CD16","CD161","CD20","CD1c","CD3","CD33",
    "CD4","CD45","CD45RA","CD56","CD7","CD8","HLA-DR","IgM")
}

qbc <- fread(file.path(DATA, paste0("mario_input_qbc_", species, ".csv")), data.table = FALSE)
cyto <- fread(file.path(DATA, paste0("mario_input_cytof_", species, ".csv")), data.table = FALSE)
qbc$donor <- neutral_donors(qbc$donor, species, "qbc")
cyto$donor <- neutral_donors(cyto$donor, species, "cytof")
cyto <- cyto %>% rename_with(~ "HLA-DR", .cols = matches("^HLA.DR$"))

# combined 18-anchor matrix; z-score each modality separately (different dynamic ranges)
scale_mod <- function(df, feats) scale(as.matrix(df[, feats]))
Q <- scale_mod(qbc, anchors)
C <- scale_mod(cyto, anchors)
combined <- rbind(Q, C)
modality <- factor(c(rep("QBC", nrow(Q)), rep("CyTOF", nrow(C))))

# Harmony integration via Seurat workflow (harmony2 API)
colnames(combined) <- anchors
cell_names <- c(paste0("QBC_", seq_len(nrow(Q))), paste0("CyTOF_", seq_len(nrow(C))))
rownames(combined) <- cell_names
meta <- data.frame(modality = modality, row.names = cell_names)
so <- CreateSeuratObject(counts = t(combined), assay = "Anchor", meta.data = meta)
so[["Anchor"]]$data <- so[["Anchor"]]$counts  # values already z-scored per modality
so <- ScaleData(so, features = rownames(so), verbose = FALSE)
so <- RunPCA(so, features = rownames(so), npcs = 15, approx = FALSE, verbose = FALSE)
so <- RunHarmony(so, group.by.vars = "modality", verbose = FALSE)
emb <- Embeddings(so, "harmony")

# 1-NN: for each QBC cell, nearest CyTOF cell in Harmony space
library(FNN)
qbc_emb <- emb[modality == "QBC", , drop = FALSE]
cyto_emb <- emb[modality == "CyTOF", , drop = FALSE]
nn <- get.knnx(cyto_emb, qbc_emb, k = 1)
cyto_match_idx <- as.vector(nn$nn.index)

# transferred labels + partner info
harmony_transfer <- data.frame(
  qbc_row = seq_len(nrow(Q)),
  cyto_idx = cyto_match_idx,
  harmony_broad = cyto$broad_lineage[cyto_match_idx],
  harmony_label = cyto$CyTOF_annotation[cyto_match_idx],
  harmony_donor = cyto$donor[cyto_match_idx],
  harmony_stim = cyto$stimulation[cyto_match_idx],
  nn_dist = as.vector(nn$nn.dist))

# MARIO partners for the same QBC cells
m <- fread(file.path(MATCHES, paste0(species, "_matches.csv")), data.table = FALSE)
mario_map <- data.frame(qbc_row = m$row_idx_0 + 1, mario_cyto_idx = m$row_idx_1 + 1)  # 1-based
mario_map$mario_broad <- cyto$broad_lineage[mario_map$mario_cyto_idx]

cmp <- harmony_transfer %>% left_join(mario_map, by = "qbc_row")
cmp <- cmp %>% filter(!is.na(mario_broad), !is.na(harmony_broad))  # matched by MARIO + valid Harmony label
cmp$harmony_broad <- as.factor(cmp$harmony_broad)
cmp$mario_broad <- as.factor(cmp$mario_broad)

# agreement metrics
agree <- mean(as.character(cmp$harmony_broad) == as.character(cmp$mario_broad))
ari <- ARI(cmp$harmony_broad, cmp$mario_broad)
nmi <- NMI(cmp$harmony_broad, cmp$mario_broad)
cat(sprintf("\n=== %s: Harmony-1NN vs MARIO partner broad-lineage agreement ===\n", species))
cat(sprintf("n QBC cells compared: %d | pct agreement: %.1f%% | ARI: %.3f | NMI: %.3f\n",
            nrow(cmp), 100 * agree, ari, nmi))
cat("\nConfusion (row=Harmony, col=MARIO), %:\n")
print(round(prop.table(table(cmp$harmony_broad, cmp$mario_broad)) * 100, 1))

fwrite(cmp, file.path(OUT, paste0(species, "_harmony_vs_mario_percell.csv")))
saveRDS(list(emb = emb, modality = modality, harmony_transfer = harmony_transfer),
        file.path(OUT, paste0(species, "_harmony_embedding.rds")))
