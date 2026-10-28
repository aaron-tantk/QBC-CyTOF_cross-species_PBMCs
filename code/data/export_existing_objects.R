#!/usr/bin/env Rscript
# Export existing fused objects without normalization changes or relabeling.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
suppressPackageStartupMessages({library(Seurat); library(data.table); library(jsonlite)})
source_dir <- args[1]; audit_dir <- args[2]; output <- args[3]
dir.create(output, recursive=TRUE, showWarnings=FALSE)
receipt <- list()
for (sp in c("human", "rhesus")) {
  prefix <- if (sp == "human") "H" else "M"
  files <- file.path(source_dir, paste0("fused_", sp, c("_pbmc", ""), ".csv"))
  donors <- lapply(files, function(p) fread(p, select=c("QBC_donor", "CyTOF_donor")))
  qdon <- sort(unique(unlist(lapply(donors, function(x) as.character(x$QBC_donor)))))
  cdon <- sort(unique(unlist(lapply(donors, function(x) as.character(x$CyTOF_donor)))))
  for (design in c("primary", "whole_blood_sensitivity")) {
    tag <- if (design == "primary") "_pbmc" else ""
    out <- file.path(output, design); dir.create(out, recursive=TRUE, showWarnings=FALSE)
    f <- fread(file.path(source_dir, paste0("fused_", sp, tag, ".csv")), data.table=FALSE)
    obj <- readRDS(file.path(source_dir, paste0("fused_", sp, tag, "_clustered.rds")))
    md <- obj[[]]
    features <- grep("-(QBC|CyTOF)$", names(f), value=TRUE)
    stopifnot(length(features)==174L, sum(endsWith(features, "-QBC"))==134L,
              !anyDuplicated(f$qbc_row_id), nrow(f)==nrow(md))
    counts <- as.matrix(GetAssayData(obj, assay="Protein", layer="counts"))
    joint <- as.matrix(GetAssayData(obj, assay="Protein", layer="data"))
    stopifnot(identical(rownames(counts), features),
              identical(dim(counts), dim(joint)),
              max(abs(t(as.matrix(f[,features]))-counts)) < 1e-12)
    recomputed <- log1p(sweep(counts, 2L, exp(colSums(log1p(counts))/nrow(counts)), "/"))
    stopifnot(max(abs(joint-recomputed)) < 1e-12)
    for (key in c("qbc_row_id", "cytof_row_id", "stimulation", "CyTOF_annotation",
                  "CyTOF_broad_lineage", "QBC_donor", "CyTOF_donor"))
      stopifnot(identical(as.character(md[[key]]), as.character(f[[key]])))
    cells <- sprintf("%s-QBC-%06d", prefix, f$qbc_row_id)
    metadata <- data.frame(
      cell_id=cells, reference_id=sprintf("%s-CYTOF-%06d", prefix, f$cytof_row_id),
      species=sp, design=design,
      qbc_donor=sprintf("%s-Q%02d", prefix, match(as.character(f$QBC_donor), qdon)),
      cytof_donor=sprintf("%s-C%02d", prefix, match(as.character(f$CyTOF_donor), cdon)),
      condition=f$stimulation, source_reference_population=f$CyTOF_annotation,
      source_reference_broad_lineage=f$CyTOF_broad_lineage,
      qbc_row_id=f$qbc_row_id, cytof_row_id=f$cytof_row_id,
      within_pair_anchor_r=f$within_pair_anchor_r,
      fused_cluster_res0.5=as.character(md[["Protein_snn_res.0.5"]]),
      fused_cluster_res1.2=as.character(md[["Protein_snn_res.1.2"]]),
      check.names=FALSE)
    mapping <- fread(file.path(audit_dir, paste0("Table_S5_regenerated_", sp, tag, "_res0.5.csv")))
    metadata$audit_plurality_population <- mapping$transferred_CyTOF_label[
      match(metadata$fused_cluster_res0.5, as.character(mapping$cluster))]
    stopifnot(!anyNA(metadata$audit_plurality_population))
    fwrite(metadata, file.path(out, paste0(sp, "_metadata.csv.gz")))
    fwrite(data.frame(cell_id=cells, t(counts), check.names=FALSE),
           file.path(out, paste0(sp, "_pre_joint_CLR.csv.gz")))
    fwrite(data.frame(cell_id=cells, t(joint), check.names=FALSE),
           file.path(out, paste0(sp, "_joint_CLR.csv.gz")))
    key <- paste(sp, design, sep="_")
    receipt[[key]] <- list(n_rows=nrow(f), n_features=length(features),
      condition_counts=as.list(table(metadata$condition)),
      source_population_count=length(unique(metadata$source_reference_population)),
      plurality_population_count=length(unique(metadata$audit_plurality_population)),
      source_plurality_disagreements=sum(metadata$source_reference_population != metadata$audit_plurality_population),
      pre_joint_values_exact=TRUE, joint_CLR_formula_verified=TRUE,
      metadata_source_exact=TRUE, donor_identity_lookup_included=FALSE)
    rm(f,obj,md,counts,joint,recomputed);gc()
  }
}
write_json(receipt, file.path(output,"export_validation.json"), auto_unbox=TRUE, pretty=TRUE)
cat(toJSON(receipt, auto_unbox=TRUE, pretty=TRUE), "\n")
