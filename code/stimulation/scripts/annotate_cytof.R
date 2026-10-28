# Attach SAVED SOURCE per-cell labels; deliberately do not reconstruct annotation
# from newly obtained cluster IDs, or replace labels with Table S5 plurality.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
source(file.path(script_dir, "portable.R"))
opt <- parse_cli("Attach saved reference-cell annotations to a row-aligned CyTOF RDS.",
                 c("input", "annotations", "output", "output-table", "species"))
input <- option_path(opt, "input")
annotations <- option_path(opt, "annotations")
if (is.null(opt$output)) stop("--output RDS is required", call. = FALSE)
sp <- species_options(opt, allow_both = FALSE)
library(Seurat)
library(data.table)
so <- readRDS(input)
ann <- fread(annotations, data.table = FALSE, colClasses = c(donor = "character"))
stopifnot(nrow(ann) == ncol(so),
          all(c("donor", "stimulation", "CyTOF_annotation", "broad_lineage") %in% names(ann)))
source_donor <- neutral_donors(ann$donor, sp, "cytof")
object_donor <- neutral_donors(so$donor, sp, "cytof")
if (!identical(source_donor, object_donor) ||
    !identical(as.character(ann$stimulation), as.character(so$stimulation)))
  stop("Source annotation donor/stimulation row order does not match the RDS")
if ("cell" %in% names(ann) && !identical(as.character(ann$cell), colnames(so)))
  stop("Source annotation cell IDs do not match RDS column order")
stopifnot(!anyNA(ann$CyTOF_annotation), !anyNA(ann$broad_lineage))
so$donor <- object_donor
so$CyTOF_annotation <- ann$CyTOF_annotation
so$broad_lineage <- ann$broad_lineage
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
saveRDS(so, opt$output)
if (!is.null(opt[["output-table"]])) {
  dir.create(dirname(opt[["output-table"]]), showWarnings = FALSE, recursive = TRUE)
  write.csv(data.frame(cell = colnames(so), donor = so$donor,
                       stimulation = so$stimulation, CyTOF_annotation = so$CyTOF_annotation,
                       broad_lineage = so$broad_lineage), opt[["output-table"]], row.names = FALSE)
}
cat("Saved source fine and broad labels retained for", ncol(so), "reference cells.\n")
