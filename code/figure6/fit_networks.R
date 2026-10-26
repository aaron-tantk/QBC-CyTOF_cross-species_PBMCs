#!/usr/bin/env Rscript
# Optional refitting, NOT required to reproduce summaries of archived networks.
# Original bootstrap RNG state is not recorded; saved networks are authoritative.
args <- commandArgs(TRUE)
if (length(args) != 2L)
  stop("Usage: Rscript fit_networks.R PROCESSED_PRIMARY OUTPUT_DIRECTORY")
suppressPackageStartupMessages({library(data.table); library(bnlearn)})
input <- args[1]; output <- args[2]
dir.create(output, recursive=TRUE, showWarnings=FALSE)
# Historical name retained: 14 phospho epitopes plus total IκBα.
phospho <- c("pSTAT1.Y701.", "pSTAT3.Y705.", "IkBa.amino.terminal", "p4EBP1.T37.T46.",
 "pZAP70.Syk.Y319.Y352.", "pS6.S235.S236.", "pCREB.S133.", "pSTAT4.Y693.",
 "pSTAT5.Y694.", "pTBK1.NAK.S172.", "pSTAT6.Y691.", "pPLCg2.Y579.",
 "pERK1.2.T202.Y204.", "pP38.MAPK.T180.Y182.", "pMAPKAPK2.T334.")
all_summary <- list()
for (species in c("human", "rhesus")) {
  md <- fread(file.path(input,paste0(species,"_metadata.csv.gz")),data.table=FALSE)
  x <- fread(file.path(input,paste0(species,"_joint_CLR.csv.gz")),data.table=FALSE)
  stopifnot(identical(md$cell_id,x$cell_id))
  features <- paste0(phospho,"-CyTOF")
  stopifnot(all(features %in% names(x)))
  dat <- x[,features,drop=FALSE]; names(dat) <- phospho
  for (population in sort(unique(md$source_reference_broad_lineage))) {
    for (condition in c("Basal","CD40L","GM-CSF")) {
      keep <- md$source_reference_broad_lineage==population & md$condition==condition
      sub <- dat[keep,,drop=FALSE]
      if (nrow(sub)<100L) next
      sub <- sub[,vapply(sub,function(z) sd(z)>0,logical(1)),drop=FALSE]
      if (ncol(sub)<3L) next
      strength <- boot.strength(sub,R=500,algorithm="hc",algorithm.args=list(score="bic-g"))
      retained <- strength[strength$strength>0.85,,drop=FALSE]
      all_summary[[paste(species,population,condition,sep="|")]] <- retained
      write.csv(retained,file.path(output,paste0("fig6_arcs_",species,"_",
        gsub("[/ ]","_",population),"_",condition,".csv")),row.names=FALSE)
    }
  }
}
saveRDS(all_summary,file.path(output,"fig6_all_networks.rds"))
writeLines(capture.output(sessionInfo()),file.path(output,"refit_sessionInfo.txt"))
