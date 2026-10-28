#!/usr/bin/env Rscript
# Original-style fused trajectory; no redesign of other manuscript figures.
args <- commandArgs(TRUE)
stopifnot(length(args)==3)
input <- args[1];out <- args[2];graph <- args[3]
.libPaths(c(Sys.getenv("TRAJECTORY_RLIB"),.Library))
suppressPackageStartupMessages({
  library(monocle3);library(igraph);library(data.table);library(jsonlite)
})
stopifnot(as.character(packageVersion("monocle3"))=="1.3.1",
          as.character(packageVersion("igraph"))=="1.6.0")
RhpcBLASctl::blas_set_num_threads(1);RhpcBLASctl::omp_set_num_threads(1)
options(mc.cores=1,expressions=500000)
seed <- 42L;set.seed(seed)
species <- sub("_.*","",graph);condition <- sub("^[^_]+_","",graph)
stopifnot(species%in%c("human","rhesus"),condition%in%c("Basal","GM-CSF"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(sub in c("objects","tables","environment"))dir.create(file.path(out,sub),showWarnings=FALSE)
meta <- as.data.frame(fread(file.path(input,paste0(species,"_metadata.csv.gz"))))
x <- as.data.frame(fread(file.path(input,paste0(species,"_matrix.csv.gz"))))
stopifnot(identical(meta$cell_id,x$cell_id),ncol(x)==175,!anyDuplicated(x$cell_id))
sel <- meta$reference_group==condition
meta <- meta[sel,,drop=FALSE];x <- x[sel,,drop=FALSE]
stopifnot(all(meta$condition==condition),!any(meta$condition=="CD40L"))
mat <- t(as.matrix(x[,-1]));colnames(mat)<-meta$cell_id
stopifnot(all(is.finite(mat)),ncol(mat)>500,all(is.finite(meta$nine_marker_score)))
seeded_function <- function(fun) {
  txt<-paste(deparse(body(fun)),collapse="\n")
  stopifnot(grepl("set.seed(2016)",txt,fixed=TRUE))
  body(fun)<-parse(text=gsub("set.seed(2016)","set.seed(42)",txt,fixed=TRUE))[[1]]
  fun
}
started<-Sys.time()
cat("START",graph,ncol(mat),"cells",nrow(mat),"features\n");flush.console()
cds<-new_cell_data_set(mat,
       cell_metadata=as.data.frame(meta,row.names=meta$cell_id),
       gene_metadata=data.frame(gene_short_name=rownames(mat),row.names=rownames(mat)))
cds<-estimate_size_factors(cds)
preprocess_seeded<-seeded_function(preprocess_cds)
reduce_seeded<-seeded_function(reduce_dimension)
cds<-preprocess_seeded(cds,num_dim=50,norm_method="none")
cat("PCA complete\n");flush.console()
cds<-reduce_seeded(cds,reduction_method="UMAP",preprocess_method="PCA",cores=1L)
cat("UMAP complete\n");flush.console()
cds<-cluster_cells(cds,reduction_method="UMAP",random_seed=seed)
cds<-learn_graph(cds,use_partition=TRUE)
threshold<-unname(quantile(meta$nine_marker_score,.95))
root_cells<-meta$cell_id[meta$nine_marker_score>=threshold]
stopifnot(length(root_cells)>0,all(root_cells%in%colnames(cds)))
cds<-order_cells(cds,root_cells=root_cells)
pt<-pseudotime(cds)[colnames(cds)]
xy<-reducedDims(cds)$UMAP
celltable<-cbind(meta,UMAP1=xy[,1],UMAP2=xy[,2],
                 pseudotime=as.numeric(pt),finite_pseudotime=is.finite(pt),
                 root_cell=meta$cell_id%in%root_cells,partition=as.character(partitions(cds)))
stopifnot(any(is.finite(pt)))
g<-principal_graph(cds)[["UMAP"]]
nodes<-t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
edges<-as_edgelist(g,names=TRUE)
fwrite(celltable,file.path(out,"tables",paste0(graph,"_cells.csv.gz")))
fwrite(data.frame(node=rownames(nodes),UMAP1=nodes[,1],UMAP2=nodes[,2]),
       file.path(out,"tables",paste0(graph,"_nodes.csv")))
fwrite(data.frame(from=edges[,1],to=edges[,2]),
       file.path(out,"tables",paste0(graph,"_edges.csv")))
# Remove optional external caches, not graph/coordinates/order, for serialization.
for(reduction in names(cds@reduce_dim_aux)) {
  cds@reduce_dim_aux[[reduction]][["nn_index"]]<-NULL
  if(reduction=="UMAP")cds@reduce_dim_aux[[reduction]][["model"]][["umap_model"]]<-NULL
}
path<-file.path(out,"objects",paste0(graph,".rds"))
saveRDS(list(cds=cds,roots=root_cells,threshold=threshold,seed=seed),path)
reopened<-readRDS(path)
again<-order_cells(reopened$cds,root_cells=reopened$roots)
stopifnot(identical(as.numeric(pseudotime(again)[colnames(cds)]),as.numeric(pt)))
finite<-is.finite(pt)
rho<-if(sum(finite)>2) cor(meta$nine_marker_score[finite],pt[finite],method="spearman") else NA_real_
params<-list(graph=graph,n_cells=ncol(cds),n_features=nrow(cds),num_dim=50,
  root_rule="top_5_percent_of_common_species_nine_marker_score",root_threshold=threshold,
  n_roots=length(root_cells),n_ordered=sum(finite),finite_fraction=mean(finite),
  partitions=length(unique(partitions(cds))),n_nodes=vcount(g),n_edges=ecount(g),
  rooting_consistency_spearman=unname(rho),
  elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),
  monocle3=as.character(packageVersion("monocle3")),igraph=as.character(packageVersion("igraph")),
  seed=seed,serialization_and_ordering_readback="PASS",UMAP_n_neighbors=15,UMAP_min_dist=.1,
  UMAP_metric="cosine",use_partition=TRUE,close_loop=TRUE)
write_json(params,file.path(out,"tables",paste0(graph,"_summary.json")),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,"environment",paste0(graph,"_sessionInfo.txt")))
for(fun in c("preprocess_cds","reduce_dimension","cluster_cells","learn_graph"))
  writeLines(capture.output(formals(get(fun))),file.path(out,"environment",paste0(fun,"_defaults.txt")))
cat(toJSON(params,auto_unbox=TRUE,pretty=TRUE),"\n")
