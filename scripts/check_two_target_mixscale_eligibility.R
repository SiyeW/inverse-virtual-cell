# Diagnose the published DE-feature gate without weakening it or fitting scores.
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
future::plan(future::sequential)
options(Seurat.object.assay.version="v3",mc.cores=2)

find_features <- function(object,target_cells,control_cells) {
  stopifnot(length(target_cells)>=5L,length(control_cells)>=20L,
    length(intersect(target_cells,control_cells))==0L)
  Seurat:::TopDEGenesMixscape(object,ident.1=target_cells,ident.2=control_cells,
    de.assay="RNA",logfc.threshold=0,labels="analysis_group",verbose=FALSE,pval.cutoff=0.05)
}

self_test <- function() {
  set.seed(91)
  labels <- c(rep("TARGET",24),rep("NTC",40),rep("INTERGENIC",24))
  raw <- matrix(rpois(300*length(labels),10),nrow=300,
    dimnames=list(paste0("gene",1:300),paste0("cell",seq_along(labels))))
  raw[1:8,labels=="TARGET"] <- raw[1:8,labels=="TARGET"]+35
  meta <- data.frame(analysis_group=labels,row.names=colnames(raw))
  object <- NormalizeData(CreateSeuratObject(raw,meta.data=meta),verbose=FALSE)
  Idents(object) <- "analysis_group"
  for(control in c("NTC","INTERGENIC")) {
    features <- find_features(object,colnames(raw)[labels=="TARGET"],colnames(raw)[labels==control])
    stopifnot(length(intersect(features,paste0("gene",1:8)))>=5L)
  }
  cat("UNCHANGED_PUBLIC_DE_GATE_AND_SEPARATE_CONTROLS_VERIFIED\n")
}

run <- function(input,output,primary) {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  stage <- function(value) {
    writeLines(value,file.path(output,"phase.txt"))
    cat(value," ",format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),"\n",sep="");flush.console()
  }
  stage("read_accepted_derived_counts_once")
  prepared <- read_two_target_counts(input,seurat_keys=TRUE)
  genes <- prepared$genes
  object <- CreateSeuratObject(prepared$counts,meta.data=prepared$cells,min.cells=3,min.features=0)
  # The global expression mask is frozen before subsetting, as in the main run.
  object <- subset(object,cells=prepared$cells$cell_barcode[prepared$cells$pilot_record])
  rm(prepared);gc()
  object <- NormalizeData(object,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
  meta <- object@meta.data
  Idents(object) <- "analysis_group"
  counts <- table(meta$guide_id)
  eligible <- names(counts[counts>=100L])
  prior <- read.csv(file.path(primary,"mixscale_summary.csv"))
  results <- list();features <- list()
  settings <- data.frame(scenario=c("NTC_all_guides","NTC_guides_ge100",
      "intergenic_all_guides","intergenic_guides_ge100","no_guide_all_guides"),
    control=c("non_targeting_control","non_targeting_control","intergenic_control",
      "intergenic_control","no_guide_detected"),
    target_ge100=c(FALSE,TRUE,FALSE,TRUE,FALSE),
    control_ge100=c(FALSE,FALSE,FALSE,TRUE,FALSE))
  for(target in c("TYMSOS","UBA1")) for(index in seq_len(nrow(settings))) {
    setting <- settings[index,]
    stage(paste(target,setting$scenario,sep="_"))
    treated <- meta$target_gene_name %in% target
    control <- meta$role==setting$control
    if(setting$target_ge100) treated <- treated & meta$guide_id %in% eligible
    if(setting$control_ge100) control <- control & meta$guide_id %in% eligible
    selected <- find_features(object,rownames(meta)[treated],rownames(meta)[control])
    if(setting$scenario=="NTC_all_guides") {
      expected <- prior$initial_DE_features[prior$target_gene_name==target]
      stopifnot(length(selected)==expected)
    }
    key <- paste(target,setting$scenario,sep="_")
    results[[key]] <- data.frame(target_gene_name=target,scenario=setting$scenario,
      target_records=sum(treated),control_records=sum(control),
      target_guides=length(unique(meta$guide_id[treated])),
      DE_features=length(selected),required_features=5L,
      passes_feature_gate=length(selected)>=5L)
    if(length(selected)>0L) features[[key]] <- data.frame(target_gene_name=target,
      scenario=setting$scenario,seurat_key=selected,
      gene_id=genes$gene_id[match(selected,genes$seurat_key)])
    write.csv(do.call(rbind,results),file.path(output,"feature_eligibility.csv"),row.names=FALSE)
    write.csv(do.call(rbind,features),file.path(output,"eligible_RNAs.csv"),row.names=FALSE)
  }
  writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
  stage("feature_eligibility_complete")
}

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=2L,args[[1]] %in% c("--self-test","--run"))
source(args[[2]])
if(args[[1]]=="--self-test") self_test() else {
  stopifnot(length(args)==5L)
  run(args[[3]],args[[4]],args[[5]])
}
