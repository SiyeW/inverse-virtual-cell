# Change only the minimum DE-feature count, retaining frozen feature selection.
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(matrixStats))
future::plan(future::sequential)
options(Seurat.object.assay.version="v3",mc.cores=2)

load_projection_helper <- function(path) {
  definitions <- parse(file=path)
  selected <- Filter(function(value) is.call(value) &&
    identical(value[[1L]],as.name("<-")) &&
    identical(value[[2L]],as.name("extract_projection")),definitions)
  stopifnot(length(selected)==1L)
  eval(selected[[1L]],envir=parent.frame())
}

selected_normalized_expression <- function(counts,features) {
  # Reproduce CreateSeuratObject(min.cells=3), then LogNormalize(10000),
  # without materializing the normalized full-gene matrix again.
  retained <- Matrix::rowSums(counts>0)>=3L
  stopifnot(all(features %in% rownames(counts)[retained]))
  totals <- Matrix::colSums(counts)-Matrix::colSums(counts[!retained,,drop=FALSE])
  stopifnot(all(totals>0))
  normalized <- log1p(sweep(as.matrix(counts[features,,drop=FALSE]),2L,totals,"/")*10000)
  list(expression=normalized,retained_features=sum(retained),totals=totals)
}

fit_at_minimum <- function(signatures,meta,target,features,minimum) {
  base <- subset(signatures,cells=rownames(meta)[meta$analysis_group %in% c(target,"NTC")])
  warnings <- character()
  fit <- withCallingHandlers(
    RunMixscale(base,assay="PRTB",slot="scale.data",labels="analysis_group",nt.class.name="NTC",
      min.de.genes=minimum,logfc.threshold=0,de.assay="RNA",DE.gene=setNames(list(features),target),
      max.de.genes=100,new.class.name="mixscale_score",fine.mode=FALSE,verbose=FALSE,split.by=NULL),
    warning=function(value) {warnings <<- c(warnings,conditionMessage(value));invokeRestart("muffleWarning")})
  result <- extract_projection(fit,signatures,meta,target,features)
  saved <- fit@tools[["RunMixscale"]][[target]]
  loo_nonfinite <- if(is.null(saved) || !length(saved)) NA_integer_ else
    sum(!is.finite(as.matrix(saved[[1L]][,features,drop=FALSE])))
  list(scores=result,warnings=unique(warnings),loo_nonfinite=loo_nonfinite)
}

self_test <- function() {
  set.seed(813)
  ids <- sprintf("RNA%06d",seq_len(50))
  group <- c(rep("NTC",40),rep("GENE",20),rep("NO_GUIDE",10),rep("INTERGENIC",10))
  cells <- paste0("cell",seq_along(group))
  raw <- matrix(rpois(50*length(group),8),nrow=50,dimnames=list(ids,cells))
  raw[1,group=="GENE"] <- raw[1,group=="GENE"]+30
  raw[50,] <- 0;raw[50,1:2] <- 4
  raw <- as(raw,"dgCMatrix")
  normalized <- selected_normalized_expression(raw,ids[1:2])
  original <- NormalizeData(CreateSeuratObject(raw,min.cells=3),verbose=FALSE)
  stopifnot(isTRUE(all.equal(unname(normalized$expression),
    unname(as.matrix(GetAssayData(original,layer="data")[ids[1:2],])),tolerance=1e-12)))
  meta <- data.frame(cell_barcode=cells,guide_id=ifelse(group=="NO_GUIDE","",group),
    role=ifelse(group=="NTC","non_targeting_control",ifelse(group=="NO_GUIDE","no_guide_detected",
      ifelse(group=="INTERGENIC","intergenic_control","targeting"))),analysis_group=group,row.names=cells)
  object <- CreateSeuratObject(raw[ids[1:2],],meta.data=meta)
  object <- SetAssayData(object,assay="RNA",layer="data",new.data=normalized$expression)
  embedding <- matrix(rnorm(length(group)*10),ncol=10,
    dimnames=list(cells,paste0("PC_",1:10)))
  object[["pca"]] <- CreateDimReducObject(embeddings=embedding,key="PC_",assay="RNA")
  signatures <- CalcPerturbSig(object,assay="RNA",slot="data",features=ids[1:2],
    gd.class="analysis_group",nt.cell.class="NTC",reduction="pca",ndims=10,num.neighbors=20,verbose=FALSE)
  for(minimum in c(5L,3L,2L,1L)) {
    got <- fit_at_minimum(signatures,meta,"GENE",ids[1],minimum)
    if(minimum==1L) {
      stopifnot(nrow(got$scores)==80L,all(is.finite(got$scores$score)),
        any(grepl("larger than 1",got$warnings)),got$loo_nonfinite==60L)
      nt <- got$scores$score[got$scores$role=="non_targeting_control"]
      stopifnot(abs(mean(nt))<1e-10,abs(sd(nt)-1)<1e-10)
    } else stopifnot(all(is.na(got$scores$score)))
  }
  cat("ONE_FEATURE_PUBLIC_SCORE_AND_FROZEN_NORMALIZATION_VERIFIED\n")
}

run <- function(input,baseline,output) {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  stage <- function(value) {
    writeLines(value,file.path(output,"phase.txt"))
    cat(value," ",format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),"\n",sep="");flush.console()
  }
  checkpoint <- file.path(output,"single_gene_signatures.rds")
  if(file.exists(checkpoint)) {
    stage("resume_saved_signatures")
    saved <- readRDS(checkpoint)
    signatures <- saved$signatures;meta <- saved$meta;markers <- saved$markers
  } else {
    stage("read_accepted_derived_counts")
    prepared <- read_two_target_counts(input,seurat_keys=TRUE)
    stopifnot(nrow(prepared$cells)==81201L,sum(prepared$cells$pilot_record)==29850L)
    genes <- prepared$genes
    markers <- read.csv(file.path(baseline,"mixscale/scoring_RNAs.csv"),stringsAsFactors=FALSE)
    stopifnot(nrow(markers)==2L,setequal(markers$target_gene_name,c("TYMSOS","UBA1")),
      all(genes$gene_id[match(markers$seurat_key,genes$seurat_key)]==markers$gene_id))
    normalized <- selected_normalized_expression(prepared$counts,markers$seurat_key)
    stage("reuse_original_normalization_and_shared_PCA")
    selected <- prepared$cells$pilot_record
    meta <- prepared$cells[selected,,drop=FALSE];rownames(meta) <- meta$cell_barcode
    pca <- readRDS(file.path(baseline,"mixscale/shared_PCA.rds"))
    stopifnot(identical(rownames(pca$embeddings),prepared$cells$cell_barcode))
    object <- CreateSeuratObject(prepared$counts[markers$seurat_key,selected,drop=FALSE],meta.data=meta)
    object <- SetAssayData(object,assay="RNA",layer="data",new.data=normalized$expression[,selected,drop=FALSE])
    object[["pca"]] <- CreateDimReducObject(embeddings=pca$embeddings[meta$cell_barcode,,drop=FALSE],key="PC_",assay="RNA")
    write.csv(data.frame(input_records=nrow(prepared$cells),pilot_records=nrow(meta),
      retained_normalization_features=normalized$retained_features,selected_scoring_features=2L),
      file.path(output,"normalization_receipt.csv"),row.names=FALSE)
    rm(prepared,normalized,pca);gc()
    stage("same_twenty_NTC_neighbors_selected_gene_signatures")
    signatures <- CalcPerturbSig(object,assay="RNA",slot="data",features=markers$seurat_key,
      gd.class="analysis_group",nt.cell.class="NTC",reduction="pca",ndims=10,num.neighbors=20,verbose=FALSE)
    saveRDS(list(signatures=signatures,meta=meta,markers=markers),checkpoint)
    rm(object);gc()
  }
  summary <- list();warnings <- list()
  for(name in c("TYMSOS","UBA1")) {
    target <- unique(meta$analysis_group[meta$target_gene_name %in% name])
    feature <- markers$seurat_key[markers$target_gene_name==name]
    stopifnot(length(target)==1L,length(feature)==1L)
    for(minimum in c(5L,3L,2L,1L)) {
      stage(paste0(name,"_minimum_",minimum))
      got <- fit_at_minimum(signatures,meta,target,feature,minimum)
      scores <- got$scores
      active <- scores$analysis_group==target
      expected <- if(name=="TYMSOS") 681L else 538L
      stopifnot(sum(active)==expected,sum(scores$role=="non_targeting_control")==19342L,
        sum(scores$role=="no_guide_detected")==9094L,sum(scores$role=="intergenic_control")==195L)
      if(minimum==1L) stopifnot(all(is.finite(scores$score)))
      else stopifnot(all(is.na(scores$score)))
      key <- paste(name,minimum,sep="_")
      summary[[key]] <- data.frame(target_gene_name=name,minimum_DE_genes=minimum,available_DE_genes=1L,
        marker_gene_id=markers$gene_id[markers$target_gene_name==name],target_records=expected,
        scored_target_records=sum(is.finite(scores$score[active])),
        NTC_mean=mean(scores$score[scores$role=="non_targeting_control"]),
        NTC_sd=sd(scores$score[scores$role=="non_targeting_control"]),
        target_median=median(scores$score[active]),
        leave_one_gene_out_nonfinite=got$loo_nonfinite,
        assessment=if(minimum==1L) "exploratory_single_gene_projection" else "insufficient_DE_features")
      warnings[[key]] <- data.frame(condition=rep(key,length(got$warnings)),warning=got$warnings)
      if(minimum==1L) {
        scores$target_gene_name <- name;scores$minimum_DE_genes <- minimum
        scores$marker_gene_id <- markers$gene_id[markers$target_gene_name==name]
        expression <- GetAssayData(signatures,assay="RNA",layer="data")[feature,scores$cell_barcode]
        scores$marker_log_normalized_expression <- as.numeric(expression)
        write_gzip_csv(scores,file.path(output,paste0(name,"_minimum1_scores.csv.gz")))
      }
      write.csv(do.call(rbind,summary),file.path(output,"threshold_summary.csv"),row.names=FALSE)
      write.csv(do.call(rbind,warnings),file.path(output,"warnings.csv"),row.names=FALSE)
    }
  }
  writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
  stage("lower_minimum_complete")
}

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=4L,args[[1]] %in% c("--self-test","--run"))
source(args[[2]])
source(args[[3]])
load_projection_helper(args[[4]])
if(args[[1]]=="--self-test") self_test() else {
  stopifnot(length(args)==7L)
  run(args[[5]],args[[6]],args[[7]])
}
