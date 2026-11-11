# Run the published Mixscale scoring function with complete RNA features.
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(matrixStats))
future::plan(future::sequential)
options(future.globals.maxSize=8*1024^3,Seurat.object.assay.version="v3",mc.cores=2)

extract_projection <- function(fit, signatures, meta, target, features, control="NTC") {
  used <- meta$analysis_group==target | meta$role %in%
    c("non_targeting_control","intergenic_control","no_guide_detected")
  result <- meta[used,c("cell_barcode","guide_id","role","analysis_group"),drop=FALSE]
  result$target_gene_id <- target
  result$score <- NA_real_
  result$raw_projection <- NA_real_
  result$assessment <- "not_estimable"
  saved <- fit@tools[["RunMixscale"]][[target]]
  if(is.null(saved) || length(saved)==0L) return(result)
  original <- saved[[1L]]
  stopifnot(setequal(setdiff(colnames(original),c("pvec","gene")),features))
  training <- rownames(original)
  nt <- original$gene==control
  targeted <- original$gene==target
  stopifnot(sum(nt)>20L,sum(targeted)>0L,all(is.finite(original$pvec)))
  raw <- as.matrix(GetAssayData(signatures[["PRTB"]],layer="data")[features,training,drop=FALSE])
  scaled <- ScaleData(raw,verbose=FALSE)[rownames(raw),colnames(raw),drop=FALSE]
  direction <- rowMeans(scaled[,targeted,drop=FALSE])-rowMeans(scaled[,nt,drop=FALSE])
  squared_length <- sum(direction^2)
  stopifnot(is.finite(squared_length),squared_length>0)
  reproduced <- as.numeric(crossprod(direction,scaled))/squared_length
  stopifnot(isTRUE(all.equal(reproduced,as.numeric(original$pvec),tolerance=1e-7)))
  centers <- rowMeans(raw)
  deviations <- matrixStats::rowSds(raw)
  stopifnot(all(is.finite(deviations)),all(deviations>0))
  manual <- sweep(sweep(raw,1L,centers,"-"),1L,deviations,"/")
  manual[manual>10] <- 10
  stopifnot(isTRUE(all.equal(unname(manual),unname(scaled),tolerance=1e-7)))
  baseline_mean <- mean(original$pvec[nt])
  baseline_sd <- sd(original$pvec[nt])
  stopifnot(is.finite(baseline_sd),baseline_sd>0)
  expected <- (original$pvec[targeted]-baseline_mean)/baseline_sd
  reported <- fit@meta.data[training[targeted],"mixscale_score"]
  stopifnot(isTRUE(all.equal(as.numeric(reported),as.numeric(expected),tolerance=1e-7)))
  location <- match(training,result$cell_barcode)
  stopifnot(!anyNA(location))
  result$raw_projection[location] <- original$pvec
  result$score[location] <- (original$pvec-baseline_mean)/baseline_sd
  result$assessment[location] <- "computed"
  query <- which(!result$cell_barcode %in% training)
  if(length(query)>0L) {
    query_raw <- as.matrix(GetAssayData(signatures[["PRTB"]],layer="data")[features,result$cell_barcode[query],drop=FALSE])
    query_scaled <- sweep(sweep(query_raw,1L,centers,"-"),1L,deviations,"/")
    query_scaled[query_scaled>10] <- 10
    projection <- as.numeric(crossprod(direction,query_scaled))/squared_length
    result$raw_projection[query] <- projection
    result$score[query] <- (projection-baseline_mean)/baseline_sd
    result$assessment[query] <- "projected_additional_control"
  }
  stopifnot(all(is.finite(result$score)))
  result
}

self_test <- function() {
  previous <- options(warn=1)
  on.exit(options(previous))
  set.seed(71)
  ids <- sprintf("RNA%06d",seq_len(300))
  labels <- c(rep("NTC",40),rep("RESPONSE",20),rep("NULL",12),rep("NO_GUIDE",8))
  names <- paste0("cell",seq_along(labels))
  meta <- data.frame(cell_barcode=names,guide_id=ifelse(labels=="NO_GUIDE","",labels),
    role=ifelse(labels=="NTC","non_targeting_control",ifelse(labels=="NO_GUIDE","no_guide_detected","targeting")),
    analysis_group=labels,row.names=names)
  counts <- matrix(rpois(300*length(labels),rep(seq(2,20,length.out=300),length(labels))),
                   nrow=300,dimnames=list(ids,names))
  counts[seq_len(6),labels=="RESPONSE"] <- counts[seq_len(6),labels=="RESPONSE"]+25
  object <- CreateSeuratObject(counts,meta.data=meta)
  object <- NormalizeData(object,verbose=FALSE)
  classic <- ScaleData(object,features=rownames(object),verbose=FALSE)
  classic <- FindVariableFeatures(classic,selection.method="vst",nfeatures=100,verbose=FALSE)
  optimized <- FindVariableFeatures(object,selection.method="vst",nfeatures=100,verbose=FALSE)
  optimized <- ScaleData(optimized,features=VariableFeatures(optimized),verbose=FALSE)
  hvg <- VariableFeatures(optimized)
  stopifnot(identical(hvg,VariableFeatures(classic)),
    isTRUE(all.equal(GetAssayData(classic,layer="scale.data")[hvg,],
                     GetAssayData(optimized,layer="scale.data")[hvg,],tolerance=1e-7)))
  optimized <- RunPCA(optimized,npcs=10,verbose=FALSE)
  Idents(optimized) <- "analysis_group"
  de <- Seurat:::TopDEGenesMixscape(optimized,ident.1=names[labels=="RESPONSE"],ident.2=names[labels=="NTC"],
    de.assay="RNA",logfc.threshold=0,labels="analysis_group",verbose=FALSE,pval.cutoff=0.05)
  stopifnot(length(de)>=5L)
  nearest <- CalcPerturbSig(optimized,assay="RNA",features=union(hvg,head(de,100)),slot="data",
    gd.class="analysis_group",nt.cell.class="NTC",reduction="pca",ndims=10,num.neighbors=20,verbose=FALSE)
  stopifnot(all(head(de,100) %in% rownames(nearest[["PRTB"]])))
  values <- matrix(rnorm(300*length(labels)),nrow=300,dimnames=list(ids,names))
  values[seq_len(6),labels=="RESPONSE"] <- values[seq_len(6),labels=="RESPONSE"]+2
  object[["PRTB"]] <- CreateAssayObject(data=values)
  markers <- ids[seq_len(6)]
  base <- subset(object,cells=names[labels %in% c("NTC","RESPONSE")])
  fit <- suppressWarnings(RunMixscale(base,assay="PRTB",slot="scale.data",labels="analysis_group",nt.class.name="NTC",
    min.de.genes=5,logfc.threshold=0,DE.gene=list(RESPONSE=markers),verbose=FALSE))
  got <- extract_projection(fit,object,meta,"RESPONSE",markers)
  control <- got$score[got$role=="non_targeting_control"]
  stopifnot(nrow(got)==68L,abs(mean(control))<1e-10,abs(sd(control)-1)<1e-10,
            all(got$assessment[got$role=="no_guide_detected"]=="projected_additional_control"))
  base <- subset(object,cells=names[labels %in% c("NTC","NULL")])
  fit <- suppressWarnings(RunMixscale(base,assay="PRTB",slot="scale.data",labels="analysis_group",nt.class.name="NTC",
    min.de.genes=5,logfc.threshold=0,DE.gene=setNames(list(character()),"NULL"),verbose=FALSE))
  unscored <- extract_projection(fit,object,meta,"NULL",character())
  stopifnot(all(is.na(unscored$score)),all(unscored$assessment=="not_estimable"))
  cat("PUBLIC_SCORE_PROJECTIONS_AND_MISSING_VALUES_VERIFIED\n")
}

run <- function(input,output) {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  stage <- function(value) {
    writeLines(value,file.path(output,"phase.txt"))
    cat(value," ",format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),"\n",sep="");flush.console()
  }
  stage("read_whole_gene_counts")
  prepared <- read_two_target_counts(input,seurat_keys=TRUE)
  meta <- prepared$cells
  genes <- prepared$genes
  stopifnot(nrow(meta)==81201L,sum(meta$pilot_record)==29850L)
  rownames(meta) <- meta$cell_barcode
  object <- CreateSeuratObject(prepared$counts,meta.data=meta,min.cells=3,min.features=0)
  rm(prepared);gc()
  stage("shared_normalization_variable_features_and_PCA")
  object <- NormalizeData(object,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
  object <- FindVariableFeatures(object,selection.method="vst",nfeatures=2000,verbose=FALSE)
  # Scaling is gene-wise. Scaling the selected genes here gives their same
  # values as scaling the whole assay before selecting them, with bounded RAM.
  object <- ScaleData(object,features=VariableFeatures(object),verbose=FALSE)
  object <- RunPCA(object,npcs=50,seed.use=42,verbose=FALSE)
  saveRDS(list(embeddings=Embeddings(object,"pca"),loadings=Loadings(object,"pca"),
               variable_features=VariableFeatures(object)),file.path(output,"shared_PCA.rds"))
  object <- subset(object,cells=meta$cell_barcode[meta$pilot_record])
  meta <- meta[colnames(object),,drop=FALSE]
  gc()
  ntc <- rownames(meta)[meta$role=="non_targeting_control"]
  targets <- c("TYMSOS","UBA1")
  target_ids <- vapply(targets,function(name) {
    value <- unique(meta$analysis_group[meta$target_gene_name==name])
    stopifnot(length(value)==1L);value
  },character(1))
  markers <- list()
  initial_counts <- integer()
  stage("whole_gene_DE_features_for_public_Mixscale")
  Idents(object) <- "analysis_group"
  for(target in target_ids) {
    cells <- rownames(meta)[meta$analysis_group==target]
    value <- Seurat:::TopDEGenesMixscape(object,ident.1=cells,ident.2=ntc,de.assay="RNA",
      logfc.threshold=0,labels="analysis_group",verbose=FALSE,pval.cutoff=0.05)
    initial_counts[target] <- length(value)
    markers[[target]] <- head(value,100)
  }
  signature_features <- union(VariableFeatures(object),unique(unlist(markers,use.names=FALSE)))
  stage("original_20_neighbor_perturbation_signatures")
  signatures <- CalcPerturbSig(object,assay="RNA",slot="data",features=signature_features,
    gd.class="analysis_group",nt.cell.class="NTC",reduction="pca",ndims=10,
    num.neighbors=20,split.by=NULL,verbose=FALSE)
  rm(object);gc()
  all_scores <- list()
  summaries <- list()
  marker_tables <- list()
  for(index in seq_along(target_ids)) {
    target <- target_ids[[index]]
    name <- targets[[index]]
    stage(paste0("published_Mixscale_",name))
    base <- subset(signatures,cells=rownames(meta)[meta$analysis_group %in% c(target,"NTC")])
    provided <- setNames(list(markers[[target]]),target)
    fit <- RunMixscale(base,assay="PRTB",slot="scale.data",labels="analysis_group",nt.class.name="NTC",
      min.de.genes=5,logfc.threshold=0,de.assay="RNA",DE.gene=provided,max.de.genes=100,
      new.class.name="mixscale_score",fine.mode=FALSE,verbose=FALSE,split.by=NULL)
    scores <- extract_projection(fit,signatures,meta,target,markers[[target]])
    scores$target_gene_name <- name
    all_scores[[name]] <- scores
    active <- scores$analysis_group==target
    summaries[[name]] <- data.frame(target_gene_name=name,target_gene_id=target,
      target_records=sum(active),NTC_records=sum(scores$role=="non_targeting_control"),
      initial_DE_features=initial_counts[[target]],scoring_features=length(markers[[target]]),
      scored_target_records=sum(is.finite(scores$score[active])),
      status=if(all(is.finite(scores$score[active]))) "computed" else "insufficient_DE_features")
    marker_tables[[name]] <- data.frame(target_gene_name=rep(name,length(markers[[target]])),seurat_key=markers[[target]],
      gene_id=genes$gene_id[match(markers[[target]],genes$seurat_key)])
    write_gzip_csv(scores,file.path(output,paste0(name,"_scores.csv.gz")))
    rm(base,fit);gc()
  }
  write.csv(do.call(rbind,summaries),file.path(output,"mixscale_summary.csv"),row.names=FALSE)
  write.csv(do.call(rbind,marker_tables),file.path(output,"scoring_RNAs.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
  stage("mixscale_complete")
}

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=3L,args[[1]] %in% c("--self-test","--run"))
source(args[[2]])
source(args[[3]])
stopifnot(exists("RunMixscale",mode="function"),exists("read_two_target_counts",mode="function"))
if(args[[1]]=="--self-test") self_test() else {
  stopifnot(length(args)==5L)
  run(args[[4]],args[[5]])
}
