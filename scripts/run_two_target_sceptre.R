# Whole-gene, per-guide SCEPTRE tests for the declared two-target pilot.
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(sceptre))
stopifnot(packageVersion("sceptre")=="0.10.3")

configure_object <- function(y,g,meta,guides) {
  mapping <- data.frame(grna_id=rownames(g),
    grna_target=ifelse(rownames(g) %in% unique(meta$guide_id[meta$role=="non_targeting_control"]),
                      "non-targeting",rownames(g)))
  pairs <- expand.grid(grna_target=guides,response_id=rownames(y),stringsAsFactors=FALSE)
  object <- import_data(response_matrix=y,grna_matrix=g,grna_target_data_frame=mapping,moi="low")
  object <- set_analysis_parameters(object,discovery_pairs=pairs,grna_integration_strategy="singleton",
    side="both",control_group="nt_cells",resampling_mechanism="permutations",multiple_testing_alpha=0.1)
  object <- assign_grnas(object,method="maximum",min_grna_n_umis_threshold=1L,
                        umi_fraction_threshold=0.8,print_progress=FALSE)
  actual <- rep(NA_character_,nrow(meta))
  for(guide in names(object@initial_grna_assignment_list)) {
    actual[object@initial_grna_assignment_list[[guide]]] <- guide
  }
  stopifnot(identical(actual,meta$guide_id))
  object <- run_qc(object,n_nonzero_trt_thresh=0L,n_nonzero_cntrl_thresh=7L,
    response_n_umis_range=c(0,1),response_n_nonzero_range=c(0,1),p_mito_threshold=1)
  stopifnot(identical(object@cells_in_use,seq_len(nrow(meta))),
    setequal(object@grna_assignments$all_nt_idxs,which(meta$role=="non_targeting_control")))
  object
}

validate_calibration_controls <- function(calibration,nt_guides) {
  # Singleton results store the actual control guide in grna_id. grna_target
  # deliberately holds the shared biological label "non-targeting".
  stopifnot(nrow(calibration)>0L,
    all(c("grna_id","grna_target") %in% colnames(calibration)),
    all(calibration$grna_id %in% nt_guides),
    all(calibration$grna_target=="non-targeting"))
  invisible(TRUE)
}

validate_checkpoint_records <- function(object,meta) {
  # import_data intentionally drops response-matrix column names. Verify the
  # ordered RNA and guide measurements instead of expecting absent barcodes.
  covariates <- object@covariate_data_frame
  actual <- rep(NA_character_,nrow(meta))
  for(guide in names(object@initial_grna_assignment_list)) {
    actual[object@initial_grna_assignment_list[[guide]]] <- guide
  }
  stopifnot(ncol(sceptre:::get_response_matrix(object))==nrow(meta),
    identical(as.numeric(covariates$response_n_umis),as.numeric(meta$total_counts)),
    identical(as.numeric(covariates$grna_n_umis),as.numeric(meta$guide_umis)),
    identical(actual,meta$guide_id),
    identical(object@cells_in_use,seq_len(nrow(meta))),
    setequal(object@grna_assignments$all_nt_idxs,which(meta$role=="non_targeting_control")))
  invisible(TRUE)
}

self_test <- function() {
  set.seed(711)
  group <- c(rep("N1",20),rep("N2",20),rep("T1",16))
  meta <- data.frame(cell_barcode=paste0("cell",seq_along(group)),guide_id=group,
    role=ifelse(group=="T1","targeting","non_targeting_control"))
  rates <- c(8,seq(0.25,5,length.out=99))
  y <- Matrix(matrix(rpois(100*length(group),rep(rates,length(group))),nrow=100),sparse=TRUE)
  rownames(y) <- paste0("gene",seq_len(nrow(y)))
  colnames(y) <- meta$cell_barcode
  y[1,group=="T1"] <- 0
  meta$total_counts <- Matrix::colSums(y)
  meta$guide_umis <- 1
  g <- sparseMatrix(i=match(group,c("N1","N2","T1")),j=seq_along(group),x=1,
    dims=c(3,length(group)),dimnames=list(c("N1","N2","T1"),meta$cell_barcode))
  object <- configure_object(y,g,meta,"T1")
  object <- run_calibration_check(object,n_calibration_pairs=4L,calibration_group_size=1L,
    parallel=FALSE,n_processors=2L,print_progress=FALSE,output_amount=2L)
  calibration <- get_result(object,"run_calibration_check")
  validate_calibration_controls(calibration,c("N1","N2"))
  incorrect <- as.data.frame(calibration)
  incorrect$grna_id[1] <- "T1"
  stopifnot(inherits(try(validate_calibration_controls(incorrect,c("N1","N2")),silent=TRUE),"try-error"))
  checkpoint <- tempfile(fileext=".rds")
  saveRDS(object,checkpoint)
  object <- readRDS(checkpoint)
  unlink(checkpoint)
  stopifnot(identical(as.data.frame(get_result(object,"run_calibration_check")),as.data.frame(calibration)))
  validate_checkpoint_records(object,meta)
  reversed <- meta[rev(seq_len(nrow(meta))),,drop=FALSE]
  stopifnot(inherits(try(validate_checkpoint_records(object,reversed),silent=TRUE),"try-error"))
  object <- run_discovery_analysis(object,parallel=FALSE,n_processors=2L,print_progress=FALSE,output_amount=2L)
  result <- get_result(object,"run_discovery_analysis")
  stopifnot(nrow(result)>0L,all(result$grna_target=="T1"),
            all(result$p_value[is.finite(result$p_value)]>=0),
            all(result$p_value[is.finite(result$p_value)]<=1),
            any(result$response_id=="gene1" & result$n_nonzero_trt==0))
  cat("RAW_COUNTS_FROZEN_CALLS_COMPLETE_LOSS_AND_CALIBRATION_CHECKPOINT_VERIFIED\n")
}

run <- function(input,output,checkpoint=NULL) {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  stage <- function(value) {
    writeLines(value,file.path(output,"phase.txt"))
    cat(value," ",format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),"\n",sep="");flush.console()
  }
  if(is.null(checkpoint)) {
    stage("read_exact_whole_gene_counts")
    prepared <- read_two_target_counts(input)
    selected <- prepared$cells$role=="non_targeting_control" |
      prepared$cells$target_gene_name %in% c("TYMSOS","UBA1")
    meta <- prepared$cells[selected,,drop=FALSE]
    y <- prepared$counts[,selected,drop=FALSE]
    rm(prepared);gc()
  } else {
    stage("resume_saved_configured_object")
    meta <- read.csv(gzfile(file.path(input,"cells.csv.gz")),stringsAsFactors=FALSE)
    meta <- meta[meta$role=="non_targeting_control" |
      meta$target_gene_name %in% c("TYMSOS","UBA1"),,drop=FALSE]
    object <- readRDS(checkpoint)
    validate_checkpoint_records(object,meta)
  }
  stopifnot(nrow(meta)==20561L,sum(meta$role=="non_targeting_control")==19342L)
  design <- read.csv(file.path(input,"guide_design.csv"),stringsAsFactors=FALSE)
  guide_ids <- sort(unique(meta$guide_id))
  targets <- sort(unique(meta$guide_id[meta$role %in% c("targeting","positive_control")]))
  stopifnot(length(guide_ids)==116L,length(targets)==10L,all(meta$guide_umis>=1))
  if(is.null(checkpoint)) {
    g <- sparseMatrix(i=match(meta$guide_id,guide_ids),j=seq_len(nrow(meta)),x=meta$guide_umis,
      dims=c(length(guide_ids),nrow(meta)),dimnames=list(guide_ids,meta$cell_barcode))
    stage("official_pairwise_QC_and_frozen_assignments")
    object <- configure_object(y,g,meta,targets)
    saveRDS(object,file.path(output,"precalibration_object.rds"))
  }
  if(!object@functs_called[["run_calibration_check"]]) {
    stage("official_NTC_calibration")
    set.seed(10282021)
    object <- run_calibration_check(object,n_calibration_pairs=5000L,calibration_group_size=1L,
      parallel=TRUE,n_processors=2L,print_progress=TRUE,output_amount=2L)
    saveRDS(object,file.path(output,"postcalibration_object.rds"))
  }
  calibration <- get_result(object,"run_calibration_check")
  write_gzip_csv(calibration,file.path(output,"NTC_calibration_tests.csv.gz"))
  nt_guides <- unique(meta$guide_id[meta$role=="non_targeting_control"])
  validate_calibration_controls(calibration,nt_guides)
  finite <- calibration$p_value[is.finite(calibration$p_value)]
  write.csv(data.frame(tests=nrow(calibration),finite_p=length(finite),
    fraction_p_below_05=mean(finite<0.05),fraction_p_below_01=mean(finite<0.01),
    BH_discoveries_at_10=sum(p.adjust(finite,"BH")<0.1)),
    file.path(output,"calibration_summary.csv"),row.names=FALSE)
  if(!object@functs_called[["run_discovery_analysis"]]) {
    stage("official_per_guide_whole_gene_discovery")
    set.seed(10282021)
    object <- run_discovery_analysis(object,parallel=TRUE,n_processors=2L,print_progress=TRUE,output_amount=2L)
    saveRDS(object,file.path(output,"postdiscovery_object.rds"))
  }
  result <- get_result(object,"run_discovery_analysis")
  write_gzip_csv(result,file.path(output,"guide_RNA_tests.csv.gz"))
  guide_column <- if("grna_id" %in% colnames(result)) "grna_id" else "grna_target"
  stopifnot(all(result[[guide_column]] %in% targets))
  summaries <- lapply(targets,function(guide) {
    part <- result[result[[guide_column]]==guide,,drop=FALSE]
    measured <- is.finite(part$p_value)
    gene <- design$target_gene_name[match(guide,design$guide)]
    data.frame(guide_id=guide,target_gene_name=gene,records=sum(meta$guide_id==guide),
      tested_RNAs=nrow(part),finite_RNA_tests=sum(measured),
      DE_RNAs=sum(part$significant %in% TRUE),
      assessment=if(sum(measured)>0L) "evaluated" else "not_estimable")
  })
  write.csv(do.call(rbind,summaries),file.path(output,"guide_DE_counts.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
  stage("sceptre_complete")
}

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=2L,args[[1]] %in% c("--self-test","--run"))
source(args[[2]])
if(args[[1]]=="--self-test") self_test() else {
  stopifnot(length(args) %in% c(4L,5L))
  run(args[[3]],args[[4]],if(length(args)==5L) args[[5]] else NULL)
}
