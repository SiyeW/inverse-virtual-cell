# Official SCEPTRE v0.10.3. This runner tests real NTCs before target discovery.
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(sceptre))
stopifnot(packageVersion("sceptre") == "0.10.3")

construct_frozen_discovery <- function(features,design,candidate_guides) {
  stopifnot(all(candidate_guides %in% design$guide),!anyDuplicated(design$guide))
  pairs <- expand.grid(grna_target=candidate_guides,response_id=features$response_id,stringsAsFactors=FALSE)
  intended <- design[match(candidate_guides,design$guide),,drop=FALSE]
  excluded <- lapply(seq_along(candidate_guides),function(i) {
    matched <- features$response_id == intended$target_gene_id[[i]] |
      features$gene_name == intended$target_gene_name[[i]]
    data.frame(grna_target=rep(candidate_guides[[i]],sum(matched)),response_id=features$response_id[matched])
  })
  excluded <- do.call(rbind,excluded)
  if(nrow(excluded)>0) pairs <- pairs[!paste(pairs$grna_target,pairs$response_id,sep="\t") %in%
    paste(excluded$grna_target,excluded$response_id,sep="\t"),,drop=FALSE]
  list(pairs=pairs,excluded=excluded)
}

configure_frozen_object <- function(y, g, meta, nt_guides, mechanism,
                                    discovery=data.frame(grna_target=character(),response_id=character())) {
  stopifnot(identical(colnames(y), meta$cell_barcode), identical(colnames(g), meta$cell_barcode))
  stopifnot(all(Matrix::colSums(g > 0) == 1), all(meta$total_RNA_full > 0), all(meta$detected_RNA_full > 0))
  mapping <- data.frame(grna_id=rownames(g), grna_target=ifelse(rownames(g) %in% nt_guides,"non-targeting",rownames(g)))
  covariates <- meta[,c("total_RNA_full","detected_RNA_full"),drop=FALSE]
  rownames(covariates) <- meta$cell_barcode
  s <- import_data(y,g,mapping,moi="low",extra_covariates=covariates)
  s <- set_analysis_parameters(s,discovery_pairs=discovery,side="both",control_group="nt_cells",
    formula_object=~log(total_RNA_full)+log(detected_RNA_full),
    resampling_mechanism=mechanism,resampling_approximation="skew_normal",multiple_testing_alpha=.05)
  # One frozen call per record: this preserves the accepted caller. It does not
  # re-estimate exposure from a complete raw guide-capture matrix.
  s <- assign_grnas(s,method="maximum",min_grna_n_umis_threshold=1L,umi_fraction_threshold=.8,print_progress=FALSE)
  called <- rep(NA_character_,ncol(y))
  for (id in names(s@initial_grna_assignment_list)) called[s@initial_grna_assignment_list[[id]]] <- id
  stopifnot(identical(called,meta$guide_id))
  s <- run_qc(s,n_nonzero_trt_thresh=7L,n_nonzero_cntrl_thresh=7L,
    response_n_umis_range=c(0,1),response_n_nonzero_range=c(0,1),p_mito_threshold=1)
  stopifnot(identical(s@cells_in_use,seq_len(ncol(y))))
  actual_nt <- s@grna_assignments$all_nt_idxs
  expected_nt <- which(meta$guide_id %in% nt_guides)
  stopifnot(setequal(actual_nt,expected_nt))
  stopifnot(setequal(names(s@grna_assignments$indiv_nt_grna_idxs),nt_guides))
  stopifnot(length(unique(unlist(s@grna_assignments$indiv_nt_grna_idxs)))==length(expected_nt))
  s
}

check_calibration_references <- function(result, y, meta, nt_guides) {
  stopifnot(all(result$grna_target %in% nt_guides))
  nt_positive <- Matrix::rowSums(y[,meta$guide_id %in% nt_guides,drop=FALSE] > 0)
  expected <- unname(nt_positive[match(result$response_id,rownames(y))])
  # An evaluated NTC guide is the treatment; every other NTC guide is control.
  stopifnot(all(result$n_nonzero_trt+result$n_nonzero_cntrl == expected))
  invisible(TRUE)
}

run_real_calibration <- function(input, output, mechanism, pairs, seed,
                                expected=c(71664L,19342L,106L,678L)) {
  dir.create(output,recursive=FALSE,showWarnings=FALSE)
  stage <- function(x) {
    cat("STAGE ",x," ",format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),"\n",sep="");flush.console()
    writeLines(x,file.path(output,"stage.txt"))
  }
  stage("read_verified_count_subset")
  meta <- read.delim(file.path(input,"cells.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
  features <- read.delim(file.path(input,"features.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
  design <- read.delim(file.path(input,"guide_design.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
  y <- as(readMM(gzfile(file.path(input,"response_counts.mtx.gz"))),"CsparseMatrix")
  g <- as(readMM(gzfile(file.path(input,"assigned_guide_counts.mtx.gz"))),"CsparseMatrix")
  stopifnot(length(expected)==4L,nrow(meta)==expected[[1]],nrow(features)==2000L,nrow(design)==1196L)
  stopifnot(identical(dim(y),c(2000L,expected[[1]])),identical(dim(g),c(1196L,expected[[1]])))
  dimnames(y) <- list(features$response_id,meta$cell_barcode)
  dimnames(g) <- list(design$guide,meta$cell_barcode)
  stopifnot(!anyDuplicated(rownames(y)),!anyDuplicated(colnames(y)),all(y@x==round(y@x)))
  nt_guides <- sort(unique(meta$guide_id[meta$role=="non_targeting_control"]))
  stopifnot(length(nt_guides)==expected[[3]],sum(meta$role=="non_targeting_control")==expected[[2]])
  stopifnot(setequal(unique(meta$role),c("targeting","positive_control","non_targeting_control","intergenic_control")))
  candidate_guides <- sort(unique(meta$guide_id[meta$role %in% c("targeting","positive_control")]))
  # Role labels are those of the frozen design, not inferred biological hits.
  stopifnot(length(candidate_guides)==expected[[4]])
  discovery <- construct_frozen_discovery(features,design,candidate_guides)
  write.csv(discovery$excluded,file.path(output,"excluded_direct_target_pairs.csv"),row.names=FALSE)
  stage("configure_frozen_calls_and_full_RNA_covariates")
  s <- configure_frozen_object(y,g,meta,nt_guides,mechanism,discovery$pairs)
  write.table(data.frame(records=length(s@cells_in_use),NTC_records=length(s@grna_assignments$all_nt_idxs),
    NTC_guides=length(s@grna_assignments$indiv_nt_grna_idxs),response_RNAs=nrow(y),pair_min_expressing_records=7L,
    direct_target_RNA_decrease_filter=FALSE,predeclared_candidate_guides=length(candidate_guides),
    discovery_pairs_after_direct_target_exclusion=nrow(discovery$pairs),
    pairwise_QC_eligible_discovery_pairs=s@n_ok_discovery_pairs),file.path(output,"cohort_contract.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  saveRDS(s,file.path(output,"precalibration_object.rds"))
  stage("official_NTC_calibration")
  set.seed(seed)
  s <- run_calibration_check(s,n_calibration_pairs=pairs,calibration_group_size=1L,
    parallel=FALSE,print_progress=TRUE,output_amount=2L)
  result <- get_result(s,"run_calibration_check")
  check_calibration_references(result,y,meta,nt_guides)
  write.csv(result,file.path(output,"NTC_expression_tests.csv"),row.names=FALSE)
  p <- result$p_value[is.finite(result$p_value)]
  stopifnot(length(p)>0,all(p>=0 & p<=1))
  summary <- data.frame(method="official_sceptre_v0.10.3",resampling_mechanism=mechanism,
    requested_pairs=pairs,tested_pairs=nrow(result),finite_p=length(p),
    represented_NTC_guides=length(unique(result$grna_target)),represented_RNAs=length(unique(result$response_id)),
    nominal_p_below_05=sum(p<.05),fraction_p_below_05=mean(p<.05),
    nominal_p_below_01=sum(p<.01),fraction_p_below_01=mean(p<.01),
    BH_below_05=sum(p.adjust(p,"BH")<.05),minimum_p=min(p),
    calibration_review_passed=(length(p)>=.95*nrow(result) && mean(p<.05)<=.10 && mean(p<.01)<=.025),
    biological_positive_control_power_assessed=FALSE,targeting_tests_performed=FALSE)
  write.csv(summary,file.path(output,"calibration_summary.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
  saveRDS(s,file.path(output,"calibrated_object.rds"))
  stage("complete_pending_scientific_review")
  print(summary)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0 && args[[1]]=="--run") {
  stopifnot(length(args) %in% c(6,7),args[[4]] %in% c("permutations","crt"))
  expected <- if(length(args)==7) as.integer(strsplit(args[[7]],",",fixed=TRUE)[[1]]) else c(71664L,19342L,106L,678L)
  run_real_calibration(args[[2]],args[[3]],args[[4]],as.integer(args[[5]]),as.integer(args[[6]]),expected)
}
