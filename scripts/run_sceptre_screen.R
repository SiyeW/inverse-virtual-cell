suppressPackageStartupMessages(library(sceptre))
suppressPackageStartupMessages(library(Matrix))
stopifnot(packageVersion("sceptre")=="0.10.3")
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==5L)
source(args[[1]])
input <- args[[2]];prior <- args[[3]];output <- args[[4]];seed <- as.integer(args[[5]])
dir.create(output,showWarnings=FALSE)
stage <- function(x) {writeLines(x,file.path(output,"stage.txt"));cat("STAGE ",x,"\n",sep="");flush.console()}
stage("reuse_accepted_calibrated_object")
s <- readRDS(prior)
stopifnot(length(s@cells_in_use)==12591L,length(s@grna_assignments$all_nt_idxs)==3769L,
  length(s@grna_assignments$indiv_nt_grna_idxs)==97L,s@functs_called[["run_calibration_check"]])
meta <- read.delim(file.path(input,"cells.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
features <- read.delim(file.path(input,"features.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
design <- read.delim(file.path(input,"guide_design.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(meta)==12591L,nrow(features)==2000L)
nt_guides <- sort(unique(meta$guide_id[meta$role=="non_targeting_control"]))
guides <- sort(unique(meta$guide_id[meta$role %in% c("targeting","positive_control")]))
stopifnot(length(guides)==368L,length(nt_guides)==97L)
stage("complete_NTC_guide_RNA_families")
set.seed(seed)
s <- run_calibration_check(s,n_calibration_pairs=194000L,calibration_group_size=1L,
  parallel=FALSE,print_progress=TRUE,output_amount=2L)
nt <- get_result(s,"run_calibration_check")
stopifnot(!anyDuplicated(nt[,c("grna_target","response_id")]),all(nt$grna_target %in% nt_guides))
write.csv(nt,file.path(output,"NTC_complete_tests.csv"),row.names=FALSE)
nt_summary <- guide_omnibus(nt,nt_guides,setNames(rep(2000L,97L),nt_guides))
write.csv(nt_summary,file.path(output,"NTC_guide_summary.csv"),row.names=FALSE)
p <- nt$p_value[is.finite(nt$p_value)]
gate <- length(p)>=.95*nrow(nt) && mean(p<.05)<=.10 && mean(p<.01)<=.025
write.table(data.frame(proceed=gate,NTC_tests=nrow(nt),fraction_p_lt05=mean(p<.05),fraction_p_lt01=mean(p<.01),
  NTC_guides_passing_same_rule=sum(nt_summary$RNA_association_candidate),
  this_is_not_independent_biological_validation=TRUE),file.path(output,"gate_status.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
if(!gate) {stage("complete_NTC_review_failed_no_targeting_tests");quit(status=0)}
stage("official_discovery_on_predeclared_pairs")
s <- run_discovery_analysis(s,parallel=FALSE,print_progress=TRUE,output_amount=2L)
t <- get_result(s,"run_discovery_analysis")
stopifnot(all(t$grna_target %in% guides),!anyDuplicated(t[,c("grna_target","response_id")]))
write.csv(t,file.path(output,"targeting_RNA_tests.csv"),row.names=FALSE)
families <- table(s@discovery_pairs$grna_target)
summary <- guide_omnibus(t,guides,families)
joint <- rbind(summary,nt_summary)
joint$guide_q_BY <- p.adjust(joint$guide_p,"BY")
joint$RNA_association_candidate <- joint$guide_q_BY<.05
joint$multiplicity_family_guides <- nrow(joint)
summary <- joint[match(guides,joint$guide),,drop=FALSE]
nt_summary <- joint[match(nt_guides,joint$guide),,drop=FALSE]
write.csv(nt_summary,file.path(output,"NTC_guide_summary.csv"),row.names=FALSE)
summary$role <- meta$role[match(summary$guide,meta$guide_id)]
summary$records <- as.integer(table(meta$guide_id)[summary$guide])
summary$target <- design$target_gene_name[match(summary$guide,design$guide)]
write.csv(summary,file.path(output,"guide_screen_summary.csv"),row.names=FALSE)
stage("same_target_guide_comparison")
source(file.path(dirname(args[[1]]),"bind_sceptre_counts.R"))
prepared <- readMM(gzfile(file.path(input,"response_counts.mtx.gz")))
raw <- bind_sceptre_counts(s,prepared,meta,features)
x <- raw %*% Diagonal(x=10000/meta$total_RNA_full);x@x <- log1p(x@x)
all_guides <- sort(unique(meta$guide_id));sizes <- table(meta$guide_id)
membership <- sparseMatrix(i=match(meta$guide_id,all_guides),j=seq_len(nrow(meta)),
  x=1/as.numeric(sizes[meta$guide_id]),dims=c(length(all_guides),nrow(meta)))
means <- as.matrix(x %*% t(membership));colnames(means)<-all_guides
reference <- Matrix::rowMeans(x[,meta$role=="non_targeting_control",drop=FALSE])
effects <- means-reference
retained <- summary$guide[summary$RNA_association_candidate]
pairs <- profile_pairs(effects,meta,features,design,retained)
write.csv(pairs,file.path(output,"same_target_pairs.csv"),row.names=FALSE)
result <- data.frame(records=12591L,NTC_records=3769L,NTC_guides=97L,eligible_targeting_guides=368L,
  RNA_family=2000L,tested_target_RNA_pairs=sum(is.finite(t$p_value)),
  RNA_association_candidates=length(retained),NTC_guides_passing_same_rule=sum(nt_summary$RNA_association_candidate),
  targets_with_two_candidate_guides=length(unique(pairs$target)),candidate_same_target_pairs=nrow(pairs),
  same_target_r_median=if(nrow(pairs))median(pairs$same_target_r,na.rm=TRUE) else NA_real_,
  different_target_r_median=if(nrow(pairs))median(pairs$different_target_r,na.rm=TRUE) else NA_real_,
  functional_KO_validated=FALSE,company_screen_substituted=FALSE)
write.csv(result,file.path(output,"screen_summary.csv"),row.names=FALSE)
stage("complete_pending_scientific_review")
print(result)
