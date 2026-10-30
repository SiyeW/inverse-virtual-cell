# Resume summaries and profiles from completed tests; never rerun associations.
suppressPackageStartupMessages(library(sceptre))
suppressPackageStartupMessages(library(Matrix))
stopifnot(packageVersion("sceptre")=="0.10.3")
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==6L)
source(args[[1]]);source(file.path(dirname(args[[1]]),"bind_sceptre_counts.R"))
input<-args[[2]];prior<-args[[3]];output<-args[[4]];completed<-args[[6]]
dir.create(output,showWarnings=FALSE)
stage <- function(x){writeLines(x,file.path(output,"stage.txt"));cat("STAGE ",x,"\n",sep="");flush.console()}
stage("verify_saved_complete_expression_results")
s<-readRDS(prior)
meta<-read.delim(file.path(input,"cells.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
features<-read.delim(file.path(input,"features.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
design<-read.delim(file.path(input,"guide_design.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(meta)==12591L,nrow(features)==2000L,length(s@cells_in_use)==12591L)
t<-read.csv(file.path(completed,"targeting_RNA_tests.csv"),stringsAsFactors=FALSE,check.names=FALSE)
nt<-read.csv(file.path(completed,"NTC_complete_tests.csv"),stringsAsFactors=FALSE,check.names=FALSE)
key<-function(d)paste(d$grna_target,d$response_id,sep="\t")
stopifnot(!anyDuplicated(key(t)),nrow(t)==nrow(s@discovery_pairs),setequal(key(t),key(s@discovery_pairs)))
guides<-sort(unique(meta$guide_id[meta$role %in% c("targeting","positive_control")]))
nt_guides<-sort(unique(meta$guide_id[meta$role=="non_targeting_control"]))
stopifnot(length(guides)==368L,length(nt_guides)==97L,!anyDuplicated(key(nt)),all(nt$grna_target %in% nt_guides),nrow(nt)==48181L)
gate<-read.delim(file.path(completed,"gate_status.tsv"));stopifnot(isTRUE(gate$proceed))
summary<-guide_omnibus(t,guides,table(s@discovery_pairs$grna_target))
nt_summary<-guide_omnibus(nt,nt_guides,setNames(rep(2000L,97L),nt_guides))
joint<-rbind(summary,nt_summary);joint$guide_q_BY<-p.adjust(joint$guide_p,"BY")
joint$RNA_association_candidate<-joint$guide_q_BY<.05;joint$multiplicity_family_guides<-nrow(joint)
summary<-joint[match(guides,joint$guide),,drop=FALSE];nt_summary<-joint[match(nt_guides,joint$guide),,drop=FALSE]
old<-read.csv(file.path(completed,"guide_screen_summary.csv"))
stopifnot(identical(old$guide,summary$guide),max(abs(old$guide_q_BY-summary$guide_q_BY))<1e-12,
  identical(old$RNA_association_candidate,summary$RNA_association_candidate))
summary$role<-meta$role[match(summary$guide,meta$guide_id)]
summary$records<-as.integer(table(meta$guide_id)[summary$guide]);summary$target<-design$target_gene_name[match(summary$guide,design$guide)]
stage("verify_cached_count_axes_and_bind_names")
prepared<-readMM(gzfile(file.path(input,"response_counts.mtx.gz")))
raw<-bind_sceptre_counts(s,prepared,meta,features)
write.table(data.frame(integer_nonzero_counts_equal=TRUE,ordered_covariates_equal=TRUE,frozen_calls_equal=TRUE,
  cells=12591L,RNAs=2000L,associations_retested=FALSE,original_counts_or_BAM_read=FALSE),
  file.path(output,"count_binding_proof.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.csv(summary,file.path(output,"guide_screen_summary.csv"),row.names=FALSE)
write.csv(nt_summary,file.path(output,"NTC_guide_summary.csv"),row.names=FALSE)
write.table(gate,file.path(output,"gate_status.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
stage("same_target_guide_comparison_from_saved_selection")
x<-raw %*% Diagonal(x=10000/meta$total_RNA_full);x@x<-log1p(x@x)
all_guides<-sort(unique(meta$guide_id));sizes<-table(meta$guide_id)
membership<-sparseMatrix(i=match(meta$guide_id,all_guides),j=seq_len(nrow(meta)),x=1/as.numeric(sizes[meta$guide_id]),dims=c(length(all_guides),nrow(meta)))
means<-as.matrix(x %*% t(membership));colnames(means)<-all_guides
reference<-Matrix::rowMeans(x[,meta$role=="non_targeting_control",drop=FALSE]);effects<-means-reference
retained<-summary$guide[summary$RNA_association_candidate];pairs<-profile_pairs(effects,meta,features,design,retained)
write.csv(pairs,file.path(output,"same_target_pairs.csv"),row.names=FALSE)
result<-data.frame(records=12591L,NTC_records=3769L,NTC_guides=97L,eligible_targeting_guides=368L,RNA_family=2000L,
  tested_target_RNA_pairs=sum(is.finite(t$p_value)),RNA_association_candidates=length(retained),
  NTC_guides_passing_same_rule=sum(nt_summary$RNA_association_candidate),targets_with_two_candidate_guides=length(unique(pairs$target)),
  candidate_same_target_pairs=nrow(pairs),same_target_r_median=if(nrow(pairs))median(pairs$same_target_r,na.rm=TRUE) else NA_real_,
  different_target_r_median=if(nrow(pairs))median(pairs$different_target_r,na.rm=TRUE) else NA_real_,functional_KO_validated=FALSE,company_screen_substituted=FALSE)
write.csv(result,file.path(output,"screen_summary.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
stage("complete_pending_scientific_review");print(result)
