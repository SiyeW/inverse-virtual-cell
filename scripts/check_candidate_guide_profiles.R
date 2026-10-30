# A guide failing a DE cutoff can still be used to check a same-target pattern.
# All eligible guides are retained here; RNA significance tests are not rerun.
compare_anchor_profiles <- function(effects,info,features,design,max_background=20L) {
  rows<-list();background<-list()
  for(a in info$guide[info$passing]) {
    ai<-match(a,info$guide);peers<-info$guide[info$target_id==info$target_id[[ai]] & info$guide!=a]
    for(b in sort(peers)) {
      bi<-match(b,info$guide)
      if(info$passing[[bi]] && b<a)next
      allowed<-which(info$target_id!=info$target_id[[ai]] & info$role==info$role[[bi]] &
        info$records>=info$records[[bi]]/2 & info$records<=2*info$records[[bi]])
      order_by<-order(abs(log(info$records[allowed]/info$records[[bi]])),info$guide[allowed])
      bg<-info$guide[head(allowed[order_by],max_background)]
      exclude<-design[match(c(a,b,bg),design$guide),,drop=FALSE]
      keep<-!(features$response_id %in% exclude$target_gene_id | features$gene_name %in% exclude$target_gene_name)
      corr<-function(c) {
        x<-effects[keep,a];y<-effects[keep,c]
        if(length(x)<3 || sd(x)==0 || sd(y)==0)NA_real_ else cor(x,y,method="pearson")
      }
      r<-corr(b);br<-vapply(bg,corr,numeric(1));finite<-is.finite(br)
      rows[[length(rows)+1L]]<-data.frame(target=info$target[[ai]],anchor=a,other_same_target=b,
        anchor_records=info$records[[ai]],other_records=info$records[[bi]],other_passed_screen=info$passing[[bi]],
        RNAs_compared=sum(keep),same_target_r=r,background_guides=sum(finite),
        different_target_median_r=if(any(finite))median(br[finite]) else NA_real_,
        fraction_background_below=if(any(finite))mean(br[finite]<r) else NA_real_)
      if(length(bg))background[[length(background)+1L]]<-data.frame(anchor=a,other_same_target=b,
        different_target_guide=bg,different_target_r=br,
        records=info$records[match(bg,info$guide)])
    }
  }
  empty<-data.frame(target=character(),anchor=character(),other_same_target=character(),anchor_records=integer(),
    other_records=integer(),other_passed_screen=logical(),RNAs_compared=integer(),same_target_r=numeric(),
    background_guides=integer(),different_target_median_r=numeric(),fraction_background_below=numeric())
  list(pairs=if(length(rows))do.call(rbind,rows) else empty,
       background=if(length(background))do.call(rbind,background) else data.frame(anchor=character()))
}

args<-commandArgs(trailingOnly=TRUE)
if(length(args)>0 && args[[1]]=="--run") {
  stopifnot(length(args)==4L)
  suppressPackageStartupMessages(library(Matrix))
  input<-args[[2]];screen<-args[[3]];output<-args[[4]];dir.create(output,showWarnings=FALSE)
  meta<-read.delim(file.path(input,"cells.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
  features<-read.delim(file.path(input,"features.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
  design<-read.delim(file.path(input,"guide_design.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
  info<-read.csv(file.path(screen,"guide_screen_summary.csv"),stringsAsFactors=FALSE)
  stopifnot(nrow(meta)==12591L,nrow(features)==2000L,nrow(info)==368L,sum(info$RNA_association_candidate)==6L)
  raw<-readMM(gzfile(file.path(input,"response_counts.mtx.gz")))
  stopifnot(identical(dim(raw),c(2000L,12591L)))
  x<-raw %*% Diagonal(x=10000/meta$total_RNA_full);x@x<-log1p(x@x)
  ids<-sort(unique(meta$guide_id));sizes<-table(meta$guide_id)
  membership<-sparseMatrix(i=match(meta$guide_id,ids),j=seq_len(nrow(meta)),
    x=1/as.numeric(sizes[meta$guide_id]),dims=c(length(ids),nrow(meta)))
  means<-as.matrix(x %*% t(membership));colnames(means)<-ids
  reference<-Matrix::rowMeans(x[,meta$role=="non_targeting_control",drop=FALSE]);effects<-means-reference
  saveRDS(list(effects=effects,guide=ids,RNA=features$response_id,reference=reference),file.path(output,"profile_checkpoint.rds"))
  info$passing<-info$RNA_association_candidate
  id<-match(info$guide,design$guide);info$target_id<-design$target_gene_id[id]
  stopifnot(!anyNA(id),all(info$records==as.integer(sizes[info$guide])),sum(meta$role=="non_targeting_control")==3769L)
  out<-compare_anchor_profiles(effects,info,features,design)
  write.csv(out$pairs,file.path(output,"candidate_same_target_pairs.csv"),row.names=FALSE)
  write.csv(out$background,file.path(output,"different_target_background.csv"),row.names=FALSE)
  selected<-unique(info$target_id[info$passing])
  genes<-do.call(rbind,lapply(selected,function(g) {
    names<-info$target[info$target_id==g];pairs<-out$pairs[out$pairs$target==names[[1]],,drop=FALSE]
    data.frame(target_id=g,target=names[[1]],eligible_guides=sum(info$target_id==g),
      passing_guides=sum(info$target_id==g & info$passing),compared_pairs=nrow(pairs),
      status=if(nrow(pairs))"descriptive_comparison_available" else "no_second_eligible_guide")
  }))
  write.csv(genes,file.path(output,"candidate_target_coverage.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
  print(out$pairs);print(genes)
}
