# Standard Bonferroni within each fixed RNA family, then BY across guides.
# These are expression-association candidates, never DNA or functional KO calls.
guide_omnibus <- function(tests,guides,family_sizes,alpha=.05) {
  stopifnot(!anyDuplicated(guides),all(guides %in% names(family_sizes)))
  rows <- lapply(guides,function(id) {
    t <- tests[tests$grna_target==id,,drop=FALSE];p <- t$p_value[is.finite(t$p_value)]
    stopifnot(all(p>=0 & p<=1),length(p)<=family_sizes[[id]])
    data.frame(guide=id,RNA_family=family_sizes[[id]],testable_RNAs=length(p),
      minimum_RNA_p=if(length(p))min(p) else 1,
      guide_p=min(1,if(length(p))min(p)*family_sizes[[id]] else 1))
  })
  result <- do.call(rbind,rows);result$guide_q_BY <- p.adjust(result$guide_p,"BY")
  result$RNA_association_candidate <- result$guide_q_BY<alpha
  result
}

profile_pairs <- function(effects,meta,features,design,retained) {
  empty <- data.frame(target=character(),guide_A=character(),guide_B=character(),n_A=integer(),n_B=integer(),
    same_target_r=numeric(),background_guide=character(),background_n=integer(),different_target_r=numeric())
  if(length(retained)<2L)return(empty)
  d <- design[match(retained,design$guide),,drop=FALSE]
  d$target <- ifelse(nzchar(d$target_gene_id),d$target_gene_id,d$target_gene_name)
  sizes <- table(meta$guide_id);roles <- setNames(meta$role,meta$guide_id)
  pearson <- function(a,b) {
    ad <- design[match(c(a,b),design$guide),,drop=FALSE]
    keep <- !(features$response_id %in% ad$target_gene_id | features$gene_name %in% ad$target_gene_name)
    x <- effects[keep,a];y <- effects[keep,b]
    if(length(x)<3 || sd(x)==0 || sd(y)==0)return(NA_real_)
    cor(x,y,method="pearson")
  }
  rows <- list()
  for(target in unique(d$target)) {
    ids <- sort(d$guide[d$target==target]);if(length(ids)<2L)next
    pairs <- combn(ids,2,simplify=FALSE)
    for(pair in pairs) {
      a <- pair[[1]];b <- pair[[2]]
      alternatives <- d$guide[d$target!=target]
      alternatives <- alternatives[roles[alternatives]==roles[[b]]]
      bg <- if(length(alternatives))alternatives[order(abs(log(as.numeric(sizes[alternatives])/as.numeric(sizes[[b]]))),alternatives)][[1]] else NA_character_
      rows[[length(rows)+1L]] <- data.frame(target=target,guide_A=a,guide_B=b,n_A=as.integer(sizes[[a]]),n_B=as.integer(sizes[[b]]),
        same_target_r=pearson(a,b),background_guide=bg,
        background_n=if(is.na(bg))NA_integer_ else as.integer(sizes[[bg]]),
        different_target_r=if(is.na(bg))NA_real_ else pearson(a,bg))
    }
  }
  if(length(rows))do.call(rbind,rows) else empty
}
