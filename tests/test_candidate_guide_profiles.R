source(commandArgs(trailingOnly=TRUE)[[1]])
effects<-cbind(a=c(1,2,4,8,16),b=c(2,4,8,16,32),c=c(-1,-2,-4,-8,-16),d=c(1,2,4,8,16),e=1:5)
info<-data.frame(guide=letters[1:5],target_id=c("T1","T1","T2","T3","T4"),target=c("A","A","B","C","D"),
  passing=c(TRUE,FALSE,FALSE,FALSE,TRUE),role="targeting",records=c(30,20,22,100,15))
features<-data.frame(response_id=paste0("R",1:5),gene_name=paste0("X",1:5))
design<-data.frame(guide=letters[1:5],target_gene_id=info$target_id,target_gene_name=info$target)
r<-compare_anchor_profiles(effects,info,features,design)
stopifnot(nrow(r$pairs)==1L,r$pairs$other_same_target=="b",!r$pairs$other_passed_screen,
  abs(r$pairs$same_target_r-1)<1e-12,!"d" %in% r$background$different_target_guide,
  "c" %in% r$background$different_target_guide,!any(r$background$different_target_guide %in% c("a","b")))
cat("Passed non-significant peer inclusion, correct profiles, quantity-matched other genes and single-guide omission\n")
