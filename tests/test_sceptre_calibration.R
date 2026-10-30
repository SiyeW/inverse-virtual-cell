# End-to-end glue check only. This synthetic signal is not a biological control.
args <- commandArgs(trailingOnly=TRUE)
source(args[[1]],local=TRUE)
test_features <- data.frame(response_id=c("R1","R2","R3"),gene_name=c("A","B","C"))
test_design <- data.frame(guide=c("a","b"),target_gene_id=c("R1","R2"),target_gene_name=c("A","B"))
test_pairs <- construct_frozen_discovery(test_features,test_design,c("a","b"))
stopifnot(nrow(test_pairs$pairs)==4L,nrow(test_pairs$excluded)==2L,
  !any(with(test_pairs$pairs,(grna_target=="a" & response_id=="R1") | (grna_target=="b" & response_id=="R2"))))
set.seed(3101)
n <- 1600L; p <- 60L
ids <- paste0("g",seq_len(8));calls <- rep(ids,each=200)
depth <- exp(rnorm(n,0,.25))
y <- matrix(rpois(n*p,rep(depth,each=p)*rep(seq(1,5,length.out=p),n)),nrow=p)
y[1,calls=="g7"] <- rpois(200,50*depth[calls=="g7"])
barcodes <- paste0("fixture-cell-",seq_len(n))
y <- as(Matrix(y,sparse=TRUE),"CsparseMatrix");dimnames(y)<-list(paste0("RNA",seq_len(p)),barcodes)
g <- sparseMatrix(i=match(calls,ids),j=seq_len(n),x=20L,dims=c(8L,n),dimnames=list(ids,barcodes))
meta <- data.frame(cell_barcode=barcodes,guide_id=calls,total_RNA_full=round(8000*depth),
  detected_RNA_full=round(3000*sqrt(depth)+rnorm(n,0,30)))
discovery <- data.frame(grna_target="g7",response_id="RNA1")
for(mechanism in c("permutations","crt")) {
  s <- configure_frozen_object(y,g,meta,ids[1:6],mechanism,discovery)
  set.seed(3102)
  s <- run_calibration_check(s,n_calibration_pairs=100L,calibration_group_size=1L,
    parallel=FALSE,print_progress=FALSE)
  r <- get_result(s,"run_calibration_check")
  check_calibration_references(r,y,meta,ids[1:6])
  stopifnot(mean(r$p_value<.05,na.rm=TRUE)<.15)
  s <- run_discovery_analysis(s,parallel=FALSE,print_progress=FALSE)
  d <- get_result(s,"run_discovery_analysis")
  stopifnot(nrow(d)==1L,d$p_value<1e-4,d$log_2_fold_change>2)
  cat("PASSED: ",mechanism," preserves calls, cells, full NTC reference and detects injected signal\n",sep="")
}
