args<-commandArgs(trailingOnly=TRUE);source(args[[1]])
suppressPackageStartupMessages(library(sceptre));suppressPackageStartupMessages(library(Matrix))
prepared<-as(Matrix(matrix(c(1,0,2,3,4,0,0,5,6),3,3),sparse=TRUE),"CsparseMatrix")
dimnames(prepared)<-list(c("R1","R2","R3"),c("cell-a","cell-b","cell-c"))
g<-sparseMatrix(i=c(1,2,1),j=1:3,x=20,dims=c(2,3),dimnames=list(c("guide-a","guide-b"),colnames(prepared)))
meta<-data.frame(cell_barcode=colnames(prepared),guide_id=c("guide-a","guide-b","guide-a"),total_RNA_full=c(100,120,90),detected_RNA_full=c(50,68,45))
cv<-meta[,c("total_RNA_full","detected_RNA_full")]
s<-import_data(prepared,g,data.frame(grna_id=rownames(g),grna_target=c("T1","non-targeting")),moi="low",extra_covariates=cv)
s<-set_analysis_parameters(s,formula_object=~log(total_RNA_full)+log(detected_RNA_full))
s<-assign_grnas(s,method="maximum",min_grna_n_umis_threshold=1L,print_progress=FALSE)
stopifnot(is.null(colnames(sceptre:::get_response_matrix(s))))
features<-data.frame(response_id=rownames(prepared));fixed<-bind_sceptre_counts(s,prepared,meta,features)
stopifnot(identical(colnames(fixed),meta$cell_barcode),identical(fixed@x,prepared@x))
bad<-prepared[,c(2,1,3)]
stopifnot(inherits(try(bind_sceptre_counts(s,bad,meta,features),silent=TRUE),"try-error"))
badmeta<-meta;badmeta$guide_id[[1]]<-"guide-b"
stopifnot(inherits(try(bind_sceptre_counts(s,prepared,badmeta,features),silent=TRUE),"try-error"))
cat("Passed official lost-column-name case, exact count reconstruction, column swap and wrong-call rejection\n")
