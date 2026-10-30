# SCEPTRE 0.10.3 deliberately drops column names when creating row-sparse data.
# Restore names only after exact count, covariate-order and frozen-call checks.
bind_sceptre_counts <- function(s,prepared,meta,features) {
  cached <- sceptre:::get_response_matrix(s)
  stopifnot(identical(dim(cached),dim(prepared)),ncol(cached)==nrow(meta),nrow(cached)==nrow(features))
  stopifnot(identical(rownames(cached),features$response_id))
  if(!is.null(colnames(cached)))stopifnot(identical(colnames(cached),meta$cell_barcode))
  a <- as(cached,"CsparseMatrix");b <- as(prepared,"CsparseMatrix")
  stopifnot(identical(a@p,b@p),identical(a@i,b@i),identical(a@x,b@x))
  cv <- s@covariate_data_frame
  stopifnot(identical(as.numeric(cv$total_RNA_full),as.numeric(meta$total_RNA_full)),
    identical(as.numeric(cv$detected_RNA_full),as.numeric(meta$detected_RNA_full)))
  called <- rep(NA_character_,nrow(meta))
  for(id in names(s@initial_grna_assignment_list))called[s@initial_grna_assignment_list[[id]]] <- id
  stopifnot(identical(called,meta$guide_id))
  dimnames(b) <- list(features$response_id,meta$cell_barcode)
  b
}
