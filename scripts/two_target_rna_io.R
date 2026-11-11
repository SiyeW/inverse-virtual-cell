# Read the exact integer count vectors exported by prepare_two_target_rna.py.
write_gzip_csv <- function(value,path) {
  stream <- gzfile(path,"wt")
  on.exit(close(stream))
  write.csv(value,stream,row.names=FALSE)
}

read_two_target_counts <- function(input, seurat_keys=FALSE) {
  suppressPackageStartupMessages(library(Matrix))
  dimensions <- read.delim(file.path(input,"matrix_dimensions.tsv"))
  cells <- read.csv(gzfile(file.path(input,"cells.csv.gz")),check.names=FALSE,stringsAsFactors=FALSE)
  genes <- read.csv(gzfile(file.path(input,"genes.csv.gz")),check.names=FALSE,stringsAsFactors=FALSE)
  stopifnot(nrow(dimensions)==1L,nrow(cells)==dimensions$cells,nrow(genes)==dimensions$genes,
            !anyDuplicated(cells$cell_barcode),!anyDuplicated(genes$gene_id))
  read_vector <- function(name, number) {
    stream <- gzfile(file.path(input,name),"rb")
    on.exit(close(stream))
    result <- readBin(stream,integer(),n=number,size=4L,endian="little")
    stopifnot(length(result)==number,length(readBin(stream,integer(),n=1,size=4L))==0L)
    result
  }
  values <- read_vector("counts.data.i32.gz",dimensions$nnz)
  indices <- read_vector("counts.indices.i32.gz",dimensions$nnz)
  pointers <- read_vector("counts.indptr.i32.gz",dimensions$cells+1L)
  stopifnot(all(values>=0),pointers[[1]]==0L,tail(pointers,1)==dimensions$nnz)
  matrix <- new("dgCMatrix",i=indices,p=pointers,x=as.numeric(values),
                Dim=as.integer(c(dimensions$genes,dimensions$cells)))
  rm(values,indices,pointers)
  rownames(matrix) <- if(seurat_keys) genes$seurat_key else genes$gene_id
  colnames(matrix) <- cells$cell_barcode
  stopifnot(!anyDuplicated(rownames(matrix)),isTRUE(validObject(matrix)),
            identical(as.numeric(Matrix::colSums(matrix)),as.numeric(cells$total_counts)))
  for(column in intersect(c("guide_id","target_gene_id","target_gene_name"),colnames(cells))) {
    cells[[column]][is.na(cells[[column]])] <- ""
  }
  cells$pilot_record <- as.logical(cells$pilot_record)
  stopifnot(!anyNA(cells$pilot_record))
  list(counts=matrix,cells=cells,genes=genes)
}
