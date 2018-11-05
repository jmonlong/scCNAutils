##' @title Merge expression of multiple samples
##' @param ge_list a list of ge_df (e.g. read from \code{\link{read_mtx}}).
##' @param sample_names the names of each sample. If NULL, tries to use ge_list's names.
##' @return a list with
##' \item{ge}{the merged gene expression data.frame}
##' \item{info}{a data.frame with new and original cell names, and
##' corresponding sample name}
##' @author Jean Monlong
##' @export
merge_samples <- function(ge_list, sample_names=NULL){
  if(is.null(sample_names)){
    sample_names = names(ge_list)
  }
  if(is.null(sample_names)){
    warning('No sample names provided, will use s1, s2, ... etc')
    sample_names = paste0('s', 1:length(ge_list))
  }
  names(ge_list) = sample_names
  ## Prepare first sample
  ge_df = ge_list[[1]]
  orig.cells = colnames(ge_df)[-1]
  new.cells = paste0(sample_names[1], '_', orig.cells)
  colnames(ge_df)[-1] = new.cells
  info.df = data.frame(cell=new.cells, orig.cell=orig.cells,
                       sample=sample_names[1], stringsAsFactors=FALSE)
  ge_df = data.table::as.data.table(ge_df)
  ## Merge other samples, one at a time
  for(samp in sample_names[-1]){
    ge.tmp = ge_list[[samp]]
    orig.cells = colnames(ge.tmp)[-1]
    new.cells = paste0(samp, '_', orig.cells)
    colnames(ge.tmp)[-1] = new.cells
    ge.tmp = data.table::as.data.table(ge.tmp)
    ge_df = merge(ge_df, ge.tmp, all=TRUE)
    info.df = rbind(info.df,
                    data.frame(cell=new.cells, orig.cell=orig.cells,
                               sample=samp, stringsAsFactors=FALSE))
  }
  ge_df = as.data.frame(ge_df)
  ## Change NAs to 0
  if(any(is.na(ge_df))){  
    ge_df[is.na(ge_df)] = 0
  }
  return(list(ge=ge_df, info=info.df))
}
