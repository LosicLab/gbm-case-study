# functions for analysis


extract.vaf <- function(var.table, AD.col='gt_AD') {
  var.table <- tidyr::separate(var.table, AD.col, c("AD_N", "AD_T"), sep=",")
  var.table[, c("AD_N", "AD_T")] <- lapply(var.table[, c("AD_N", "AD_T")], as.numeric)
  var.table$VAF_AD <- 100 * ( var.table$AD_T / (var.table$AD_N + var.table$AD_T) )
  return(var.table)
}


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

parse_extra <- function(df, annot){
  res <- substr(x=df$Extra, start=regexpr(pattern=paste(annot, '=', sep=''),
                                          text=df$Extra), stop=nchar(df$Extra) )
  res <- substr(x=res,  start=1, stop=regexpr(pattern=';', text=res))
  res <- gsub(x=gsub(pattern=paste(annot, '=', sep=''), replacement='', x=res),
              pattern=';', replacement='')
  res[!grepl(annot, df$Extra)] <- 'NA'
  return(res)
}



# GSEA lister
# inputs: 
#   rank_list - named list of genes (deg_ranks above). names should correspond to sample comparison in the DESeq2 test.
#   sig_list - gene signature being tested
# output: a dataframe (tibble) containing the fgsea results across all sample comparisons for the tested gene signatures
gsea_lister <- function(sig_list, rank_list, sig_name){
  require(fgsea)
  require(tidyverse)
  # perform the gsea, include the sample comparison
  y <- lapply(names(rank_list), function(i){
    
    x <- fgsea(pathways=sig_list, stats=rank_list[[i]], nperm=1000) %>%
      as_tibble() %>%
      dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
      mutate(sample_comparison=i, signature_name=sig_name)
    return(x)
  })
  
  y <-  do.call(rbind, y) %>%
    as_tibble()
  
  return(y)
}