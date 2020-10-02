#' tx2gene wrapper
#'
#' Convert kallisto output to gene level
#' @param path_to_target (string): path to the target_id.txt. file that converts ENSEMBL ids to HGVS. Filtering by trusight genes
#' @param path_to_kallisto_output (string): path to output of kallisto
#' @return genetxi
#'
#' @export


tx2gene <- function(path_to_target="target_id.txt", path_to_kallisto_output="/Volumes/exahead1/lustre1/HeiserLab/patterja/tatlow_SMMART/output"){
  target = read.csv(path_to_target, sep="\t", stringsAsFactors = FALSE)
  tx_to_gene=tibble(TXNAME=target$TXNAME, GENEID=target$X.1)

  fdir=list.files(path_to_kallisto_output)
  files= file.path(path_to_kallisto_output, fdir, "abundance.h5")
  #filter here if necessary
  #files=files[2:17]
  print(paste0("Processing:\n",files))
  names(files) = fdir
  print(paste0("Processing Kallisto output: ", files))
  txi <- tximport(files = files, type = "kallisto", txOut = TRUE)
  genetxi = summarizeToGene(txi, tx2gene = tx_to_gene)
  return(genetxi)

}
