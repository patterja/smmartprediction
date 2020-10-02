#' adjust_EffLength
#'
#' Adjusts effective length measurements 
#' @param ruvs (data frame) rows are gene name
#' @param efflength_pathdir (string) efflengths.txt
#' @return dataframe
#'
#' @export

adjust_EffLength <- function(ruvs, efflength_pathdir="/Users/patterja/Box/Heiser_Atwater/SMMART_Normalization/DATA_matrices/efflengths.txt", includeBCCL=FALSE){
  #deprecated. makes no sense
    efflength = read.csv(efflength_pathdir, sep="\t", stringsAsFactors = FALSE, row.names = 1)
  if (includeBCCL == TRUE){
    filt_efflength = efflength[rownames(ruvs),]
    reordruvs = ruvs[rownames(filt_efflength), colnames(filt_efflength)]
  } else {
    filt_efflength = efflength[rownames(ruvs), grepl("NS500681|NS500642", colnames(efflength))]
    filt_efflength = filt_efflength[, !grepl("X170424_NS500681", colnames(filt_efflength))]
    reordruvs = ruvs[rownames(filt_efflength), colnames(filt_efflength)]
  }
  adjruvs = reordruvs/filt_efflength
  return(adjruvs)

}
