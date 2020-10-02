#' processRUV 
#'
#' Filter by UHR reps by default. Filtmat must include UHR replicates
#'  For use in normalize RUV functions
#' @param filtmat (matrix): filtered matrix of sample=columns, genes=rows
#' @param includeBCCL (logic): TRUE/FALSE
#' @param setk (numeric): number of factors to remove.
#' if includeBCCL=T use CPM + setk=7
#' if includeBCCL=F use TPM + setk=1
#' @return ruvset
#' @import RUVSeq
#' @export

processRUV <- function(filtmat, includeBCCL=TRUE, setk){
  filtmat= as.matrix(round(filtmat, digits=0))

  #set <- EDASeq::newSeqExpressionSet(as.matrix(filtmat))

  if (includeBCCL == TRUE){
    replicate_matx = colnames(filtmat)
    replicate_matx[grep("*UHR*", (colnames(filtmat)), value=FALSE)] = "UHR"
    replicate_matx[grep("*SkBr3*", (colnames(filtmat)), value=FALSE)] = "SkBr3"
    replicate_matx[grep("T47D_", (colnames(filtmat)), value=FALSE)]="T47D"
    replicate_matx[grep("*HCC1954*", (colnames(filtmat)), value=FALSE)]="HCC1954"
  } else {
    print("no BCCL")
    replicate_matx = colnames(filtmat)
    replicate_matx[grep("*UHR*", (colnames(filtmat)), value=FALSE)] = "UHR"
  }
  scIdx = makeGroups(replicate_matx)
  dim(scIdx)
  ruvs_set <- RUVSeq::RUVs(filtmat, cIdx=rownames(filtmat), k=setk, scIdx=scIdx)
  return(ruvs_set)

}
