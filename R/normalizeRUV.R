#' normalize RUV
#'
#' Takes a matrix and filters by counts greater than 10 in filtperc \% of SMMART samples
#' @param exp_matx: matrix of expression, rownames=genes
#' @param includeBCCL: TRUE OR FALSE
#' @param normtype (string): "NONE", "CPM"
#' @param filtperc (int): 10\% of SMMART samples will have counts greater than mincnt
#' @param mincnt (int): TPM use 3, counts to CPM use 10
#' @param setk (numeric): number of factors to remove.
#' if includeBCCL=T use CPM, setk=7, filtperc=10, mincnt = 10,
#' if includeBCCL=F use TPM, setk=1, filtperc=10, mincnt = 3,
#' @return list of matrices
#'
#' @export

normalizeRUV <- function(exp_matx=read.csv("/Users/patterja/Box\ Sync/Heiser_Atwater/SMMART_NORMALIZATION/DATA_MATRICES/PROCESSED/tpm_breast_20190214.txt", sep="\t", row.names = 1),
                         includeBCCL=FALSE, filtperc=10, mincnt=3, normtype="NONE", setk) {

  #smmart samples differentiated by sequencer name
  smmartsamps = colnames(exp_matx)[!grepl("ACXX|AAXX", colnames(exp_matx))]
  smmartsamps = smmartsamps[!grepl("RNA170328LH|RNA170329LH|^X170322_NS500642", smmartsamps)]
  filtprop = filtperc/100

  if (includeBCCL == TRUE){
    cnts = exp_matx
    } else {
      smmartcnts = exp_matx[,smmartsamps]
      cnts = smmartcnts
      }

  if (normtype =="NONE"){
    #get SMMART only
    smmartcnts = cnts[,smmartsamps]
    #FILTER
    perc_smmart = floor(filtprop*(ncol(smmartcnts)))
    perc_smmart = ifelse(perc_smmart == 0, yes=1, no=perc_smmart)
    gene_smmart = rownames(smmartcnts[apply(data.frame(smmartcnts), 1,
                                            function(x) length(x[x>mincnt])>=perc_smmart),])
    filt_cnts = cnts[gene_smmart,]
    #ruv
    ruvs = processRUV(filt_cnts, includeBCCL=includeBCCL, setk=setk)

  }

  if (normtype =="CPM"){
    #cpm
    cpm = apply(cnts,2, function(x) (x/sum(x))*1000000)
    #get SMMART only
    smmartcpm = cpm[,smmartsamps]
    #FILTER
    perc_smmart = floor(filtprop*(ncol(smmartcpm)))
    genelistcpm = rownames(smmartcpm[apply(data.frame(smmartcpm), 1,
                                           function(x) length(x[x>mincnt])>=perc_smmart),])
    filtcpm = cpm[genelistcpm,]
    #ruv cpm
    ruvs = processRUV(filtcpm, includeBCCL=includeBCCL, setk=setk)

  }

return(ruvs)
}

