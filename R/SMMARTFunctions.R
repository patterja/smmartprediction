

tx2gene <- function(path_to_target="/Users/patterja/Workspace/SMMART_DrugPrediction/DrugPrediction_v1.0.0/data/target_id.txt", path_to_kallisto_output="/Volumes/exahead1/lustre1/HeiserLab/patterja/tatlow_SMMART/output"){
  #' tx2gene wrapper
  #' 
  #' Convert kallisto output to gene level
  #' @param path_to_target (string): path to the target_id.txt. file that converts ENSEMBL ids to HGVS. Filtering by trusight genes
  #' @param path_to_kallisto_output (string): path to output of kallisto
  #' @return genetxi
  #' 
  #' @export
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

processRUV <- function(filtmat, includeBCCL=FALSE){
  #' Filter by UHR reps by default. Filtmat must include UHR replicates
  #' 
  #'  For use in normalize RUV functions
  #' @param filtmat (matrix): filtered matrix of sample=columns, genes=rows
  #' @param includeBCCL (logic): TRUE/FALSE
  #' @return ruvset
  #' @import RUVSeq
  #' @export

  filtmat= as.matrix(round(filtmat, digits=0))
  
  #set <- EDASeq::newSeqExpressionSet(as.matrix(filtmat))
  
  if (includeBCCL == TRUE){
    longestrow=length(grep("*UHR*", (colnames(filtmat)), value=FALSE))
    reps=rbind(grep("*UHR*", (colnames(filtmat)), value=FALSE), 
               c(-1,grep("*SkBr3*", (colnames(filtmat)), value=FALSE), 
                 rep(-1,longestrow-length(grep("*SkBr3*", (colnames(filtmat))))-1)),
               c(grep("*T47D_6*", (colnames(filtmat)), value=FALSE), 
                 rep(-1,longestrow-length(grep("*T47D_6*", (colnames(filtmat)))))),
               c(grep("*HCC1954*", (colnames(filtmat)), value=FALSE),
                 rep(-1,longestrow-length(grep("*HCC1954*", (colnames(filtmat)))))))
    
    nonreps = setdiff(c(1:ncol(filtmat)), c(reps[1,], reps[2,which(reps[2,]>0)], 
                                            reps[3,which(reps[3,]>0)],reps[4,which(reps[4,]>0)]))
    } else {
      print("no BCCL")
      longestrow=length(grep("*UHR*", (colnames(filtmat)), value=FALSE))
      reps=rbind(grep("*UHR*", (colnames(filtmat)), value=FALSE))#, 
      nonreps = setdiff(c(1:ncol(filtmat)), c(reps[1,]))
    }
  
  newreps=reps
  idx=nrow(reps)+1
  #idx=5
  
  for (i in nonreps){
    first3 = regmatches(colnames(filtmat)[i], regexpr("^.{0,3}", colnames(filtmat)[i]))
    if (first3 != "X17"){
      newreps=rbind(newreps, c(rep(-1, longestrow-1), i))
      idx=idx+1
    } else {
      print("smmart")
      rgx=paste0(regmatches(colnames(filtmat)[i], regexpr("^([^_]*_[^_]*_[^_]*_[^_]*_[^_]*_)", 
                                                          colnames(filtmat)[i])),".*UHR*")
      j=which(grepl(rgx, colnames(filtmat))) #get index corresponding UHR
      jidx=which(newreps[1,]==j) #get index on the rep matrix
      newreps=rbind(newreps, -1) #add an "empty" row
      newreps[idx,jidx]=i
      idx=idx+1
    }
  }
  
  if (includeBCCL == TRUE){
    setk=7
  } else {
    setk=1
  }
  ruvs_set <- RUVSeq::RUVs(filtmat, cIdx=rownames(filtmat), k=setk, scIdx=newreps)
  return(ruvs_set)

}



normalizeRUV <- function(pathtocnts="/Users/patterja/Workspace/SMMART_HeiserNormalization/data/FINAL_whole.txt", includeBCCL=FALSE, filtperc=10, normtype="CPM") {
  #' normalize RUV
  #' 
  #' Takes a matrix and filters by counts greater than 10 in filtperc % of SMMART samples
  #' @param pathtocnts: matrix of expression, rownames=genes
  #' @param normtype (string): "NONE", "CPM"
  #' @return list of matrices
  #' 
  #' @export

  
  allcnts = read.csv(pathtocnts, sep="\t", row.names = 1)
  smmartsamps = colnames(allcnts[,grepl("NS500642|NS500681", colnames(allcnts))==TRUE])
  smmartsamps = smmartsamps[grepl("RNA170328LH|RNA170329LH", smmartsamps)==FALSE]
  filtprop = filtperc/100
  
  if (includeBCCL == TRUE){
    cnts = allcnts
    } else {
      smmartcnts = allcnts[,smmartsamps]
      cnts = smmartcnts
      }

  if (normtype =="NONE"){
    #get SMMART only
    smmartcnts = cnts[,smmartsamps]
    #FILTER
    perc_smmart = floor(filtprop*(ncol(smmartcnts)))
    gene_smmart = rownames(smmartcnts[apply(data.frame(smmartcnts), 1, 
                                            function(x) length(x[x>5])>=perc_smmart),])
    filt_cnts = cnts[gene_smmart,]
    #ruv
    ruvs = processRUV(filt_cnts, includeBCCL=includeBCCL)
    
  }

  if (normtype =="CPM"){
    #cpm 
    cpm = apply(cnts,2, function(x) (x/sum(x))*1000000)
    #get SMMART only
    smmartcpm = cpm[,smmartsamps]
    #FILTER
    perc_smmart = floor(filtprop*(ncol(smmartcpm)))
    genelistcpm = rownames(smmartcpm[apply(data.frame(smmartcpm), 1, 
                                           function(x) length(x[x>10])>=perc_smmart),])
    filtcpm = cpm[genelistcpm,]
    #ruv cpm
    ruvs = processRUV(filtcpm, includeBCCL=includeBCCL)
    
  }

return(ruvs)
}


updateSMMART_datasets <- function(proj_output="/Volumes/exahead1/lustre1/HeiserLab/patterja/tatlow_SMMART/output",
                         target_filename = "/Users/patterja/Workspace/SMMART_HeiserNormalization/data/target_id.txt") {
  #' updateSMMART_datasets
  #' 
  #' TODO: needs updating from generate_smmart_predictions
  #' @param proj_output (character string) full path to all kallisto output directory
  #' @param target_filename (character string) full path to the target_id.txt file.
  #' @return matrix aggregated to gene gene counts.
  #' 
  #' @export
  target = read.csv(target_filename, sep="\t")
  tx2gene=data.frame(TXNAME=target$TXNAME, GENEID=target$X.1)
  proj_outputdir = proj_output
  fdir=list.files(proj_outputdir)
  files= file.path(proj_outputdir, fdir, "abundance.h5")
  names(files) = fdir
  print(files)
  txi <- tximport(files = files, type = "kallisto", txOut = TRUE)
  genetxi = summarizeToGene(txi, tx2gene = tx2gene)
  write.table(genetxi$counts, "/Users/patterja/Workspace/SMMART_HeiserNormalization/data/smmart_counts_bygene.txt", sep="\t", quote=F, col.names=NA)
  write.table(genetxi$abundance, "/Users/patterja/Workspace/SMMART_HeiserNormalization/data/smmart_abund_bygene.txt", sep="\t", quote=F, col.names=NA)
  write.table(genetxi$abundance, "/Users/patterja/Workspace/SMMART_HeiserNormalization/data/smmart_efflngth_bygene.txt", sep="\t", quote=F, col.names=NA)

  
}


adjust_EffLength <- function(ruvs, includeBCCL=FALSE){
  #deprecated. makes no sense
    efflength = read.csv("/Users/patterja/Box Sync/Heiser_Atwater/SMMART_Normalization/DATA_matrices/efflengths.txt", sep="\t", stringsAsFactors = FALSE, row.names = 1)
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
