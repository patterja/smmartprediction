#Functions

theme_jsnick <- function () { 
  #' ggplot2 themejsnick
  #' @return theme
  #' @export
  theme_bw(base_size=12) %+replace% 
    theme(panel.grid.major=element_line(colour="gray", size=0.2),
          plot.title = element_text(hjust = 0.5, size=7),
          panel.border = element_rect(colour = "gray", fill=NA, size=1),
          legend.text=element_text(size=7),
          legend.title = element_blank(),
          legend.position="bottom",
          legend.justification = "top",
          legend.key.size = unit(0.25, "cm"),
          axis.text.x = element_text(size=7, angle=90, hjust=1),
          axis.text.y = element_text(size=8)
    )
}




bprle <- function (mat, title="RLE", ...){
  #' Prettier relative log expression (RLE) boxplots
  #' 
  #' @param mat: make it matrix.
  #' @param title (string): title of the box plot.
  #' @return boxplot
  #' @export
  #par(mar=c(20,2,1,1)) 
  x=mat
  #x=counts(set)
  y <- log2(x + 1)
  median <- apply(y, 1, median)
  rle <- apply(y, 2, function(x) x - median)
  #colors=brewer.pal(length(levels(pData(ruvs_set)[["batch"]])), "Set3")
  #boxplot(rle, ..., outline=FALSE, las=2, col=colors[pData(set)[["batch"]]], main=title, cex.axis=0.6, legend="bottom")
  boxplot(rle, ..., las=2, main=title)
  abline(h = 0, lty = 2)
}




pcaplt<- function(mat, title="PCA Plot"){
  #' PCA plot function
  #' @param mat: dataframe, datatable. First column is gene names, rows are genes, columns are samples.
  #' @param title (string): Title of PCA plot.
  #' @return pcaplot
  #' @export
  lg=log(mat[,-1]+1)
  #lg=mat
  var=lg[apply(lg, 1, var, na.rm=TRUE) != 0,]
  cc.var=var[complete.cases(var),]
  
  pca_prcomp = prcomp(t(var), center = T, scale = F)
  
  PC1_and_PC2 = data.frame(PC1=pca_prcomp$x[,1], PC2= pca_prcomp$x[,2], type = rownames(pca_prcomp$x))
  #PC1_and_PC2md =merge(PC1_and_PC2, md, by.x="type", by.y="target_id", all.x)
  
  perc=(pca_prcomp$sdev^2)/sum(pca_prcomp$sdev^2) *100
  labs <- sapply(seq_along(perc), function(i) {paste("PC ", i, " (", round(perc[i], 2), "%)", sep = "")})
  
  
  p=ggplot(PC1_and_PC2,aes(PC1, PC2)) + 
    geom_point(size=4) +
    geom_text(aes(label=type), vjust=-1) +
    labs(title=title, x=labs[1], y=labs[2]) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5),
          legend.text=element_text(size=4),
          legend.position="right")
  return(p)
}



pcaplt_grid_smmart <- function(mat, metadat=md, shape="diagnosis", title="PCA Plots"){
  #' PCA plot of multiple principle components. 1vs2 2vs3 1vs3
  #' Filtered by SMMART samples only
  #' 
  #' @param mat: full matrix of expression
  #' @param md is the metadata, a row for each column is ideal.
  #' @param title (string): title of the scatter.
  #' @return pca plots
  #' @export
  lg=log(mat+1)
  var=lg[apply(lg, 1, var, na.rm=TRUE) != 0,]
  cc.var=var[complete.cases(var),]
  
  pca_prcomp = prcomp(t(var), center = T, scale = F)
  
  PCs = data.frame(PC1=pca_prcomp$x[,1], PC2= pca_prcomp$x[,2], PC3= pca_prcomp$x[,3],type = rownames(pca_prcomp$x))
  md$adjnames = make.names(md$fastqname)
  
  PCsmd =merge(PCs, md, by.x="type", by.y="adjnames", all.x=TRUE)
  perc=(pca_prcomp$sdev^2)/sum(pca_prcomp$sdev^2) *100
  labs <- sapply(seq_along(perc), function(i) {paste("PC ", i, " (", round(perc[i], 2), "%)", sep = "")})
  
  lenshape=length(unique(PCsmd[,which(colnames(PCsmd)==shape)]))
  
  p1=ggplot(PCsmd,aes_string("PC1", "PC2", col="batch", shape=shape)) + 
    geom_point(size=4) +
    geom_text(aes(label=ifelse(diagnosis=="uhr",as.character(patid),'')),hjust=0,vjust=0) +
    scale_colour_brewer(palette = "Set1") +
    labs(title="PC1 vs PC2", x=labs[1], y=labs[2]) +
    guides(colour=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "gray", fill=NA, size=1),
          legend.text=element_text(size=7),
          legend.position="bottom")
    #legend.text=element_text(size=7),
    #legend.position="bottom")
  p2=ggplot(PCsmd,aes_string("PC1", "PC3", col="batch", shape=shape))  + 
    geom_point(size=4) +
    geom_text(aes(label=ifelse(diagnosis=="uhr",as.character(patid),'')),hjust=0,vjust=0) +
    scale_colour_brewer(palette = "Set1") +
    labs(title="PC1 vs PC3", x=labs[1], y=labs[3]) +
    guides(colour=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "gray", fill=NA, size=1),
          legend.text=element_text(size=7),
          legend.position="bottom")
  #legend.text=element_text(size=7),
  #legend.position="right")
  p3=ggplot(PCsmd,aes_string("PC2", "PC3", col="batch", shape=shape)) + 
    geom_point(size=4) +
    geom_text(aes(label=ifelse(diagnosis=="uhr",as.character(patid),'')),hjust=0,vjust=0) +
    scale_colour_brewer(palette = "Set1") +
    labs(title="PC2 vs PC3", x=labs[2], y=labs[3]) +
    guides(colour=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "gray", fill=NA, size=1),
          legend.text=element_text(size=7),
          legend.position="bottom")
  #legend.text=element_text(size=7),
  #legend.position="right")
  p=grid.arrange(p1,p2,p3, nrow=1, top=title)
  return(p)
}

pcaplt_grid <- function(mat, md, title="PCA Plots"){
  #' PCA plot of multiple principle components. 1vs2 2vs3 1vs3
  #' 
  #' @param mat: full matrix of expression
  #' @param md is the metadata, a row for each column is ideal.
  #' @param title (string): title of the scatter.
  #' @return pca plots
  #' @export
  lg=mat
  #lg=log(mat+1)
  
  var=lg[apply(lg, 1, var, na.rm=TRUE) != 0,]
  cc.var=var[complete.cases(var),]
  
  pca_prcomp = prcomp(t(var), center = T, scale = F)
  
  PCs = data.frame(PC1=pca_prcomp$x[,1], PC2= pca_prcomp$x[,2], PC3= pca_prcomp$x[,3],type = rownames(pca_prcomp$x))
  PCs$batch = "BCCL"
  PCs$batch[grep("RNA170328LH", PCs$type)] = "BCCL_UHR"
  PCs$batch[grep("RNA170210CC|RNA171102CC|NS500642", PCs$type)] = "SMMART"
  PCs$batch[grep("UHR", PCs$type)] = "uhr"
  PCs$batch[PCs$type %in% make.names(md$fastqname[md$diagnosis=="Breast"])] = "SMMART Breast Cancer"
  #PCs$batch[PCs$type %in% make.names(md$fastqname[md$diagnosis=="Pancreatic"])] = "Pancreatic_Cancer"
  
  perc=(pca_prcomp$sdev^2)/sum(pca_prcomp$sdev^2) *100
  labs <- sapply(seq_along(perc), function(i) {paste("PC ", i, " (", round(perc[i], 2), "%)", sep = "")})
  
  
  p1=ggplot(PCs,aes(PC1, PC2, col=batch)) + 
    geom_point(size=4) +
    #geom_text(aes(label=ifelse(grepl("RNA180618RS", type),as.character(type),'')),hjust=1,vjust=0) +
    geom_text(aes(label=ifelse(batch=="uhr",as.character(batch),'')),hjust=0,vjust=0) +
    scale_colour_brewer(palette = "Set1") +
    labs(title="PC1 vs PC2", x=labs[1], y=labs[2]) +
    guides(colour=guide_legend(ncol=1)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5),
          legend.position="bottom", 
          panel.border = element_rect(colour = "gray", fill=NA, size=1))
  #legend.text=element_text(size=7),
  p2=ggplot(PCs,aes(PC1, PC3, col=batch)) + 
    geom_point(size=4) +
    #geom_text(aes(label=ifelse(grepl("RNA180618RS", type),as.character(type),'')),hjust=1,vjust=0) +
    geom_text(aes(label=ifelse(batch=="uhr",as.character(batch),'')),hjust=0,vjust=0) +
    scale_colour_brewer(palette = "Set1") +
    labs(title="PC1 vs PC3", x=labs[1], y=labs[3]) +
    guides(colour=guide_legend(ncol=1)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5),
          legend.position="bottom", panel.border = element_rect(colour = "gray", fill=NA, size=1))
  #legend.text=element_text(size=7),
  p3=ggplot(PCs,aes(PC2, PC3, col=batch)) + 
    geom_point(size=4) +
    #geom_text(aes(label=ifelse(grepl("RNA180618RS", type),as.character(type),'')),hjust=1,vjust=0) +
    geom_text(aes(label=ifelse(batch=="uhr",as.character(batch),'')),hjust=0,vjust=0) +
    scale_colour_brewer(palette = "Set1") +
    labs(title="PC2 vs PC3", x=labs[2], y=labs[3]) +
    guides(colour=guide_legend(ncol=1)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5),
          legend.position="bottom")
  #legend.text=element_text(size=7),
  p=grid.arrange(p1,p2,p3, nrow=1, top=title)
  return(p)
}

pcaplt_grid_smmart_tissue <- function(mat, metadat=md, shape="tissue", title="PCA Plots"){
  #' PCA plot of multiple PCs with tissue not batch as shape condition
  #' 
  #' @param mat: matrix
  #' @param md is the metadata, a row for each column is ideal
  #' @param title (string): title of the scatter
  #' @return pca plot
  #' @export
  lg=log(mat+1)
  var=lg[apply(lg, 1, var, na.rm=TRUE) != 0,]
  cc.var=var[complete.cases(var),]
  
  pca_prcomp = prcomp(t(var), center = T, scale = F)
  
  PCs = data.frame(PC1=pca_prcomp$x[,1], PC2= pca_prcomp$x[,2], PC3= pca_prcomp$x[,3],type = rownames(pca_prcomp$x))
  md$adjnames = make.names(md$fastqname)
  
  PCsmd =merge(PCs, md, by.x="type", by.y="adjnames", all.x=TRUE)
  perc=(pca_prcomp$sdev^2)/sum(pca_prcomp$sdev^2) *100
  labs <- sapply(seq_along(perc), function(i) {paste("PC ", i, " (", round(perc[i], 2), "%)", sep = "")})
  
  lenshape=length(unique(PCsmd[,which(colnames(PCsmd)==shape)]))
  
  p1=ggplot(PCsmd,aes_string("PC1", "PC2", col="batch", shape=shape)) + 
    geom_point(size=4) +
    scale_shape_manual(values=c(16:18,7:10,12:14)[1:lenshape]) +
    geom_text(aes(label=ifelse(diagnosis=="uhr",as.character(patid),'')),hjust=0,vjust=0) +
    scale_colour_brewer(palette = "Set1") +
    labs(title="PC1 vs PC2", x=labs[1], y=labs[2]) +
    guides(colour=guide_legend(ncol=1)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "gray", fill=NA, size=1),
          legend.text=element_text(size=7),
          legend.position="bottom")
  
  p2=ggplot(PCsmd,aes_string("PC1", "PC3", col="batch", shape=shape))  + 
    geom_point(size=4) +
    scale_shape_manual(values=c(16:18,7:10,12:14)[1:lenshape]) +
    geom_text(aes(label=ifelse(diagnosis=="uhr",as.character(patid),'')),hjust=0,vjust=0) +
    scale_colour_brewer(palette = "Set1") +
    labs(title="PC1 vs PC3", x=labs[1], y=labs[3]) +
    guides(colour=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "gray", fill=NA, size=1),
          legend.text=element_text(size=7),
          legend.position="bottom")
  p3=ggplot(PCsmd,aes_string("PC2", "PC3", col="batch", shape=shape)) + 
    geom_point(size=4) +
    scale_shape_manual(values=c(16:18,7:10,12:14)[1:lenshape]) +
    geom_text(aes(label=ifelse(diagnosis=="uhr",as.character(patid),'')),hjust=0,vjust=0) +
    scale_colour_brewer(palette = "Set1") +
    labs(title="PC2 vs PC3", x=labs[2], y=labs[3]) +
    guides(colour=guide_legend(ncol=1), shape=guide_legend(ncol=1)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "gray", fill=NA, size=1),
          legend.text=element_text(size=7),
          legend.position="bottom")
  
  p=grid.arrange(p1,p2,p3, nrow=1, top=title)
  return(p)
}


hmap <- function(mat, title="heatmap",...){
  #' Heatmap function
  #' 
  #' Wrapper around pheatmap
  #' @param mat (matrix): RXC matrix with rownames and colnames
  #' @param title (string): name of plot title
  #' @return heatmap plot
  #' @examples
  #' 
  #' Additional parameters
  #' 
  #' #color            = c("blue", "red"), 
  #' #color            = c("#FFFFFF", colorRampPalette(rev(brewer.pal(9, "RdBu")))(200)),
  #' #cluster_cols     = FALSE,
  #' #cluster_rows     = FALSE,
  #' #label_rows       = 
  #' @export
	
	mat = mat[!rowSums(mat, na.rm=TRUE)==0, ]
	p=pheatmap(
	mat               = mat,
	color             = c("#FFFFFF", colorRampPalette(brewer.pal(9, "OrRd"))(200)),
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  drop_levels       = TRUE,
  fontsize          = 9,
	display_numbers = T,number_format = "%.1f", number_color="black",
  main              = title,...
  )
	return(p)
}


