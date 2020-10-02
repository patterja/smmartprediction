#' DrugPrediction Test
#'
#' version 2
#' @param path_to_bcclsmmart (string):path to FINAL_norm.txt run through normalizeRUV in SMMARTFunctions.R
#' @param state (string): path to "smmart_trained_machine_learning_model.RData" from drugPredTrain
#' @param path_to_targetid (string): path to target_id.txt
#' @param path_to_trusight (string): path to gene list Trusight.csv
#' @param kbmtl_train (string): "kbmtl_semisupervised_classification_variational_train.R"
#' @param kbmtl_test (string): "kbmtl_semisupervised_classification_variational_test.R"
#' @param threshold (numeric): 0.25, 0.5, 0.75
#' @return Y_predicted(data.frame): prediction
#'
#' @export

drugPredTestonTrain<- function(path_to_bcclsmmart,
                        cell_line_response = "cell_line_response_threshold_0.50_large_and_small_screen.RData",
                        state = "smmart_trained_machine_learning_model.RData",
                        kbmtl_train ="kbmtl_semisupervised_classification_variational_train.R",
                        kbmtl_test = "kbmtl_semisupervised_classification_variational_test.R",
                        path_to_targetid = "target_id.txt",
                        path_to_trusight = "Trusight_genes.csv",
                        threshold = 0.5)
{
  source(kbmtl_train)
  source(kbmtl_test)

  ###########TESTING
  normalization_type <- "none"
  classification_type <- "trusight"
  threshold <- threshold

  #Gene Names
  target_id <- read.csv(path_to_targetid, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  target_id <- target_id[which(target_id$X.7 == "protein_coding"),]
  target_id <- unique(target_id[, c("X.1", "X.5")])
  rownames(target_id) <- target_id$X.1

  #Drug Screen data, creates cell_line_response matrix
  load(cell_line_response)
  JWGray_GR50 <- cell_line_response

  #Training
  JWGray_RNASeq <- read.csv(path_to_bcclsmmart, row.names = 1, check.names = FALSE, sep = "\t")

  JWGray_RNASeq <- JWGray_RNASeq[,1:70]
  colnames(JWGray_RNASeq) <- toupper(sapply(colnames(JWGray_RNASeq), function(name) {strsplit(name, split = "_", fixed = TRUE)[[1]][1]}))
  colnames(JWGray_RNASeq)[1:7] <- substr(colnames(JWGray_RNASeq)[1:7], 2, nchar(colnames(JWGray_RNASeq)[1:7]))
  JWGray_RNASeq <- t(JWGray_RNASeq)
  JWGray_RNASeq <- matrix(as.numeric(JWGray_RNASeq), nrow = nrow(JWGray_RNASeq), ncol = ncol(JWGray_RNASeq), dimnames = list(rownames(JWGray_RNASeq), colnames(JWGray_RNASeq)))

  #Commmon Drug Screen and Cell lines
  common_cell_lines <- intersect(rownames(JWGray_RNASeq), rownames(JWGray_GR50))
  JWGray_RNASeq <- JWGray_RNASeq[common_cell_lines,]
  JWGray_GR50 <- JWGray_GR50[common_cell_lines,]

  #Test Data
  SMMART_RNASeq <- read.csv(path_to_bcclsmmart, row.names = 1, check.names = FALSE, sep = "\t")
  SMMART_RNASeq = SMMART_RNASeq[, grepl("NS500642|NS500681", colnames(SMMART_RNASeq))]
  SMMART_RNASeq = SMMART_RNASeq[, !grepl("RNA170328LH|RNA170329LH", colnames(SMMART_RNASeq))]
  SMMART_RNASeq <- t(SMMART_RNASeq)
  SMMART_RNASeq <- matrix(as.numeric(SMMART_RNASeq), nrow = nrow(SMMART_RNASeq), ncol = ncol(SMMART_RNASeq), dimnames = list(rownames(SMMART_RNASeq), colnames(SMMART_RNASeq)))
  SMMART_RNASeq <- SMMART_RNASeq[grepl(rownames(SMMART_RNASeq), pattern = "UHR|horizon") == FALSE,]

  hits <- which(rownames(target_id) %in% colnames(JWGray_RNASeq))
  JWGray_RNASeq <- JWGray_RNASeq[,rownames(target_id)[hits]]
  SMMART_RNASeq <- SMMART_RNASeq[,rownames(target_id)[hits]]
  colnames(JWGray_RNASeq) <- target_id$X.5[hits]
  colnames(SMMART_RNASeq) <- target_id$X.5[hits]

  JWGray_RNASeq <- log2(JWGray_RNASeq + 1)
  SMMART_RNASeq <- log2(SMMART_RNASeq + 1)

  X_train <- scale(JWGray_RNASeq)
  Y_train <- JWGray_GR50
  frequencies <- colSums(Y_train == 1, na.rm = TRUE) / colSums(is.na(Y_train) == FALSE)
  Y_train <- Y_train[,which(frequencies > 0.05 & frequencies < 0.95)]
  X_test <- X_train

  valid_genes <- which(colSums(is.na(X_train)) == 0)

  X_train <- X_train[, valid_genes]
  X_test <- X_test[, valid_genes]
  X_test[is.na(X_test)] <- 0

  if (classification_type == "trusight") {
    trusight_genes <- read.csv(path_to_trusight, header = FALSE, stringsAsFactors = FALSE)[,1]
    X_train <- X_train[, colnames(X_train) %in% trusight_genes]
    X_test <- X_test[, colnames(X_test) %in% trusight_genes]
  }

  K_train <- X_train %*% t(X_train)
  K_test <- X_train %*% t(X_test)
  normalizer <- max(abs(K_train))
  K_train <- K_train / normalizer
  K_test <- K_test / normalizer

  ###########TRAINING

  load(file=state)

  prediction <- kbmtl_semisupervised_classification_variational_test(K_test, state)

  Y_predicted <- as.data.frame(prediction$P)
  colnames(Y_predicted) <- colnames(Y_train)
  return(Y_predicted)

}

