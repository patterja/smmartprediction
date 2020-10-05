#' DrugPrediction TRAINING
#'
#' @param path_to_bcclsmmart (string):path to FINAL_norm.txt run through normalizeRUV in SMMARTFunctions.R
#' @param cell_line_response (string):cell_line_response_threshold_0.50_large_and_small_screen.RData
#' @param path_to_targetid (string): path to target_id.txt
#' @param path_to_trusight (string): path to gene list Trusight.csv
#' @param threshold (numeric): 0.25, 0.5, 0.75
#' @return smmart_trained_machine_learning_model.RData
#' @export
#'

drugPredTrain <- function(path_to_bcclsmmart,
                          cell_line_response = "cell_line_response_threshold_0.50_large_and_small_screen.RData",
                          path_to_targetid = "target_id.txt",
                          path_to_trusight = "Trusight_genes.csv",
                          threshold = 0.5) {

  normalization_type <- "none"
  classification_type <- "trusight"

  threshold <- threshold

  target_id <- read.csv(path_to_targetid, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  target_id <- target_id[which(target_id$X.7 == "protein_coding"),]
  target_id <- unique(target_id[, c("X.1", "X.5")])
  rownames(target_id) <- target_id$X.1

  load(cell_line_response)
  JWGray_GR50 <- cell_line_response

  JWGray_RNASeq <- read.csv(path_to_bcclsmmart, row.names = 1, check.names = FALSE, sep = "\t")
  JWGray_RNASeq <- JWGray_RNASeq[,1:70]
  colnames(JWGray_RNASeq) <- toupper(sapply(colnames(JWGray_RNASeq), function(name) {strsplit(name, split = "_", fixed = TRUE)[[1]][1]}))
  colnames(JWGray_RNASeq)[1:7] <- substr(colnames(JWGray_RNASeq)[1:7], 2, nchar(colnames(JWGray_RNASeq)[1:7]))
  JWGray_RNASeq <- t(JWGray_RNASeq)
  JWGray_RNASeq <- matrix(as.numeric(JWGray_RNASeq), nrow = nrow(JWGray_RNASeq), ncol = ncol(JWGray_RNASeq), dimnames = list(rownames(JWGray_RNASeq), colnames(JWGray_RNASeq)))

  common_cell_lines <- intersect(rownames(JWGray_RNASeq), rownames(JWGray_GR50))
  JWGray_RNASeq <- JWGray_RNASeq[common_cell_lines,]
  JWGray_GR50 <- JWGray_GR50[common_cell_lines,]

  hits <- which(rownames(target_id) %in% colnames(JWGray_RNASeq))
  JWGray_RNASeq <- JWGray_RNASeq[,rownames(target_id)[hits]]
  colnames(JWGray_RNASeq) <- target_id$X.5[hits]

  JWGray_RNASeq <- log2(JWGray_RNASeq + 1)

  X_train <- scale(JWGray_RNASeq)
  Y_train <- JWGray_GR50
  frequencies <- colSums(Y_train == 1, na.rm = TRUE) / colSums(is.na(Y_train) == FALSE)
  Y_train <- Y_train[,which(frequencies > 0.05 & frequencies < 0.95)]

  valid_genes <- which(colSums(is.na(X_train)) == 0)
  X_train <- X_train[, valid_genes]

  if (classification_type == "trusight") {
    trusight_genes <- read.csv(path_to_trusight, header = FALSE, stringsAsFactors = FALSE)[,1]
    X_train <- X_train[, colnames(X_train) %in% trusight_genes]
  }

  K_train <- X_train %*% t(X_train)
  normalizer <- max(abs(K_train))
  K_train <- K_train / normalizer

  parameters <- list()
  parameters$alpha_lambda <- 1
  parameters$beta_lambda <- 1
  parameters$iteration <- 200
  parameters$margin <- 1
  parameters$R <- 20
  parameters$seed <- 1606
  parameters$sigma_h <- 0.1
  parameters$sigma_w <- 1.0

  state <- kbmtl_semisupervised_classification_variational_train(K_train, Y_train, parameters)
  save(state, file =  "smmart_trained_machine_learning_model.RData")
  print(paste0("kbmtl trained model saved ", getwd(), "/smmart_trained_machine_learning_model.RData"))
  return(state)
}

