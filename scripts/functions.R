# Author: Romane Cathelin
# Date: 2025-05-07


CombineDirectories <- function(data.dir, new.dir, new.dir.filter, sample, metadata) {
  # Create the output directory if it doesn't exist
  if (!dir.exists(new.dir)) {
    dir.create(new.dir, recursive = TRUE)
  }

  # Symlink GEX files
  gex.dir <- file.path(data.dir, "filtered_feature_bc_matrix")
  gex.files <- list.files(path = gex.dir, full.names = TRUE)
  lapply(gex.files, function(file.to.link) {
    file.symlink(from = file.to.link,
                 to   = file.path(new.dir.filter, basename(file.to.link)))
  })

  gex.dir.raw <- file.path(data.dir, "raw_feature_bc_matrix")
  gex.files.raw <- list.files(path = gex.dir.raw, full.names = TRUE)
  lapply(gex.files.raw, function(file.to.link) {
    file.symlink(from = file.to.link,
                 to   = file.path(new.dir, basename(file.to.link)))
  })

  # Symlink metadata file if it exists and is not NA, and rename it
  if (!is.na(metadata) && file.exists(metadata)) {
    # Extract the file extension
    extension <- tools::file_ext(metadata)
    new.metadata.name <- paste0("metadata.", extension)

    file.symlink(from = metadata,
                 to   = file.path(new.dir.filter, new.metadata.name))
  }
}


# Function to save metrics to CSV
save.metrics <- function(log, sample, sample.dir, status = "completed") {
  log$sample <- sample
  log$completed.time <- Sys.time()
  log$runtime.minutes <- round(difftime(log$completed.time, log$timestamp, units = "mins"), 2)
  log$status <- status
  
  # Convert log to data frame for CSV output
  log.df <- data.frame(
    sample = sample,
    timestamp = as.character(log$timestamp),
    completed.time = as.character(log$completed.time),
    runtime.minutes = log$runtime.minutes,
    last.step = log$step,
    status = log$status,
    initial.cells = if(exists("cell.count.initial", where = log)) log$cell.count.initial else NA,
    cells.after.filtering = if(exists("cell.count.postfilter", where = log)) log$cell.count.postfilter else NA,
    cells.after.doublet = if(exists("cell.count.postdoublet", where = log)) log$cell.count.postdoublet else NA,
    percent.cells.kept = if(exists("percent.cells.kept", where = log)) log$percent.cells.kept else NA,
    percent.doublets = if(exists("percent.doublets", where = log)) log$percent.doublets else NA,
    clusters.identified = if(exists("clusters.identified", where = log)) log$clusters.identified else NA,
    pipelines.run = log$pipelines.run
  )
  
  # Write metrics to CSV
  write.csv(log.df, file = file.path(sample.dir, "metrics.csv"), row.names = FALSE)
  message("Last completed step: ", log$step)
  message("Metrics saved to: ", file.path(sample.dir, "metrics.csv"))
  
  return(log.df)
}

############################################################################
################                                            ################
################    FUNCTIONS FOR QC SCRNA-SEQ/ SNRNA-SEQ   ################
################                                            ################
############################################################################


suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(DropletUtils))

#' Filter out ambient RNA using default protocol (SoupX)
#'
#' @param matrix.path String. Path to the CellRanger output directory. Must contain 
#' `filtered_feature_bc_matrix` and `raw_feature_bc_matrix`.
#'
#' @return A matrix of counts with ambient RNA filtered out.
#' @export
filter_ambient_RNA <- function(data.path, params, plot.path, VERBOSE = T){
  # tod = Seurat::Read10X(raw.path)
  # toc = Seurat::Read10X(filtered.path)
  # sc = SoupChannel(tod,toc)
  sc <- load10X(dataDir = data.path, keepDroplets = TRUE)
  tryCatch({
    pdf(file.path(plot.path,  paste0(params$project.prefix, "_", params$sample, "_qc-contamination.pdf")))
      sc <- autoEstCont(sc)
    dev.off()
    if (VERBOSE) return(list(count = adjustCounts(sc), message = "Decontamination worked")) else return(adjustCounts(sc))
    
  }, error = function(e) {
    message("autoEstCont failed, not filtering ambient RNA")
    if (VERBOSE) return(list(count = sc$toc, message = "Decontamination failed")) else return(sc$toc)
  })

}

#' Filter ambient RNA for single nuclei sequencing data
#'
#' Filter out ambient RNA for single nuclei data using mitochondrial genes as a proxy for contamination.
#'
#' @param matrix.path String. Path to the CellRanger output directory. Must contain 
#' `filtered_feature_bc_matrix` and `raw_feature_bc_matrix`.
#'
#' @return A matrix of counts with ambient RNA filtered out based on mitochondrial contamination.
#' @export
filter_ambient_RNA_sn <- function(data.path, params, plot.path, VERBOSE = T){
  # tod = Seurat::Read10X(raw.path)
  # toc = Seurat::Read10X(filtered.path)
  # sc = SoupChannel(tod,toc)
  print(data.path)
  
  sc <- load10X(dataDir = data.path, keepDroplets = TRUE)
  genes <- rownames(sc$toc)
  gene.list <- list(MT = genes[grepl("^MT-", genes)])
  non.expressing <- estimateNonExpressingCells(sc, gene.list)
  sc <- calculateContaminationFraction(sc, gene.list, useToEst = non.expressing)
  print("counts")
  pdf(file.path(plot.path,  paste0(params$project.prefix, "_", params$sample, "_qc-contamination.pdf")))
  counts <- adjustCounts(sc)
  dev.off()
  if (VERBOSE) return(list(count = counts, message = "Decontamination mt worked")) else return(counts)
}



safe_create_seurat_linear <- function(counts, min.cells = 3, project = "LUAD", VERBOSE = T) {
  feature_thresholds <- seq(200, 1, by = -10)
  for (min.features in feature_thresholds) {
    try_result <- try(CreateSeuratObject(counts = counts, project = project, min.cells = min.cells, min.features = min.features), silent = TRUE)
    if (!inherits(try_result, "try-error")) {
      message("✅ Using min.features = ", min.features)
      if (VERBOSE) return(list(data = try_result, message = paste0("min.features = ", min.features))) else return(try_result)
      return(try_result)
    }
  }
  stop("❌ Could not create Seurat object even at min.features = 1")
}

#' Filter outliers using MAD-based threshold
#'
#' Detect and filter out outliers based on the median absolute deviation (MAD) from the median value.
#'
#' @param data SeuratObject. Original Seurat object.
#' @param metric String. Name of the metadata column to use as the metric.
#' @param log.transform Logical. Whether to log-transform the data before computing outliers. Default is FALSE.
#' @param nmads Numeric. Number of MADs away from median to consider a value an outlier. Default is 5.
#'
#' @return A filtered SeuratObject with outliers removed.
#' @export
filter_outliers <- function(data, metric, max.nFeature_RNA = 7500, min.nFeature_RNA = 100, max.nCount_RNA = 50000, min.nCount_RNA = 200, max.percent.mt = 20, log.transform = FALSE, nmads = 3) {
  outliers <- isOutlier(data@meta.data[[metric]], nmads = nmads, log = log.transform)
  outlier.col <- paste0("outlier_", metric)
  data[[outlier.col]] <- outliers
  expr <- FetchData(object = data, vars = outlier.col)
  data <- data[, which(expr == FALSE)]

  if (metric == "percent.mt") {
    data <- subset(data, subset = percent.mt < max.percent.mt)
  } else if (metric == "nFeature_RNA") {
    data <- subset(data, subset = nFeature_RNA > min.nFeature_RNA & nFeature_RNA < max.nFeature_RNA)
  } else if (metric == "nCount_RNA") {
    data <- subset(data, subset = nCount_RNA > min.nCount_RNA & nCount_RNA < max.nCount_RNA)
  }

  return(data)
}

#' Filter mitochondrial outliers in single nuclei data
#'
#' Filters out cells with mitochondrial gene expression in single nuclei RNA-seq data.
#' Assumes no mitochondrial RNA should be present in high-quality single nuclei.
#'
#' @param data SeuratObject. Original Seurat object.
#' @param metric String. Name of the metric column (typically percent.mt).
#'
#' @return A filtered SeuratObject with mitochondrial outliers removed.
#' @export
filter_outliers_mt_snRNA <- function(data, metric, nmads = 3, max.percent.mt = 5){
  data$outlier_percent.mt <- isOutlier(data$percent.mt, type = "higher", min.diff = 0.5, nmads = nmads)
  expr <- FetchData(object = data, vars = "outlier_percent.mt")
  data <- data[, which(x = expr == F)]

  #Add a threshold if the distribution is really high, isOutlier might keep to much 
  data <- subset(data, subset = percent.mt < max.percent.mt)

  return(data)
}

#' Estimate 10X Genomics doublet formation rate
#'
#' Predict the expected doublet rate based on the number of cells loaded into a 10X run.
#' Reference: \url{https://kb.10xgenomics.com/hc/en-us/articles/360001378811}
#'
#' @param x Integer. Number of cells loaded.
#'
#' @return Numeric. Predicted doublet formation rate (percentage).
#' @export
get_10x_multiplets <- function(x){
  expected.df <- data.frame(
    "Cells.Recovered" = c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000),
    "Multiplet.Rate" = c(0.4, 0.8, 1.6, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4, 7.2, 8.0)
  )
  lm.res <- lm(Multiplet.Rate ~ Cells.Recovered, data = expected.df)
  test.df <- data.frame("Cells.Recovered" = x)
  return(predict(lm.res, test.df))
}


filter_doublet <- function(data, VERBOSE = T){
  # removing doublet

  
  suppressMessages({ invisible(capture.output({
  data.sweep <- paramSweep(data, PCs = 1:20, sct = T)
  data.sweep.stats <- summarizeSweep(data.sweep, GT = FALSE)
  data.bcmvn <- find.pK(data.sweep.stats)
  data.pK <- as.numeric(as.character(data.bcmvn$pK[which.max(find.pK(data.sweep.stats)$BCmetric)]))
  homotypic.prop <- modelHomotypic(data$seurat_clusters)  
  mlt.est <- get_10x_multiplets(nrow(data@meta.data))/100
  nExp_poi <- round(mlt.est*nrow(data@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  data <- doubletFinder(data, PCs = 1:20, pN = 0.25,  pK = data.pK,  nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)
  data <- doubletFinder(data, PCs = 1:20, pN = 0.25, pK = data.pK, nExp = nExp_poi.adj, sct = TRUE, 
                     reuse.pANN = grep("pANN", colnames(data@meta.data), value = TRUE))
  })) })
  doublet.col <- grep("^DF.classifications", colnames(data@meta.data), value = TRUE)
  
  expr <- FetchData(object = data, vars = doublet.col)
  data <- data[, which(x = expr == "Singlet")]
  return(data)
}



########
subset_process <-function(data, res){
  data <- NormalizeData(data) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
  data <- RunHarmony(data, group.by.vars = "ExternalSubjectID")
  data <- RunUMAP(data, reduction = "harmony", dims = 1:30)
  print(res)
  data <- FindNeighbors(data, reduction = "harmony", dims = 1:30) 
  print("clusters")
  data <- FindClusters(data,  resolution = c(res))
  return(data)
}


find_and_transfer_labels <- function(ref.data, 
                                    query.data, 
                                    output.dir,
                                    dims = 1:20,
                                    integration = FALSE,
                                    save = F,
                                    normalization = "LogNormalize", 
                                    reference.reduction = "pca",
                                    ref.metadata.col = "ann_level_1",
                                    save.anchors = F,
                                    prefix = "") {
  
  print("Loading data..")
  # ref.data <- readRDS(reference.path)


  print("FindTransferAnchors..")
  query.anchors <- FindTransferAnchors(
    reference = ref.data, 
    query = query.data, 
    dims = dims, 
    reference.reduction = reference.reduction,
    normalization.method = normalization
  )
  
  if (save.anchors) {
    anchors.path <- file.path(output.dir, paste0(prefix, "query_anchors.rds"))
    saveRDS(query.anchors, file = anchors.path)
  }
  
  print("TransferData..")
  predictions <- TransferData(
    anchorset = query.anchors, 
    refdata = ref.data@meta.data[[ref.metadata.col]]
  )
  
  query.data <- AddMetaData(query.data, metadata = predictions)
  
  # pdf(file.path(output.dir, "UMAP_annotation.pdf"))
  #     d <- DimPlot(query.data, group.by = "predicted.id")
  #     print(d)
  #     if(integration){
  #         d1 <- DimPlot(query.data, group.by = "predicted.id", reduction = "UMAP_harmony")
  #         print(d1)
  #     }
  # dev.off()
  if(save){
    print("Saving data..")
    result.path <- file.path(output.dir, paste0(prefix, "processed_data_TransferData_pred.rds"))
    saveRDS(query.data, result.path)
  }
  return(query.data)
}
