suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
options(future.globals.maxSize = 70 * 1024^3)

# Create argument parser
parser <- ArgumentParser(description="Run single-cell RNA Analysis")
parser$add_argument("sample", help="Sample of interest")
parser$add_argument("samplesheet", help="Path to sample sheet CSV")
parser$add_argument("config", help="Path to YAML config file")
# parser$add_argument("--pipeline", help="Comma-separated list of pipelines to run", 
#                     default="init,ambient,create,qc,filter,processing,doublet,visualization")
args <- parser$parse_args()
message("++++++++ Starting pipeline ++++++++++")

# Parse the YAML file
params <- yaml::read_yaml(args$config)

set_defaults <- function(params) {
  defaults <- list(
    output = "output",
    project.prefix = "scRNA",
    verbose = FALSE,
    normalization = "SCT",
    RDS.file = FALSE,
    # QC parameters
    min_features = 200,
    max.nFeature_RNA = 7500,
    min.nFeature_RNA = 100,
    max.nCount_RNA = 50000,
    min.nCount_RNA = 200,
    max.percent.mt = 20
  )
  # For each default parameter, check if it exists in params
  for (param_name in names(defaults)) {
    if (!param_name %in% names(params)) {
      params[[param_name]] <- defaults[[param_name]]
    }
  }
  
  return(params)
}

# Apply default parameters
message("++++++++ Setting paramaters ++++++++++")
params <- set_defaults(params)
print(str(params))
# Setup log tracking
init_log <- function() {
  log <- list()
  log$timestamp <- Sys.time()
  log$pipelines.run <- params$pipeline
  log$status <- "running"
  log$step <- "start"
  return(log)
}

log <- init_log()
for (param.name in names(params)) {
  log[[paste0("params:", param.name)]] <- as.character(params[[param.name]])
}
# Parse pipeline steps to run
pipeline.to.run <- unlist(strsplit(params$pipeline, split = ","))
if(pipeline.to.run == "full"){
  pipeline.to.run <- c("init", "ambient", "create", "qc", "filter", "preprocessing", "doublet", "visualization")
}
# Read sample sheet
read_samplesheet <- function() {
  samplesheet <- read.csv(args$samplesheet)
  if (!"sampleName" %in% colnames(samplesheet)) {
    stop("Sample sheet must contain a 'sampleName' column")
  }
  if (!"path" %in% colnames(samplesheet)) {
    stop("Sample sheet must contain a 'path' column")
  }

  return(samplesheet)
}

samplesheet <- read_samplesheet()

sample <- args$sample
sample.path <- samplesheet$path[samplesheet$sampleName == sample]
message(samplesheet)
message(sample, sample.path)
log$sample <- sample
params$sample <- sample

# Load helper functions
source("scripts/functions.R")

message("++++++++ Setting directories ++++++++++")

# Create output directories
setup_directories <- function(sample, params) {
  dirs <- list()
  dirs$sample.dir <- file.path(params$output, "samples", sample)
  dirs$data.dir <- file.path(dirs$sample.dir, "data")
  dirs$cellranger <- file.path(dirs$data.dir, "cellranger-outs")
  dirs$raw.dir <- file.path(dirs$cellranger, "raw")
  dirs$filtered.dir <- file.path(dirs$cellranger, "filtered")
  dirs$output <- file.path(dirs$sample.dir, "output")
  dirs$processed.dir <- file.path(dirs$output, "RDS-files")
  dirs$plot.dir <- file.path(dirs$output, "plots")
  
  for (dir in c(dirs$raw.dir, dirs$processed.dir, dirs$filtered.dir, dirs$plot.dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
  
  return(dirs)
}

dirs <- setup_directories(sample, params)

message("Processing sample: ", sample)

# Save metrics
save_metrics <- function(status, log, dirs) {
  print("save metrics")
  log$status <- status
  log$timestamp.end <- Sys.time()
  log$runtime <- difftime(log$timestamp.end, log$timestamp, units = "mins")
  log.char <- lapply(log, function(x) {
  if (length(x) == 1 && (is.atomic(x) || is.factor(x))) {
    return(as.character(x))
  } else {
    return(NA_character_)  # Replace non-scalar or non-atomic entries
  }
  })
  log.df <- as.data.frame(log.char)
  write.csv(log.df, file = file.path(dirs$output, "metrics.csv"), row.names = FALSE)
  message("metrics saved: ", file.path(dirs$output, "metrics.csv"))
  return(log)
}


# Function to run initialization step
run_init <- function(samplesheet, sample, sample.path, dirs, log) {
  message("Beginning scRNA processing")
  
  metadata <- if (!is.null(samplesheet$metadata)) samplesheet$metadata else NA
  CombineDirectories(sample.path, dirs$raw.dir, dirs$filtered.dir, sample, metadata)
  
  log$path <- dirs$raw.dir
  log$step <- "init"
  
  return(log)
}

# Function to process ambient RNA
run_ambient <- function(dirs, sample.path, params, log) {
  message("Processing ambient RNA")
  
  if (params$type.sequencing == "snRNA") {
    message("Estimating ambient RNA based on mitochondrial contamination")
    r <- filter_ambient_RNA_sn(sample.path, params, dirs$plot.dir)
  } else if (exists("filter_ambient_RNA")) {
    message("Estimating ambient RNA using default automatic estimation of the contamination")
    r <- filter_ambient_RNA(sample.path, params, dirs$plot.dir)
  } else {
    message("Ambient RNA filtering functions not found, reading counts directly")
    r <- list(Read10X(data.dir = dirs$filtered.dir), NA)
  }
  
  counts <- r[[1]]
  log$contamination <- r[[2]]
  log$step <- "ambient"
  
  return(list(counts = counts, log = log))
}

# Function to create Seurat objects
run_create <- function(counts, samplesheet, sample, params, dirs, log) {
  message("Creating Seurat object")
  
  r <- safe_create_seurat_linear(counts, project = params$project.prefix)
  seu <- r[[1]]
  log$min.features <- r[[2]]
  print(seu)
  seu$orig.ident <- sample
  seu <- RenameCells(seu, add.cell.id = sample)
  
  # Adding metadata columns to Seurat object
  path.col <- grep("path", colnames(samplesheet))
  if (path.col < ncol(samplesheet)) {
    for (md.col.to.add in colnames(samplesheet)[(path.col+1):ncol(samplesheet)]) {
      seu[[md.col.to.add]] <- samplesheet[samplesheet$sampleName == sample, md.col.to.add]
    }
  }
  # Find metadata file (starts with "metadata.")
  metadata.file <- list.files(dirs$filtered.dir, pattern = "^metadata\\.(csv|tsv|txt)$", full.names = TRUE)

  # Read it based on extension
  if (length(metadata.file) > 0) {
    ext <- tools::file_ext(metadata.file)

    m <- switch(
      ext,
      csv = read.csv(metadata.file, stringsAsFactors = FALSE),
      tsv = read.delim(metadata.file, stringsAsFactors = FALSE),
      txt = read.delim(metadata.file, stringsAsFactors = FALSE),
      stop("Unsupported file type: ", ext)
    )
    seu@meta.data <- m
  }
  
  # Assay-specific metrics
  # seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  # ignore case since mouse are mt
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "(?i)^mt-")
  # seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "(?i)^rp[sl]")

  
  log$cell.count.initial <- ncol(seu)
  log$gene.count.initial <- nrow(seu)
  log$step <- "create"
  print(colnames(seu@meta.data))
  saveRDS(seu, file = file.path(dirs$processed.dir, 
                               paste0(params$project.prefix, "_", sample, "_create.rds")))
  
  return(list(seu = seu, log = log))
}

# Function to generate QC plots
run_qc <- function(seu, params, sample, dirs, log, fname = "") {
  message("Generating QC plots")
  
  pdf(file = file.path(dirs$plot.dir, paste0(params$project.prefix, "_", sample, fname, "_qc_plots.pdf")),
      height = 8, width = 12)

    list.of.vars <- list(
    "1" = "nCount_RNA",
    "2" = "nFeature_RNA",
    "3" = "percent.mt",
    "4" = "percent.rb"
    )
    p.list <-
    lapply(list.of.vars, function(vars.to.plot) {
        p <-
            seu@meta.data %>%
            dplyr::select(orig.ident, all_of(vars.to.plot)) %>%
            reshape2::melt() %>%
            ggplot(aes(x = orig.ident, y = value, fill = orig.ident)) +
            geom_violin() +
            facet_grid(~variable) +
            theme_classic() +
            theme(axis.text = element_text(color = "black"))

        if (!("percent.mt" %in% vars.to.plot | "percent.rb" %in% vars.to.plot)) {
            p <- p + scale_y_log10()
        }

        return(p)
    })

    print(p.list)

    density.scatter <- ggplot(seu@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
        geom_point(alpha = 0.5) +
        scale_x_log10() +  # if you want log scale
        theme_minimal()

    print(density.scatter)

    dev.off()
  
  log$step <- "qc"
  
  return(list(seu = seu, log = log))
}

# Function to filter outliers
run_filter <- function(seu, dirs, params, sample, log) {
  message("Filtering outliers")
  
  log$cell.count.prefilter <- ncol(seu)
  log$median.features.prefilter <- median(seu$nFeature_RNA)
  log$median.counts.prefilter <- median(seu$nCount_RNA)
  log$median.mt.percent.prefilter <- median(seu$percent.mt)
  seu <- filter_outliers(seu, "nFeature_RNA", min.nFeature_RNA = params$min.nFeature_RNA, max.nFeature_RNA = params$max.nFeature_RNA, log.transform = TRUE)
  seu <- filter_outliers(seu, "nCount_RNA", min.nCount_RNA = params$min.nCount_RNA, max.nCount_RNA = params$max.nCount_RNA, log.transform = TRUE)
  # seu <- filter_outliers(seu, "percent.rb", max.percent.mt = params$max.percent.mt, log.transform = FALSE)
  
  if (params$type.sequencing == "snRNA") {
    seu <- filter_outliers_mt_snRNA(seu, "percent.mt", max.percent.mt = params$max.percent.mt)
  } else {
    seu <- filter_outliers(seu, "percent.mt", max.percent.mt = params$max.percent.mt)
  } 
  
  log$cell.count.postfilter <- ncol(seu)
  log$percent.cells.kept.postfilter <- round(log$cell.count.postfilter / log$cell.count.prefilter * 100, 2)
  log$step <- "filter"
  log$median.features.postfilter <- median(seu$nFeature_RNA)
  log$median.counts.postfilter <- median(seu$nCount_RNA)
  log$median.mt.percent.postfilter <- median(seu$percent.mt)
  saveRDS(seu, file = file.path(dirs$processed.dir, 
                               paste0(params$project.prefix, "_", sample, "_filtered.rds")))
  
  return(list(seu = seu, log = log))
}

# Function to process data
run_preprocessing <- function(seu, params, dirs, sample, log) {
  message("Processing data")
  
  if (params$normalization == "SCT") {
    seu <- SCTransform(seu, verbose = F)
  } else {
    seu <- NormalizeData(seu) %>% 
      FindVariableFeatures() %>% 
      ScaleData()
  }
  
  seu <- RunPCA(seu, verbose = F) %>%
    FindNeighbors(dims = 1:30, verbose = F) %>%
    FindClusters(resolution = 0.5, verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F)
  
  log$clusters.identified <- length(unique(seu$seurat_clusters))
  log$step <- "processing"


  pdf(file.path(dirs$plot.dir, paste0(params$project.prefix, "_", sample, "_qc-clusters.pdf")))
    list.of.vars <- list(
      "1" = "nCount_RNA",
      "2" = "nFeature_RNA",
      "3" = "percent.mt"
    ) 
  
    p <- VlnPlot(seu, 
        pt.size = 0, 
        features = unlist(list.of.vars), 
        group.by = "seurat_clusters",
        stack = T, 
        fill.by = "ident", 
        log = T) + 
        NoLegend() +
        ggtitle(sample)

    print(p)
  dev.off()
  saveRDS(seu, file = file.path(dirs$processed.dir, 
                               paste0(params$project.prefix, "_", sample, "_preprocessed.rds")))
  
  return(list(seu = seu, log = log))
}

# Function to detect and filter doublets
run_doublet <- function(seu, params, dirs, sample, log) {
  message("Detecting doublets")
  
  seu <- filter_doublet(seu)
  
  # Re-process after doublet removal
  if (params$normalization == "SCT") {
    seu <- SCTransform(seu, verbose = F)
  } else {
    seu <- NormalizeData(seu) %>% 
      FindVariableFeatures() %>% 
      ScaleData()
  }
  
  seu <- RunPCA(seu, verbose = F) %>%
    FindNeighbors(dims = 1:30, verbose = F) %>%
    FindClusters(resolution = 0.5, verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F)
  
  log$cell.count.postdoublet <- ncol(seu)
  log$percent.doublets <- round((log$cell.count.postfilter - log$cell.count.postdoublet) / 
                                 log$cell.count.postfilter * 100, 2)
  log$step <- "doublet"
  
  saveRDS(seu, file = file.path(dirs$processed.dir, 
                               paste0(params$project.prefix, "_", sample, "_final.rds")))
  
  DefaultAssay(seu) <- "RNA"
  # Save a minimal version for easier sharing/loading
  seu.minimal <- DietSeurat(seu, 
                           counts = TRUE, 
                           data = TRUE, 
                           scale.data = FALSE,
                           assays = "RNA",
                           dimreducs = c("pca", "umap"))
  
  saveRDS(seu.minimal, file = file.path(dirs$processed.dir, 
                                      paste0(params$project.prefix, "_", sample, "_minimal.rds")))
  
  return(list(seu = seu, log = log))
}

# Function to create visualization plots
run_visualization <- function(seu, params, dirs, sample, log) {
  message("Creating visualization plots")
  
  # Common marker genes for major cell types
  marker.genes <- c(
    # Epithelial markers
    "EPCAM", "KRT19", "KRT8", "KRT18",
    # Immune markers
    "PTPRC", "CD3E", "CD4", "CD8A", "MARCO", "CD14", "CD68",
    # Stromal markers
    "MYLK", "DCN", "COL1A1", "VIM",
    # Endothelial markers
    "PECAM1", "VWF", "CDH5"
  )
  
  # Check which markers are present in the dataset
  valid.markers <- marker.genes[marker.genes %in% rownames(seu)]
  log$markers.found <- length(valid.markers)
  
  if (length(valid.markers) > 0) {
    DefaultAssay(seu) <- ifelse(params$normalization == "SCT", "SCT", "RNA")
    print(seu)
    pdf(file.path(dirs$plot.dir, paste0(params$project.prefix, "_", sample, "_markers.pdf")), , width = 14, height = 10)
  
      # UMAP colored by clusters
      p1 <- DimPlot(seu, reduction = "umap", label = TRUE) + 
        ggtitle("Clusters")
      print(p1)
      
      p <- FeaturePlot(seu, features = marker.genes, ncol = min(4, length(marker.genes)))
        print(p)
    
    dev.off()
  } else {
    message("No valid markers found in dataset — skipping FeaturePlot")
  }
  

  log$step <- "visualization"
  
  return(list(seu = seu, log = log))
}

# Main execution function
main <- function() {
  seu <- NULL
  counts <- NULL
  if ("full" %in% pipeline.to.run || "create" %in% pipeline.to.run) {
    if (!isFALSE(params$RDS.file)) {
      stop("ERROR: You gave an RDS to read but also want the pipeline to create the object. Please choose one.")
    }
  }
  print(paste0("log: ",str(log)))

  tryCatch({
    message("\n\n++++++++ Pipelines steps ++++++++++")
    message(pipeline.to.run)

    if (params$RDS.file != FALSE){
      message("\n\n++++++++ reading RDS file ++++++++++")
      seu <- readRDS(params$RDS.file)
      log$step <- "read RDS"
    }
    # Run pipeline steps based on user selection
    if ("init" %in% pipeline.to.run) {
      message("\n\n++++++++ init ++++++++++")
      log <- run_init(samplesheet, sample, sample.path, dirs, log)
    }
    
    if ("ambient" %in% pipeline.to.run) {
      message("\n\n++++++++ ambient ++++++++++")
      result <- run_ambient(dirs, sample.path, params, log)
      counts <- result$counts
      log <- result$log
    } else {
      counts <- Read10X(data.dir = dirs$raw.dir)
      log$contamination <- NA
    }
      print(paste0("log: ",str(log)))

    if ("create" %in% pipeline.to.run) {
      message("\n\n++++++++ create ++++++++++")
      result <- run_create(counts, samplesheet, sample, params, dirs, log)
      seu <- result$seu
      log <- result$log
    }
    
    if ("qc" %in% pipeline.to.run && !is.null(seu)) {
      message("\n\n++++++++ qc ++++++++++")
      result <- run_qc(seu, params, sample, dirs, log)
      seu <- result$seu
      log <- result$log
    }
    
    if ("filter" %in% pipeline.to.run && !is.null(seu)) {
      message("\n\n++++++++ filter ++++++++++")
      message("Number of cells prefiltering: ", ncol(seu))
      result <- run_filter(seu, dirs, params, sample, log)
      seu <- result$seu
      log <- result$log
      message("Number of cells postfiltering: ", ncol(seu))
      result <- run_qc(seu, params, sample, dirs, log, "_postfiltering")
    }
      print(paste0("log: ",str(log)))

    if ("preprocessing" %in% pipeline.to.run && !is.null(seu)) {
      message("\n\n++++++++ preprocessing ++++++++++")
      result <- run_preprocessing(seu, params, dirs, sample, log)
      seu <- result$seu
      log <- result$log
    }
    
    if ("doublet" %in% pipeline.to.run && !is.null(seu)) {
      message("\n\n++++++++ doublet ++++++++++")
      message("Number of cells pre doublet filtering: ", ncol(seu))
      result <- run_doublet(seu, params, dirs, sample, log)
      seu <- result$seu
      log <- result$log
      message("Number of cells post doublet filtering: ", ncol(seu))
      log$median.features <- median(seu$nFeature_RNA)
      log$median.counts <- median(seu$nCount_RNA)
    }
      print(paste0("log: ",str(log)))

    if ("visualization" %in% pipeline.to.run && !is.null(seu)) {
      message("\n\n++++++++ visualization ++++++++++")
      result <- run_visualization(seu, params, dirs, sample, log)
      seu <- result$seu
      log <- result$log
    }
      print(paste0("log: ",str(log)))

    # Successfully completed all requested pipelines
    log <- save_metrics("completed", log, dirs)
    message("\n\nProcessing completed successfully for sample: ", sample)
    
  }, error = function(e) {
    # Handle errors
    log$error <- as.character(e)
    message("Error in processing: ", e$message)
    log$error <- e
    log <- save_metrics("failed", log, dirs)
    message("Processing failed for sample: ", sample)
  }, finally = {
    # This code will run regardless of success or failure
    if (is.null(log$status) || (log$status != "completed" & log$status != "failed")) { 
      message("Ensuring metrics are saved even if processing was interrupted")
      log$status <- "interrupted"
      log <- save_metrics("interrupted", log, dirs)
    }
  })
}

# Execute main function
main()
