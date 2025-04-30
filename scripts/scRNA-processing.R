suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
options(future.globals.maxSize = 70 * 1024^3)

# Create argument parser
parser <- ArgumentParser(description="Run single-cell RNA Analysis")
parser$add_argument("samplesheet", help="Path to sample sheet CSV")
parser$add_argument("config", help="Path to YAML config file")
args <- parser$parse_args()

# Parse the YAML file
params <- yaml::read_yaml(args$config)
message("++++++++ Starting pipeline ++++++++++")

set_defaults <- function(params) {
  defaults <- list(
    output_dir = "output",
    min_cells = 200,
    normalization = "SCT",
    var.to.regress = NULL,
    cc = FALSE,
    variable.features.n = 5000,
    n.dim = 30,
    n.dim.integration = 40,
    cluster_resolutions = c(0.1, 0.2, 0.3),
    integration = "CCA",
    ref.metadata.col = "cell_type"
  )
  
  # For each default parameter, check if it exists in params
  for (param_name in names(defaults)) {
    if (!param_name %in% names(params)) {
      params[[param_name]] <- defaults[[param_name]]
    }
  }
  return(params)
}

# Load helper functions
source("scripts/functions.R")

# Apply default parameters
params <- set_defaults(params)
# Print summary of parameters used
  cat("\nParameters used in this analysis:\n")
  for (param_name in sort(names(params))) {
    if (!is.list(params[[param_name]]) && !is.null(params[[param_name]])) {
      cat(sprintf("- %s: %s\n", param_name, paste(params[[param_name]], collapse=", ")))
    }
  }

# Parse pipeline steps to run
pipeline.to.run <- unlist(strsplit(params$pipeline, split = ","))
if("full" %in% pipeline.to.run){
  pipeline.to.run <- c("filter", "merge", "qc", "processing", "integration", "annotation")
}

# Function to run filtering step
run_filter <- function(in.path, dirs, params) {
  metrics.file <- read.csv(file.path(in.path, "qc_metrics.csv"))
  samples.id <- as.character(metrics.file[metrics.file$cell.count.postfilter > params$min_cells & 
                         metrics.file$status == "completed", "sample"])
  
  message(paste("Processing", length(samples.id), "samples"))
  sample.dirs <- file.path(in.path, samples.id)
  
  # Load and prepare Seurat objects
  seurat.list <- lapply(seq_along(sample.dirs), function(i) {
    dir <- sample.dirs[i]
    dir <- file.path(dir, "output/RDS-files")
    sample <- samples.id[i]
    # fn <- paste0(params$project.prefix, "_", sample, "_minimal.rds")
    fn <- paste0(params$project.prefix, "_", sample, "_minimal.rds")
    print(paste0("loading file: ", file.path(dir, fn)))
    obj <- readRDS(file.path(dir, fn))
    print(obj)
    DefaultAssay(obj) <- "RNA"
    obj$sample <- sample
    return(obj)
  })
  
  saveRDS(seurat.list, file = file.path(dirs$RDS.files, 
                               paste0(params$project.prefix, "-create-obj-list.RDS")))
  message("filtering done.")
  return(seurat.list)
}

# Function to merge objects
run_merge <- function(seurat.list, dirs) {
  message("Merging samples...")
  first.obj <- seurat.list[[1]]
  remaining.objs <- seurat.list[-1]
  seu <- merge(first.obj, y = remaining.objs, project = "LUAD")
  saveRDS(seu, file = file.path(dirs$RDS.files, 
                               paste0(params$project.prefix, "-merged-obj.RDS")))
  message("merging done.")
  return(seu)
}

run_qc <- function(seu, params, dirs){
  message("QC...")
  # Join layers if necessary
  joined.seu <- if ("counts" %in% Layers(seu)) JoinLayers(seu) else seu
  Idents(joined.seu) <- "sample"
  
  # Pre-QC visualization
  p1 <- VlnPlot(joined.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  # Pre-QC cell counts
  ncells <- joined.seu@meta.data %>%
    group_by(sample) %>%
    summarize(n_cells = n())
  message("Number of cells per samples pre QC: \n", ncells)
  
  # Filter outliers
  joined.seu <- filter_outliers(joined.seu, metric = "nFeature_RNA", nmads = 3, log.transform = TRUE)
  joined.seu <- filter_outliers(joined.seu, metric = "nCount_RNA", nmads = 3, log.transform = TRUE)
  joined.seu <- filter_outliers(joined.seu, metric = "percent.mt", nmads = 3, log.transform = FALSE)
  
  # Post-QC cell counts
  ncells <- joined.seu@meta.data %>%
    group_by(sample) %>%
    summarize(n_cells = n())
  message("Number of cells per samples post QC: \n", ncells)
  
  # Post-QC visualization
  Idents(joined.seu) <- "sample"
  p2 <- VlnPlot(joined.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  # Save QC plots
  pdf(file.path(dirs$plots, paste0(params$project.prefix, "-qc.pdf")))
  print(p1)
  print(p2)
  dev.off()
  
  # Save QC'd object
  saveRDS(joined.seu, file = file.path(dirs$RDS.files,
                                      paste0(params$project.prefix, "-postqc-obj.RDS")))
  
  # Here's the key fix: Use joined.seu, not data in the subset
  seu <- subset(seu, cells = colnames(joined.seu))
  
  message("qc done.")
  return(seu)
}

# Function to run data processing
run_processing <- function(seu, params, dirs) {
  message("Processing data...")
  # Add cell cycle scoring if needed
  add_cell_cycle <- function(obj, assay_type) {
    if(params$cc) {
      message("\t - Cell cycle scoring")
      s.genes <- cc.genes$s.genes
      g2m.genes <- cc.genes$g2m.genes
      obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, 
                             assay = assay_type, set.ident = TRUE, verbose = F)
    }
    
    return(obj)
  }
  
  # Process based on normalization method
  if(params$normalization == "SCT") {
    message("\t - SCT normalization")
    seu <- SCTransform(seu, vars.to.regress = params$var.to.regress, 
                      assay = 'RNA', new.assay.name = 'SCT', 
                      variable.features.n = params$variable.features.n, verbose = F)
    
    seu <- add_cell_cycle(seu, 'SCT')
    
    # Re-normalize with cell cycle if needed
    if(params$cc) {
      message("\t - SCT normalization - regressing out cell cycling genes")
      seu <- SCTransform(seu, variable.features.n = params$variable.features.n, 
                        assay = 'RNA', new.assay.name = 'SCT', 
                        vars.to.regress = c(params$var.to.regress, 'S.Score', 'G2M.Score'), verbose = F)
    }
    DefaultAssay(seu) <- "SCT"
  } else {
    message("\t - Log normalization")
    seu <- NormalizeData(seu, verbose = F) %>% 
           FindVariableFeatures(nfeatures = params$variable.features.n, verbose = F) %>% 
           ScaleData(vars.to.regress = params$var.to.regress, verbose = F)
    
    seu <- add_cell_cycle(seu, DefaultAssay(seu))
    
    # Re-scale with cell cycle if needed
    if(params$cc) {
      message("\t - Log normalization  - regressing out cell cycling genes")
      seu <- ScaleData(seu, vars.to.regress = c(params$var.to.regress, "S.Score", "G2M.Score"))
    }
  }
  
  # Run dimensionality reduction and clustering
  message("\t - Running PCA > Clusters > UMAP")
  seu <- RunPCA(seu, verbose = F) %>% 
         FindNeighbors(dims = 1:params$n.dim, reduction = "pca", verbose = F) %>%  
         FindClusters(cluster.name = "unintegrated_clusters", verbose = F) %>% 
         RunUMAP(dims = 1:params$n.dim, reduction = "pca", reduction.name = "UMAP", verbose = F)
  
  pdf(file.path(dirs$plots, paste0(params$project.prefix, "-umap-non-integrated.pdf")))
    p <- DimPlot(seu, group.by = "unintegrated_clusters") + NoLegend()
    p2 <- DimPlot(seu, group.by = "sample") + NoLegend()
    print(p)
    print(p2)
  dev.off()

  saveRDS(seu, file = file.path(dirs$RDS.files, 
                               paste0(params$project.prefix, "-processed-obj.RDS")))
  message("processing done.")
  return(seu)
}

# Function to run integration
run_integration <- function(seu, params, dirs) {
  message("Running integration...")
  if(params$normalization == "SCT") {
    message("\t - Integrate Layers w/: ", params$integration)
    seu <- IntegrateLayers(object = seu, normalization.method = params$normalization, 
                          orig.reduction = "pca", method = paste0(params$integration,"Integration"), 
                          new.reduction = "integrated", verbose = F)
    message("\t - Re processed: Clusters > UMAP")
    seu <- FindNeighbors(seu, dims = 1:params$n.dim, reduction = "integrated") %>%  
           FindClusters(cluster.name = "integrated_clusters", verbose = F) %>% 
           RunUMAP(reduction = "integrated", dims = 1:params$n.dim, verbose = FALSE, 
                  reduction.name = "UMAP_integrated")
  } else {
    if(params$integration == "Harmony") {
      message("\t - Run Harmony integration")
      seu <- RunHarmony(seu, group.by.vars = "sample", reduction.save = "integrated")
    } else {
      message("\t - Run ",  params$integration, " integration")
      seu[['RNA']] <- split(seu[['RNA']], seu$sample)
      seu <- IntegrateLayers(object = seu, method = params$integration, 
                            orig.reduction = "pca", new.reduction = "integrated", 
                            verbose = FALSE)
    }
    message("\t - Re processed: Clusters > UMAP")
    seu <- FindNeighbors(seu, reduction = "integrated", dims = 1:params$n.dim.integration) %>% 
           FindClusters(cluster.name = "integrated_clusters", resolution = params$cluster_resolutions) %>%
           RunUMAP(reduction = "integrated", dims = 1:params$n.dim.integration, 
                  reduction.name = "UMAP_integrated")
  }

  pdf(file.path(dirs$plots, paste0(params$project.prefix, "-umap-integrated.pdf")))
    p <- DimPlot(seu, group.by = "integrated_clusters", reduction = "UMAP_integrated")
    p2 <- DimPlot(seu, group.by = "sample", reduction = "UMAP_integrated")
    print(p)
    print(p2)
  dev.off()
  saveRDS(seu, file = file.path(dirs$RDS.files, 
                               paste0(params$project.prefix, "-integrated-obj.RDS")))
  message("integration done.")
  return(seu)
}

# Function to run annotation
run_annotation <- function(seu, params, dirs) {
  message("Running cell type annotation...")
  if (!file.exists(params$reference)) {
    warning("Reference file not found. Skipping annotation step.")
    return(seu)
  }

  ref <- readRDS(params$reference)
  
  # Handle normalization types for reference mapping
  if(params$normalization == "SCT") {

    if(!"SCT" %in% SeuratObject::Assays(ref)) {
      message("\t- Different normalization between reference and query - need to normalize the query again")
      DefaultAssay(seu) <- "RNA"
      seu <- JoinLayers(seu)
      seu <- NormalizeData(seu)
      seu <- find_and_transfer_labels(ref, seu, out.path, ref.metadata.col = params$ref.metadata.col)
    } else {
      seu <- find_and_transfer_labels(seu, ref, out.path, ref.metadata.col = params$ref.metadata.col, 
                                    normalization = "SCT") 
    }
  } else {
    if("SCT" %in% SeuratObject::Assays(ref)) {
      message("\t- Different normalization between reference and query - need to normalize the query again")
      seu.sct <- SCTransform(seu, verbose = F)
      seu <- find_and_transfer_labels(ref, seu.sct, out.path, ref.metadata.col = params$ref.metadata.col, 
                                    normalization = "SCT")
    } else {
      seu <- find_and_transfer_labels(seu, ref, out.path, ref.metadata.col = params$ref.metadata.col)
    }
  }
  
  pdf(file.path(dirs$plots, paste0(params$project.prefix, "-annotated-obj.pdf")))
    p <- DimPlot(seu, group.by = "predicted.id")
    print(p)
    if("integrated" %in% Reductions(seu)){
      p2 <-  DimPlot(seu, group.by = "predicted.id", reduction = "UMAP_integrated")
      print(p2) 
    }
  dev.off()
  saveRDS(seu, file = file.path(dirs$RDS.files, 
                               paste0(params$project.prefix, "-annotated-obj.RDS")))
  message("annotation done.")
  return(seu)
}

# Function to run DEA
run_dea <- function(seu, params) {
  message("Running differential expression analysis...")
  # DEA code would go here
  # For example:
  # markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # write.csv(markers, file.path(params$output_dir, "all_markers.csv"))
  message("DEA done.")
  return(seu)
}

# Main execution pipeline
main <- function() {
  # Ensure output directory exists
  in.path <- file.path(params$output, "samples")
  out.path <- file.path(params$output, "merged")
  dirs <- list(plots = file.path(out.path, "plots"), RDS.files = file.path(out.path, "RDS-files"))

  print(paste0("output dir: ",out.path))
  dir.create(out.path, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirs$plots, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirs$RDS.files, recursive = TRUE, showWarnings = FALSE)

  # Initialize objects
  seurat.list <- NULL
  seu <- NULL
  
  # Run pipeline steps based on user selection
  if("filter" %in% pipeline.to.run) {
    message("\n\n++++++++ Loading data and filtering samples ++++++++++")
    seurat.list <- run_filter(in.path, dirs, params)
  }
  
  if("merge" %in% pipeline.to.run && !is.null(seurat.list)) {
    message("\n\n++++++++ Merging ++++++++++")
    seu <- run_merge(seurat.list, dirs)
  }
  
  if("qc" %in% pipeline.to.run && !is.null(seurat.list)) {
    message("\n\n++++++++ QC ++++++++++")
    seu <- run_qc(seu, params, dirs)
  }

  if("processing" %in% pipeline.to.run && !is.null(seu)) {
    message("\n\n++++++++ Processing ++++++++++")
    seu <- run_processing(seu, params, dirs)
  }
  
  if("integration" %in% pipeline.to.run && !is.null(seu)) {
    message("\n\n++++++++ Integration ++++++++++")
    seu <- run_integration(seu, params, dirs)
  } else if(!is.null(seu) && params$normalization == "SCT") {
    seu <- JoinLayers(seu, assay = "RNA")
  }
  
  if("annotation" %in% pipeline.to.run && !is.null(seu)) {
    message("\n\n++++++++ Annotation ++++++++++")
    seu <- run_annotation(seu, params, dirs)
  }
  
  if("DEA" %in% pipeline.to.run && !is.null(seu)) {
    message("\n\n++++++++ DEA ++++++++++")
    seu <- run_dea(seu, params)
  }
  

  
   message("Analysis complete. Final object saved in", out.path)
}

# Execute main function
main()
