# Consolidate Sample Metrics CSV Files

#
# This script scans through all sample output folders and combines
# their metrics CSV files into one master CSV file.

# Load required libraries
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)

# Define the base output path (same as in your processing script)
# base_path <- "/gpfs/data/tsirigoslab/home/cather01/projects/LUAD/scratch/pipeline-test"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Error: Please provide the output directory as an argument.")
}
base_path <- args[1]


# Function to safely read a CSV file, returning NA if file doesn't exist or has issues
safe_read_csv <- function(file_path) {
  tryCatch({
    if (file.exists(file_path)) {
      df <- read_csv(file_path, show_col_types = FALSE)
      return(df)
    } else {
      message(paste("File not found:", file_path))
      return(NULL)
    }
  }, error = function(e) {
    message(paste("Error reading file", file_path, ":", e$message))
    return(NULL)
  })
}
message(base_path)
# Get all subdirectories in the base path
sample_dirs <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
message(sample_dirs)
message(paste("Found", length(sample_dirs), "sample directories"))

# Process each sample directory
all_metrics <- data.frame()
successful_samples <- 0
failed_samples <- 0

for (dir_path in sample_dirs) {
  message(dir_path)
  sample_id <- basename(dir_path)
  message(file.path(dir_path, "output", "metrics.csv"))
  metrics_file <- file.path(dir_path, "output", "metrics.csv")
  
  metrics_df <- safe_read_csv(metrics_file)
  
  if (!is.null(metrics_df)) {
    # Add the sample ID if it's not already included
    if (!"sample_id" %in% colnames(metrics_df)) {
      metrics_df$sample_id <- sample_id
    }
    
    # Append to the master dataframe
    if (nrow(all_metrics) == 0) {
      all_metrics <- metrics_df
    } else {
      # Handle potential column mismatches
      missing_cols <- setdiff(colnames(all_metrics), colnames(metrics_df))
      for (col in missing_cols) {
        metrics_df[[col]] <- NA
      }
      
      missing_cols_in_all <- setdiff(colnames(metrics_df), colnames(all_metrics))
      for (col in missing_cols_in_all) {
        all_metrics[[col]] <- NA
      }
      
      all_metrics <- bind_rows(all_metrics, metrics_df)
    }
    
    successful_samples <- successful_samples + 1
    message(paste("Processed metrics for sample:", sample_id))
  } else {
    failed_samples <- failed_samples + 1
    message(paste("No metrics found for sample:", sample_id))
  }
}

# Write the consolidated metrics to a file
if (nrow(all_metrics) > 0) {
  output_file <- file.path(base_path, "qc_metrics.csv")
  write_csv(all_metrics, output_file)
  
  
  message(paste("Successfully qc metrics from", successful_samples, "samples"))
  message(paste("Failed to find metrics for", failed_samples, "samples"))
  message(paste("Output written to:", output_file))
  
  # Generate summary statistics
  summary_stats <- all_metrics %>%
    summarize(
      total_samples = n(),
      completed_samples = sum(status == "completed", na.rm = TRUE),
      failed_samples = sum(status == "failed", na.rm = TRUE),
      avg_cells_pre_filtering = mean(cell.count.prefilter, na.rm = TRUE),
      avg_cells_post_filtering = mean(cell.count.postfilter, na.rm = TRUE),
      avg_cells_post_doublet = mean(cell.count.postdoublet, na.rm = TRUE),
      percent_cell_kept_postfilter = mean(percent.cells.kept.postfilter , na.rm = TRUE),
      avg_percent_doublets = mean(percent.doublets, na.rm = TRUE),
      avg_median_features = mean(median.features, na.rm = TRUE),
      avg_median_counts = mean(median.counts, na.rm = TRUE),
      avg_processing_time = mean(runtime, na.rm = TRUE)
    )
  
  # Write summary statistics
  summary_file <- file.path(base_path, "metrics_summary.csv")
  write_csv(summary_stats, summary_file)
  message(paste("Summary statistics written to:", summary_file))
  
  # Create a processing status summary
  if ("completed_step" %in% colnames(all_metrics)) {
    step_summary <- all_metrics %>%
      group_by(completed_step) %>%
      summarize(count = n()) %>%
      arrange(count)
    
    step_file <- file.path(base_path, "processing_step_summary.csv")
    write_csv(step_summary, step_file)
    message(paste("Processing step summary written to:", step_file))
  }
} else {
  message("No metrics found in any sample directory")
}
print(all_metrics)
# Create a simple HTML report of the metrics
if (nrow(all_metrics) > 0) {
  html_report <- file.path(base_path, "metrics_report.html")
  
  # Function to create an HTML table from a dataframe
  df_to_html_table <- function(df) {
    header <- paste0("<tr>", paste0("<th>", colnames(df), "</th>", collapse = ""), "</tr>")
    rows <- apply(df, 1, function(row) {
      paste0("<tr>", paste0("<td>", row, "</td>", collapse = ""), "</tr>")
    })
    body <- paste(rows, collapse = "")
    return(paste0("<table border='1'>", header, body, "</table>"))
  }
  
  # Create a simple HTML report
  html_content <- paste0(
    "<!DOCTYPE html>
    <html>
    <head>
      <title>LUAD Sample Processing Metrics</title>
      <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1, h2 { color: #333366; }
        table { border-collapse: collapse; margin: 15px 0; }
        th { background-color: #f2f2f2; padding: 8px; text-align: left; }
        td { padding: 8px; }
        tr:nth-child(even) { background-color: #f9f9f9; }
      </style>
    </head>
    <body>
      <h1>LUAD Sample Processing Metrics</h1>
      <p>Report generated on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "</p>
      
      <h2>Summary Statistics</h2>
      ", df_to_html_table(summary_stats), "
      
      <h2>Sample Metrics</h2>
      <p>Total samples processed: ", nrow(all_metrics), "</p>
      <p>Click <a href='qc_metrics.csv'>here</a> to download the complete metrics CSV.</p>
      
      <h3>Processing Status</h3>
      <ul>
        <li>Completed: ", sum(all_metrics$status == "completed", na.rm = TRUE), "</li>
        <li>Failed: ", sum(all_metrics$status == "failed", na.rm = TRUE), "</li>
        <li>Running: ", sum(all_metrics$status == "running", na.rm = TRUE), "</li>
      </ul>
      
      <h2>Samples</h2>
      <table border='1'>
        <tr>
          <th>Sample ID</th>
          <th>Status</th>
          <th>Last Completed Step</th>
          <th>Cells (Final)</th>
          <th>Doublets (%)</th>
          <th>Processing Time (min)</th>
        </tr>",
        paste(apply(all_metrics, 1, function(row) {
          paste0("<tr>
            <td>", row["sample"], "</td>
            <td>", row["status"], "</td>
            <td>", row["step"], "</td>
            <td>", row["median.features"], "</td>
            <td>", row["median.counts"], "</td>
            <td>", row["runtime"], "</td>
          </tr>")
        }), collapse = ""),
      "</table>
    </body>
    </html>"
  )
  
  # Write the HTML report
  writeLines(html_content, html_report)
  message(paste("HTML report written to:", html_report))
}

message("Consolidation complete!")