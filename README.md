# scRNA-seq Processing Pipeline

This pipeline processes single-cell or single-nucleus RNA sequencing data and merges the results for further analysis.

## Table of Contents

1. [Introduction](#introduction)
2. [Setup](#setup)
3. [Execution](#execution)
4. [Parameters](#parameters)
5. [Running the Pipeline Independently](#running-the-pipeline-independently)

---

## Introduction

This pipeline is designed to process RNA-seq data from single-cell or single-nucleus samples. It allows for processing multiple samples independently and merging them into a single Seurat object for downstream analysis.

---

## Setup

### 1. Clone the Repository

To start using the pipeline, clone the repository to your local machine:

```bash
git clone <repository_url>
```

### 2. Modify the Sample Sheet

The pipeline requires a sample sheet in `.csv` format. The sample sheet should contain the following two columns:

- **Sample ID**: A unique identifier for each sample.
- **Cell Ranger Output Path**: The file path to the Cell Ranger output for each sample.

Make sure to adjust the paths in the sample sheet as needed for your environment.

### 3. Modify `params.yaml`

The `params.yaml` file contains various configuration parameters for the pipeline. The following parameters are mandatory:

- **pipeline**: Specify the steps you want to run, separated by commas (no spaces). To run the entire pipeline, use `full` (e.g., `pipeline: full`).
- **type.sequencing**: make sure you set this argument to either snRNA or scRNA

---

## Execution

Once the setup is complete, you can execute the pipeline using the following command:

```bash
sbatch scRNA-pipeline.sh
```

### Step-by-Step Execution Process

1. **Sample Processing**: The pipeline will first process each sample independently based on the information in the sample sheet. This will include quality control, filtering, and other preprocessing steps.
   
2. **Metric Summary**: After processing all samples, a metric summary will be generated to help evaluate the quality of the data for each sample.

3. **Merging Process**: Once all samples are processed successfully (i.e., without premature stops or errors), the pipeline will merge the samples. Only completed samples will be included in the merging process. The number of cells in each sample is also evaluated, and the default threshold is 200 cells. This threshold can be modified by adjusting the `min_cells` parameter in the `params.yaml`.

---

## Parameters

### Required Parameters in `params.yaml`

- **pipeline**: Define the pipeline steps you want to run (comma-separated, no spaces). Use `full` to run all steps.
- **min_cells**: The minimum number of cells that must be present in each sample to include it in the merge step. Default is `200`. You can modify this value in `params.yaml`.

### Optional Parameters for Running the Entire Pipeline

- **reference**: Path to a reference RDS file that contains a Seurat object for annotation.
- **ref.metadata.col**: The specific column in the reference metadata to use for annotation transfer.

---

## Running the Pipeline Independently

If you wish to run just the merging and processing steps (i.e., skip the individual sample processing), you can do so by following these steps:

1. Move the `processing.sh` file from the `scripts` directory into the previous directory.
   
2. Modify `params.yaml` as needed. If you skip the loading step (`load`), ensure that you provide a pre-processed Seurat object RDS file by specifying the path to it under `RDS.file` in `params.yaml`.

---

## Additional Notes

- Make sure all necessary files (e.g., sample sheet, Cell Ranger outputs, and Seurat reference) are properly linked and accessible.
