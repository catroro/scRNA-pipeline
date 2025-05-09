# Introduction
Pipeline to process single cell or single nucleus RNA seq data

# Setup
Process multiple sample independently and create a merged object to further analyze it
To start, git clone the repository
First thing is to modify the sample sheet. There is 2 columns, one with the sample ID, the other is the path to the cell ranger output.
Then modify the params.yaml. The following parameters are mandatory:
- pipeline: steps you want to run, comma separated without space. If you want to run the whole pipeline, write full


# Execute
You can now run the pipeline by doing 
sbatch scRNA-pipeline.sh
This will first run the processing for each sample independently that you have provided in the sample sheet.
Once all the samples are done, it will generate a metric summary 

It will now run the second part of the pipeline which merge and process the samples. Only the samples that have been completed (without premaatory stops or errors) will be merge. There is also a threshold for the minimum number of cells that has to be present in each samples (default = 200, can be modified by specifying min_cells in the params.yaml)
If you run the whole pipeline there is other parameter you want to modify in the params.yaml:
- reference: for the annotation step, need to provide a reference RDS file which contain a seurat object
- ref.metadata.col: which column from the reference you want to use to transfer the annotation

You can also choose to run the second part of the pipeline alone. To do so, move the processing.sh file from scripts into the previous directory.
You can then choose which step you want to run. If you don't use the load step, please provide a RDS file with a seurat object in the params.yaml under RDS.file