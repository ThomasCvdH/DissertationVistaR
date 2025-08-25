# DissertationVistaR: Visualisation of Isoform, Splicing and Transcript Abundance pipeline.

This repository is a pipeline I built for my MSc Bioinformatics project, looking at gene expression and transcript structures in Arabidopsis thaliana to identify temperature-dependant alternative splicing. 
## Included files

- **main.R**  
  This is the main script. It loads the data, does some wrangling, and makes lots of plots (dashboards, diversity comparisons, grids) for a given list of genes. You can change which dataset you want to look at at the top.

- **functions.R**  
  All the main functions live here. There’s stuff for parsing data, plotting, calculating stats, and picking the most variable transcripts using an elbow-point method (L-BowR).

- **Timeseries_cross_day.R**  
  If you want to compare transcript diversity between two days in a time series, this script makes a combined dataframe for the overlay.

- **MybDomainFinder.R**  
  This script finds Myb-like DNA-binding domains, or any other protein domains, in transcripts using the InterProScan API. It maps the protein domain back to the genome and saves the results as a GTF file.

## How To Use

1. **Unzip fasta,gtf, and csv files**  
   Unzip the files so they can be used in the pipeline. 
   They need to be at the same directory level as the scripts.

2. **Install Packages**  
   The scripts will try to install any R packages you need, but you’ll need an internet connection for this and for the InterProScan API.

3. **Get Your Data Ready**  
   - Put your transcriptome FASTA and GTF annotation files in the same folder as the scripts.
   - Make sure you have the CSV files (see the filenames in `main.R`).

4. **Run main.R**  
   - At the top, set `dataset` to `"adam"`, `"coolLL2"`, or `"timeseries"`.
   - Run the script. It’ll make figures and dashboards for the genes in the list.
   - Figures are saved in the `figures/<dataset>/` folder.

5. **Domain Annotation (MybDomainFinder.R)**  
   - Run this script to find Myb domains in your transcripts.
   - It’ll save new GTF files with the domain info.

6. **Time Series Comparison (Timeseries_cross_day.R)**  
   - The functions here are automatically used to make dataframes for comparing across days if the timeseries datset is selected.

## Features

- Works with different datasets (just change the setting at the top).
- Picks the most variable transcripts automatically (L-BowR).
- Makes dashboards with expression, structure, and diversity plots.
- Can add protein domain info to transcript structures.

## Notes

- These scripts are for R (version 4.0 or higher).
- You’ll need internet for installing packages and using InterProScan.
- Output folders are made automatically.
- To change which genes are plotted, edit the `genes_to_plot` bit in `main.R`.

 ## Author
 Thomas Cornelius van der Hoven
 thomasvanderhoven@gmail.com




