# Install required packages if not already installed
required_packages <- c(
  "dplyr", "tidyr", "stringr", "ggplot2", "tibble",
  "ggtranscript", "rtracklayer", "patchwork", "svglite", "viridis"
)
# Check which packages are not already installed
new_packages <- required_packages[!(
  required_packages %in% installed.packages()[, "Package"]
)]
# Install any missing packages
if (length(new_packages)) install.packages(new_packages)

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)
library(ggtranscript)
library(rtracklayer)
library(patchwork)
library(viridis)
# Install ggtranscript from github
remotes::install_github("dzhang32/ggtranscript") #needed for plotting transcript structures

# Load Custom Functions
# Source the custom functions file
source("functions.R")
# GTF Import
# Define the path to the GTF file
gtf_path <- file.path("ALL_combined_annotations3.gtf")
# Import the GTF file and convert it to a tibble
gtf_df <- rtracklayer::import(gtf_path) %>% as_tibble()

#  1. Set Dataset Context
# Define the dataset to be used (options: "adam", "coolLL2", "timeseries")
dataset <- "timeseries" # Options: "adam", "coolLL2", "timeseries"

#  2. Set up output directory for figures 
# Define the output directory for figures
figures_dir <- file.path("figures", dataset)
# Create the output directory if it does not exist
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir)

#  3. Load Data Based on Context
# Load the data based on the selected dataset context
data_list <- switch(dataset,
  "adam" = parse_data("Transcript_TPM_adam_TRIMMED.csv"),
  "coolLL2" = parse_data("Transcript TPMcoolLL.csv"),
  "timeseries" = parse_data("Transcript_TPM_timeSeries_trimmed_col0.csv")
)
df <- data_list$df
averaged_df <- data_list$averaged_df
print(averaged_df)

# If the dataset is "timeseries", create overlayed dataframe for Shannon Diversity
if (dataset == "timeseries") {
  source("Timeseries_cross_day.R") # Source the timeseries-specific functions

  timeseries_df <- create_timeseries_df(
    raw_df = df,
    dataset = dataset,
    # choose days to comare for overlay
    day1 = 1,
    day2 = 3
  )
}

#  4. Usage Examples

# Generate and print a dashboard for a single gene
rve2_dashboard <- plot_expression_dashboard(
  gene_id = "AT5G37260",
  gene_name = "RVE2",
  gtf_data = gtf_df,
  raw_df = df,
  avg_df = averaged_df,
  dataset = dataset
  #manual_cutoff = 19 # Uncomment if you want to set a manual cutoff
)
print(rve2_dashboard)

# Generate and print a shannon diversity comparison for a single gene
if (dataset == "timeseries") {
  comparison_plot <- plot_diversity_comparison(
    tair_code = "AT3G09600",
    gene_name = "RVE8",
    raw_df = timeseries_df,
    dataset = "timeseries_df"
  )
} else {
  comparison_plot <- plot_diversity_comparison(
    tair_code = "AT3G09600",
    gene_name = "RVE8",
    raw_df = df,
    dataset = dataset,
    day_to_compare = 3 # Specify day for coolLL2 dataset
  )
}
print(comparison_plot)

# Batch process list of genes for dashboards, grids, Statistical comparison.

genes_to_plot <- tibble::tribble(
  ~gene_name, ~tair_code,
  #"RVE1", "AT5G17300 ",
  "RVE2", "AT5G37260",
  "RVE8", "AT3G09600",
  "LHY", "AT1G01060",
  "CCA1", "AT2G46830",
  "PHYB", "AT2G18790",
  "PRR9", "AT2G46790",
  "PRR7", "AT5G02810",
  "PRR5", "AT5G24470",
  "TOC1", "AT5G61380",
  "LUX", "AT3G46640",
  "ELF3", "AT2G25930",
  "COP1", "AT2G32950",
  "LNK4", "AT5G06980",
  #"LWD1", "AT1G12910", "LWD2", "AT3G26640", "TTG1", "AT5G24520",
  #"COR27", "AT5G42900", "COR28", "AT4G33980", "PIF1", "AT2G20180",
  #"PIF3", "AT1G09530", "PIF4", "AT2G43010", "PIF5", "AT3G59060",
  #"TCP2", "AT3G27010", "TCP21", "AT5G08330",  "TCP22", "AT1G72010",
  #"RVE3", "AT1G01520", "RVE4", "AT5G02840", "RVE5", "AT4G01280",
  #"RVE6", "AT5G52660", "RVE7", "AT1G18330", "ELF4", "AT2G40080",
  #"ZTL", "AT5G57360", "GI", "AT1G22770", "BOA", "AT5G59570",
  #"LNK1", "AT5G64170", "LNK2", "AT3G54500", "LNK3", "AT3G12320",
)

#  Batch Dashboard Generation
# Loop through each gene in the list and generate expression dashboards
for (i in 1:nrow(genes_to_plot)) {
  current_gene_name <- genes_to_plot$gene_name[i]
  current_gene_id <- str_trim(genes_to_plot$tair_code[i])

  tryCatch(
    {
      cat(paste(" Generating dashboard for", current_gene_name, "\n"))
      dashboard_plot <- plot_expression_dashboard(
        gene_id = current_gene_id,
        gene_name = current_gene_name,
        gtf_data = gtf_df,
        raw_df = df,
        avg_df = averaged_df,
        dataset = dataset # Pass dataset
      )
    },
    error = function(e) {
      cat(paste(
        "Could not generate dashboard for",
        current_gene_name, ":",
        e$message, "\n\n"
      ))
    }
  )
}

# Batch Shannon Diversity Comparison Plot batch
# Loop through each gene and generate diversity comparison plots
for (i in 1:nrow(genes_to_plot)) {
  current_gene_name <- genes_to_plot$gene_name[i]
  current_gene_id <- str_trim(genes_to_plot$tair_code[i])

  tryCatch(
    {
      cat(paste(" Generating SD comparison for ", current_gene_name, "\n"))
      if (dataset == "timeseries") {
        comparison_plot <- plot_diversity_comparison(
          tair_code = current_gene_id,
          gene_name = current_gene_name,
          raw_df = timeseries_df,
          dataset = "timeseries_df"
        )
      } else {
        comparison_plot <- plot_diversity_comparison(
          tair_code = current_gene_id,
          gene_name = current_gene_name,
          raw_df = df,
          dataset = dataset,
          day_to_compare = 3
        )
      }
    },
    error = function(e) {
      cat(paste(
        "Could not generate SD comparison plot for",
        current_gene_name, ":",
        e$message, "\n\n"
      ))
    }
  )
}

# Generating Grid Plots
tryCatch(
  {
    cat(" Generating  Grid Plot ")
    # Generate total expression grid plot 
    total_expression_plot<- plot_total_expression_grid(
      df,
      averaged_df,
      genes_to_plot,
      dataset
    )
    # Save total expression grid plot as SVG
    ggsave(
      file.path(figures_dir, paste0(
        "total_expression_grid_",
        toupper(dataset),
        ".svg"
      )),
      plot = total_expression_plot,
      width = 24,
      height = 10
    )
    message(" grid plots saved successfully.")
  },
  error = function(e) {
    cat("Error generating  grid plots:", e$message, "\n")
  }
)


plot_list <- list() # Initialise an empty list to store each plot

genes_to_plot<-arrange(genes_to_plot, gene_name) # If you want to order the plots alphabetically

# Loop through each row of the genes_to_plot tibble
for (i in 1:nrow(genes_to_plot)) {
  current_gene_name <- genes_to_plot$gene_name[i]
  current_tair_code <- genes_to_plot$tair_code[i]

  cat(" Generating plot for:", current_gene_name, "\n")

  tryCatch(
    {
      # Call function to generate one plot for the current gene
      single_plot <- plot_individual_transcripts(
        tair_code = current_tair_code,
        gene_name = current_gene_name,
        avg_df = averaged_df,
        raw_df = df,
        dataset = dataset
      )
      # Add the successfully created plot to the list
      plot_list[[current_gene_name]] <- single_plot
    },
    error = function(e) {
      # If the function fails for a gene, print an error message
      cat(paste(
        "Could not generate plot for",
        current_gene_name, ":",
        e$message, "\n\n"
      ))
    }
  )
}

# After the loop, check if any plots were successfully created
if (length(plot_list) > 0) {

  # Add lettered tags (A, B, C...) to each gene plot
  tagged_plots <- purrr::map2(
    plot_list,
    seq_along(plot_list),
    function(p, index) {
      p + labs(tag = paste0(LETTERS[index], ")")) +
        theme(plot.tag = element_text(size = 18, face = "bold"))
    }
  )

  final_grid <- patchwork::wrap_plots(tagged_plots, ncol = 3, nrow = 2)

  # Print the final, stitched plot
  print(final_grid)

  # Save the final plot to a file
  ggsave(
    file.path(figures_dir, "Combined_Transcript_Structure_Grid.svg"),
    plot = final_grid,
    width = 28,
    height = 15,
    units = "in"
   )
  cat("\nPlot saved successfully as 'Combined_Transcript_Structure_Grid.svg'\n")

} else {
  cat("No plots were generated successfully.\n")
}
