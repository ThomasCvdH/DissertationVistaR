required_packages <- c(
  "rtracklayer", "GenomicRanges", "Biostrings", "GenomicFeatures", "httr",
  "readr", "dplyr", "tibble"
)
# Check which packages are not already installed
new_packages <- required_packages[!(
  required_packages %in% installed.packages()[, "Package"]
)]
# Install any missing packages
if (length(new_packages)) install.packages(new_packages)
# Load required libraries for genomic data handling and API requests
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)
library(GenomicFeatures)
library(httr)
library(readr)
library(dplyr)
library(tibble)

# Specify email for InterProScan API (required for job submission)
email <- NULL
if (is.null(email)) {
  stop("Please specify your email address in the 'email' variable to use the InterProScan API.")
}

# Specify input files
transcriptome_file <- "atRTD3_29122021.fa"# FASTA file containing transcript sequences
annotation_file <- "ALL_combined_annotations3.gtf" # GTF file containing transcript/gene annotations

# InterProScan settings
target_domain_id <- "IPR017930" # InterPro ID for Myb-like DNA-binding domain

# List of genes to analyse (gene name and gene ID)
genes_to_plot <- tibble::tribble(
  ~gene_name, ~gene_id,
  "LHY", "AT1G01060"
)

# Load data files and check existence
if (!file.exists(transcriptome_file)) stop(paste("FASTA file not found:", transcriptome_file))
if (!file.exists(annotation_file)) stop(paste("GTF file not found:", annotation_file))
transcript_seqs <- readDNAStringSet(transcriptome_file) # Load transcript sequences

cat("Importing GTF as a simple table...\n")
# Import the GTF directly as a GRanges object
gff <- rtracklayer::import(annotation_file)

# Remove any whitespace from gene IDs
genes_to_plot$gene_id <- trimws(genes_to_plot$gene_id)

# Find all transcript IDs associated with the target genes
target_tx_ids_in_gff <- gff[gff$gene_id %in% genes_to_plot$gene_id]$transcript_id
all_target_tx_ids <- unique(na.omit(target_tx_ids_in_gff))
cat("Found", length(all_target_tx_ids), "transcripts for", nrow(genes_to_plot), "genes.\n")

# Create an empty GRangesList to store results for all transcripts
all_results_list <- GRangesList()

# Loop over each transcript ID to process
for (target_transcript_id in all_target_tx_ids) {
  cat("Processing Transcript:", target_transcript_id, "\n")

  #  Create Gene Model from GTF Table 
  target_transcript_sequence <- transcript_seqs[[target_transcript_id]]
  if (is.null(target_transcript_sequence)) {
    cat("  -> INFO: Transcript not found in FASTA file. Skipping.\n")
    next
  }

  # Filter the GFF table for the current transcript
  current_tx_features <- gff[gff$transcript_id == target_transcript_id]
  if (length(current_tx_features) == 0) {
    cat("  -> INFO: No features found for this transcript in the GTF. Skipping.\n")
    next
  }

  # Extract CDS features and sort them to ensure correct order
  cds_features <- current_tx_features[current_tx_features$type == "CDS"]
  if (length(cds_features) == 0) {
    cat("  -> INFO: No CDS features found for this transcript. Skipping.\n")
    next
  }
  cds_genomic_coords <- sort(cds_features)

  # For reverse strand genes, exons are transcribed from high-to-low coordinates
  is_reverse_strand <- as.character(strand(cds_features[1])) == "-"
  cds_genomic_coords <- sort(cds_features, decreasing = is_reverse_strand)

  # Extract 5' UTR features to calculate offset for CDS extraction
  five_prime_utr_features <- current_tx_features[current_tx_features$type == "five_prime_utr"]

  # Calculate lengths of 5' UTR and CDS
  five_utr_width <- sum(width(five_prime_utr_features)) # Total length of 5' UTR
  cds_length <- sum(width(cds_genomic_coords)) # Total length of CDS

  # Check that CDS length is valid for translation
  if (cds_length == 0 || cds_length %% 3 != 0) {
    cat("  -> INFO: CDS length is zero or not a multiple of 3. Skipping.\n")
    next
  }

  # Extract CDS sequence from transcript sequence
  cds_only_sequence <- subseq(target_transcript_sequence,
    start = five_utr_width + 1,
    width = cds_length
  )
  protein_sequence <- as.character(translate(cds_only_sequence)) # Translate to protein sequence
  protein_sequence <- gsub("\\*", "", protein_sequence) # Remove stop codons

  # Submit protein sequence to InterProScan API to find domains
  iprscan_url <- "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
  submit_job <- POST(
    iprscan_url,
    body = list(sequence = protein_sequence, email= email),
    encode = "form"
  )

  # Check if submission was successful
  if (status_code(submit_job) != 200) {
    cat("WARN: Failed to submit sequence to InterProScan API. Skipping.\n")
    next
  }
  job_id <- content(submit_job, as = "text")

  # Poll InterProScan API for job status until finished or error
  repeat {
    status_url <- paste0("https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/", job_id)
    status_req <- GET(status_url)
    job_status <- content(status_req, as = "text")
    if (job_status %in% c("FINISHED", "ERROR", "FAILURE", "NOT_FOUND")) break
    Sys.sleep(5) # Wait before polling again
  }

  # Check if job finished successfully
  if (job_status != "FINISHED") {
    cat("WARN: Job failed or finished with error. Status:", job_status, ". Skipping.\n")
    next
  }

  # Retrieve results from InterProScan
  result_url <- paste0("https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/", job_id, "/tsv")
  results_text <- content(GET(result_url), as = "text")
  if (nchar(results_text) == 0) {
    cat("INFO: No results returned from InterProScan. Skipping.\n")
    next
  }

  # Parse results into a data frame
  col_names <- c(
    "protein_id", "md5", "length", "analysis", "signature_acc",
    "description", "start_pos", "end_pos", "score", "status", "date", "interpro_acc"
  )
  results_df <- read_tsv(results_text, col_names = col_names, col_types = cols())

  # Filter for the target domain
  myb_domain <- results_df %>% filter(interpro_acc == target_domain_id)

  if (nrow(myb_domain) == 0) {
    cat("  -> INFO: Domain", target_domain_id, "not found. Skipping.\n")
    next
  }

  # Get domain start/end positions in protein sequence
  domain_start_aa <- myb_domain$start_pos[1]
  domain_end_aa <- myb_domain$end_pos[1]

  # Map domain coordinates from protein sequence to genome
  domain_start_nt <- (domain_start_aa - 1) * 3 + 1
  domain_end_nt <- domain_end_aa * 3
  cumulative_exon_lengths <- cumsum(width(cds_genomic_coords))
  # Function to map CDS nucleotide position to genomic coordinates
  mapCdsToGenomic <- function(cds_nt_pos, cds_exons, cumulative_lengths) {
    exon_index <- findInterval(cds_nt_pos, cumulative_lengths) + 1
    target_exon <- cds_exons[exon_index]
    pos_in_exon <- cds_nt_pos - c(0, cumulative_lengths)[exon_index]
    if (as.character(strand(target_exon)) == "+") {
      genomic_pos <- start(target_exon) + pos_in_exon - 1
    } else {
      genomic_pos <- end(target_exon) - pos_in_exon + 1
    }
    return(list(pos = genomic_pos, exon_index = exon_index, exon = target_exon))
  }
  start_info <- mapCdsToGenomic(domain_start_nt, cds_genomic_coords, cumulative_exon_lengths)
  end_info <- mapCdsToGenomic(domain_end_nt, cds_genomic_coords, cumulative_exon_lengths)

  # Create GRanges object for the domain parts

  gene_id <- current_tx_features$gene_id[1]

  # Create an empty list to store the GRanges object for each partial domain
  partial_domain_list <- GRangesList()

  # Loop through each exon the domain touches
  for (i in start_info$exon_index:end_info$exon_index) {
    current_exon <- cds_genomic_coords[i]
    is_reverse <- as.character(strand(current_exon)) == "-"

    # Determine the start coordinate of the domain piece on this exon
    if (i == start_info$exon_index) {
      partial_start <- start_info$pos
    } else {
      partial_start <- if (is_reverse) end(current_exon) else start(current_exon)
    }

    # Determine the end coordinate of the domain piece on this exon
    if (i == end_info$exon_index) {
      partial_end <- end_info$pos
    } else {
      partial_end <- if (is_reverse) start(current_exon) else end(current_exon)
    }

    # Create the GRanges object for this specific partial domain piece
    partial_domain_gr <- GRanges(
      seqnames = as.character(seqnames(current_exon)),
      ranges = IRanges(
        start = min(partial_start, partial_end),
        end = max(partial_start, partial_end)
      ),
      strand = strand(current_exon),
      source = "MybDomainFinder", type = "mybdomain_part", score = ".", phase = ".",
      gene_id = gene_id, transcript_id = target_transcript_id, exon_number = as.character(i)
    )
    partial_domain_list[[length(partial_domain_list) + 1]] <- partial_domain_gr
  }

  # Combine all partial domain GRanges into one object
  myb_domain_gr_detailed <- unlist(partial_domain_list)
  all_results_list[[target_transcript_id]] <- myb_domain_gr_detailed
  cat("SUCCESS: Found and processed MYB domain in", target_transcript_id, "\n")
}

# Consolidate and Export All Results

cat("All transcripts processed. Consolidating results...\n")

# If no domains found, print message
if (length(all_results_list) == 0) {
  cat("No MYB domains were found in any of the provided transcripts.\n")
} else {
  # Combine all results from the list into one GRanges object
  final_combined_results <- unlist(all_results_list)

  # Export the combined GTF 
  original_gr <- import(annotation_file, format = "gtf")
  combined_gr <- c(original_gr, final_combined_results)
  combined_output_file <- "ALL_combined_annotations4.gtf"
  rtracklayer::export(combined_gr, combined_output_file, format = "gtf")
  cat("Exported combined GTF with MYB domains to:", combined_output_file, "\n")

}