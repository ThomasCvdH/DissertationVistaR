

# Get Regex Patterns for Parsing Sample Names and Gene IDs
get_regex_patterns <- function(dataset) {
  # Define base patterns that are common
  base_patterns <- list(
    tair_code = "AT.G\\d+",
    transcript_suffix = "\\.\\d+$"
  )

  if (dataset == "adam") {
    # Patterns for adam dataset
    c(list(
      condition_group = "^X\\d+",
      time_h = "(?<=\\.X)(\\d+)(?=H)",
      rep = "\\.X\\d+$",
      pattern_type = "direct"
    ), base_patterns)
  } else if (dataset == "coolLL2") {
    # Patterns for coolLL2 dataset
    c(list(
      condition_group = "^X\\d+C",
      time_h = "(?<=\\.X)(\\d+)(?=h)",
      rep = "\\.rep\\d+$",
      pattern_type = "direct"
    ), base_patterns)
  } else if (dataset == "timeseries") {
    # Patterns for timeseries dataset
    c(list(
      condition_group = "X\\d+",
      timepoint = "(?<=T)\\d+",
      rep = "\\.rep\\d+$",
      pattern_type = "mapped"
    ), base_patterns)
  } else if (dataset == "timeseries_df") {
    c(list(
      condition_group = "^X\\d+",
      time_h = "(-?\\d+)(?=h)",
      rep = "_\\.rep\\d+$",
      pattern_type = "direct"
    ), base_patterns)
  } else {
    stop("Unknown dataset for regex patterns")
  }
}

# Get Time Mapping for timeseries Dataset
get_time_mapping <- function() {
  timepoints <- 1:26
  hours <- numeric(length(timepoints))
  periods <- character(length(timepoints))

  # Set the first timepoint
  hours[1] <- 0
  # Calculate hours for each timepoint based on experimental design
  for (t in 2:length(timepoints)) {
    if (t <= 11) {
      hours[t] <- hours[t - 1] + 3
    } else if (t %in% 12:13) {
      hours[t] <- hours[t - 1] + 1.5
    } else if (t > 13 && t <= 18) {
      hours[t] <- hours[t - 1] + 3
    } else if (t == 19) {
      hours[t] <- hours[t - 1] + 48
    } else {
      hours[t] <- hours[t - 1] + 3
    }
  }
  periods[timepoints <= 18] <- "Initial Period (T1-T18)"
  periods[timepoints > 18] <- "Post-Gap Period (T19-T26)"

  tibble::tibble(
    timepoint = as.character(timepoints),
    time_h = hours,
    period = factor(periods,
      levels = c("Initial Period (T1-T18)", "Post-Gap Period (T19-T26)")
    )
  )
}

# Get Plotting Aesthetics (Shapes, colours, Linetypes)
get_plot_aesthetics <- function(dataset) {
  if (dataset == "adam") {
    # Adam dataset aesthetics
    list(
      labels = c("15°C", "20°C"),
      values = c("X15" = "#377EB8", "X20" = "#E41A1C"),
      shapes = c("X15" = 21, "X20" = 24),
      linetypes = c("X15" = "solid", "X20" = "dashed")
    )
  } else if (dataset == "coolLL2") {
    # coolLL2 dataset aesthetics
    list(
      labels = c("12°C", "20°C"),
      values = c("X12C" = "#377EB8", "X20C" = "#E41A1C"),
      shapes = c("X12C" = 21, "X20C" = 24),
      linetypes = c("X12C" = "solid", "X20C" = "dashed")
    )
  } else if (dataset == "timeseries") {
    # timeseries dataset aesthetics
    list(
      labels = c("20°C", "4°C"),
      values = c("X20" = "#E41A1C", "X4" = "#377EB8"),
      shapes = c("X20" = 24, "X4" = 21),
      linetypes = c("X20" = "dashed", "X4" = "solid")
    )
  } else if (dataset == "timeseries_df") {
    # timeseries dataset aesthetics
    list(
      labels = c("20°C", "4°C"),
      values = c("X20" = "#E41A1C", "X4" = "#377EB8"),
      shapes = c("X20" = 24, "X4" = 21),
      linetypes = c("X20" = "dashed", "X4" = "solid")
    )
  } else {
    stop("Unknown dataset for plot aesthetics")
  }
}

# Parse Data File and Calculate Averaged Expression
# Reads a CSV file and calculates mean expression values for each gene and condition.
parse_data <- function(file_path) {
  regex <- get_regex_patterns(dataset)
  # Read CSV with row names as first column
  df <- read.csv(file_path, row.names = 1, check.names = FALSE)

  # Transform data:
  averaged_df <- df %>%
    # Convert row names to column for easier manipulation
    tibble::rownames_to_column(var = "gene") %>%
    # Convert wide format to long format (each row = one observation)
    pivot_longer(cols = -gene, names_to = "sample", values_to = "value") %>%
    # Remove replicate info from sample names to group by condition
    mutate(condition = str_remove(sample, regex$rep)) %>%
    # Calculate mean expression for each gene x condition combination
    group_by(gene, condition) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Convert back to wide format for plotting
    pivot_wider(names_from = "condition", values_from = "mean_value")

  list(df = df, averaged_df = averaged_df)
}

# Calculates Standard Deviation for Individual Transcripts
indv_trans_SD <- function(raw_df) {
  regex <- get_regex_patterns(dataset)
  # Read CSV with row names as first column
  df <- raw_df

  # Transform data:
  sd_df <- df %>%
    # Convert row names to column for easier manipulation
    tibble::rownames_to_column(var = "gene") %>%
    # Convert wide format to long format (each row = one observation)
    pivot_longer(cols = -gene, names_to = "sample", values_to = "value") %>%
    # Remove replicate info from sample names to group by condition
    mutate(condition = str_remove(sample, regex$rep)) %>%
    # Calculate mean expression for each gene x condition combination
    group_by(gene, condition) %>%
    summarise(
      sd_value = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Convert back to wide format for plotting
    pivot_wider(names_from = "condition", values_from = "sd_value")

  return(sd_df)
}


# Extract Time (hours) from Condition/Sample String
extract_time <- function(x, regex_patterns, time_map = NULL) {
  sapply(x, function(val) {
    if (!is.null(time_map)) { # Mapped time (e.g., timeseries)
      tp <- str_extract(val, regex_patterns$timepoint)
      if (is.na(tp)) {
        return(NA_real_)
      }
      time_row <- time_map %>% filter(timepoint == tp)
      if (nrow(time_row) == 0) {
        return(NA_real_)
      }
      return(time_row$time_h)
    } else { # Direct time extraction
      as.numeric(str_extract(val, regex_patterns$time_h))
    }
  })
}

# Plot Total Gene Expression with Error Bars
plot_gene_expression <- function(
    tair_code, gene_name,
    raw_df, avg_df, dataset) {
  regex <- get_regex_patterns(dataset)
  aesthetics <- get_plot_aesthetics(dataset)
  time_map <- if (regex$pattern_type == "mapped") get_time_mapping() else NULL
  print(time_map)

  if (!("gene" %in% colnames(raw_df))) {
    raw_df <- tibble::rownames_to_column(raw_df, var = "gene")
  }

  # Data transformation pipeline:
  plot_data <- raw_df %>%
    dplyr::filter(str_starts(gene, tair_code)) %>%
    pivot_longer(
      cols = where(is.numeric),
      names_to = "sample",
      values_to = "abundance"
    ) %>%
    mutate(
      condition_group = str_extract(sample, regex$condition_group),
      time_h = extract_time(sample, regex, time_map),
    ) %>%
    group_by(sample, condition_group, time_h) %>%
    summarise(
      total_abundance_per_rep = sum(abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(condition_group, time_h) %>%
    summarise(
      mean_abundance = mean(total_abundance_per_rep),
      sd_abundance = sd(total_abundance_per_rep),
      .groups = "drop"
    ) %>%
    filter(!is.na(time_h))


  # Generate plot differently based on dataset
  if (dataset == "timeseries") {
    gap_start_time <- 48
    display_gap_width <- 4 # The width of the blank space
    compression_amount <- (96 - gap_start_time) - display_gap_width

    plot_data <- plot_data %>%
      mutate(
        x_val = ifelse(time_h > gap_start_time, time_h - compression_amount, time_h),
        period = ifelse(time_h <= gap_start_time, "pre", "post")
      )

    original_breaks <- c(0, 24, 48, 96, 120, 147)
    compressed_breaks <- ifelse(original_breaks > gap_start_time, original_breaks - compression_amount, original_breaks)

    p <- ggplot(plot_data, aes(x = x_val, y = mean_abundance, colour = condition_group)) +
      geom_errorbar(aes(ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), width = 1.5, alpha = 0.5) +
      geom_line(aes(group = interaction(condition_group, period), linetype = condition_group), linewidth = 1) +
      geom_point(aes(shape = condition_group), fill = "white", size = 3.5, alpha = 0.8) +
      geom_vline(xintercept = c(gap_start_time, gap_start_time + display_gap_width), linetype = "dashed", colour = "grey40") +
      scale_x_continuous(name = "Time since cHilling (Hrs)", breaks = compressed_breaks, labels = original_breaks)
    } 
    else if (dataset == "coolLL2") {

    plot_data <- plot_data %>%
    mutate(time_h = time_h - 12) 
   
    p <- ggplot(plot_data, aes(x = time_h, y = mean_abundance, colour = condition_group)) +
      geom_errorbar(aes(ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), width = 1.5, alpha = 0.5) +
      geom_line(aes(linetype = condition_group), linewidth = 1) +
      geom_point(aes(shape = condition_group), fill = "white", size = 3.5, alpha = 0.8) +
      scale_x_continuous(name = "Time since chilling (Hrs)", breaks = scales::breaks_width(4)) 
  }
   else {
    p <- ggplot(plot_data, aes(x = time_h, y = mean_abundance, colour = condition_group)) +
      geom_errorbar(aes(ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), width = 1.5, alpha = 0.5) +
      geom_line(aes(linetype = condition_group), linewidth = 1) +
      geom_point(aes(shape = condition_group), fill = "white", size = 3.5, alpha = 0.8) +
      scale_x_continuous(name = "Time since dawn (Hrs)", breaks = scales::breaks_width(24))
  }

  # Add shared aesthetics and labels
  p +
    scale_colour_manual(name = "Condition", values = aesthetics$values, labels = aesthetics$labels) +
    scale_shape_manual(name = "Condition", values = aesthetics$shapes, labels = aesthetics$labels) +
    scale_linetype_manual(name = "Condition", values = aesthetics$linetypes, labels = aesthetics$labels) +
    labs(
      title = paste("Total Expression of", gene_name),
      subtitle = paste("TAIR ID:", tair_code),
      y = "Mean Expression (TPM)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, colour = "gray40"),
      legend.position = "bottom"
    )
}

# Plot Top Variable Individual Gene Transcripts
plot_individual_transcripts <- function(tair_code, gene_name, avg_df, raw_df, dataset, manual_cutoff = NULL) {
  regex <- get_regex_patterns(dataset)
  aesthetics <- get_plot_aesthetics(dataset)
  time_map <- if (regex$pattern_type == "mapped") get_time_mapping() else NULL

  error_bar_df <- indv_trans_SD(raw_df) %>%
    dplyr::filter(str_starts(gene, tair_code)) %>%
    pivot_longer(
      cols = where(is.numeric),
      names_to = "condition",
      values_to = "sd_value"
    ) %>%
    mutate(
      condition_group = str_extract(condition, regex$condition_group),
      time_h = extract_time(condition, regex, time_map),
      transcript_id = str_extract(gene, regex$transcript_suffix)
    ) %>%
    filter(!is.na(time_h) & !is.na(transcript_id))

  print(error_bar_df)

  plot_data_full <- avg_df %>%
    dplyr::filter(str_starts(gene, tair_code)) %>%
    pivot_longer(
      cols = where(is.numeric),
      names_to = "condition",
      values_to = "mean_value"
    ) %>%
    mutate(
      condition_group = str_extract(condition, regex$condition_group),
      time_h = extract_time(condition, regex, time_map),
      transcript_id = str_extract(gene, regex$transcript_suffix)
    ) %>%
    filter(!is.na(time_h) & !is.na(transcript_id))
  print(plot_data_full)

  all_variances <- plot_data_full %>%
    group_by(transcript_id) %>%
    summarise(variance = var(mean_value, na.rm = TRUE), mean_exp = mean(mean_value, na.rm = TRUE)) %>%
    dplyr::filter(!is.na(variance)) %>%
    arrange(desc(variance))

  print(all_variances, n = 100)

  plot_data_full <- plot_data_full %>%
    full_join(error_bar_df, by = c("transcript_id", "condition_group", "time_h"))

  print(plot_data_full)
  

  print(all_variances, n = 100)
  auto_cutoff <- find_elbow_point(all_variances$variance) # Sets the automatic cutoff
  cat("The automatically determined cutoff is at:", auto_cutoff, "transcripts.\n")

  plot_elbow_geometry(all_variances$variance) #plot elbow geometry for visual inspection

  num_to_select <- auto_cutoff

  # Checks for manual cutoff
  if (!is.null(manual_cutoff)) {
    num_to_select <- manual_cutoff
    cat("Using manual cutoff of:", num_to_select, "transcripts.\n")
  }
  #Selects top variable transcripts based on cutoff
  top_transcripts <- all_variances %>%
    slice_head(n = num_to_select) %>%
    pull(transcript_id)

  plot_data <- plot_data_full %>%
    dplyr::filter(transcript_id %in% top_transcripts)

  if (nrow(plot_data) == 0) {
    warning(paste("No transcript data to plot for", gene_name))
    return(ggplot() +
      theme_void() +
      labs(title = paste(gene_name, "\n(No transcript data)")))
  }

  # Generate plot differently based on dataset
  if (dataset == "timeseries") {
    gap_start_time <- 48
    display_gap_width <- 4
    compression_amount <- (96 - gap_start_time) - display_gap_width

    plot_data <- plot_data %>%
      mutate(
        x_val = ifelse(time_h > gap_start_time, time_h - compression_amount, time_h),
        period = ifelse(time_h <= gap_start_time, "pre", "post")
      )
    print(plot_data)

    original_breaks <- c(0, 24, 48, 96, 120, 147)
    compressed_breaks <- ifelse(original_breaks > gap_start_time, original_breaks - compression_amount, original_breaks)

    p <- ggplot(plot_data, aes(x = x_val, y = mean_value, colour = transcript_id)) +
      geom_line(aes(group = interaction(transcript_id, condition_group, period), linetype = condition_group), linewidth = 1) +
      geom_point(aes(shape = condition_group), fill = "white", size = 3.5, alpha = 0.8) +
      geom_vline(xintercept = c(gap_start_time, gap_start_time + display_gap_width), linetype = "dashed", colour = "grey40") +
      scale_x_continuous(name = "Time since chilling (Hrs)", breaks = compressed_breaks, labels = original_breaks) +
      geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 1.5, alpha = 0.5)
  } else if (dataset == "coolLL2") {
    plot_data <- plot_data %>%
    mutate(time_h = time_h - 12) 
    p <- ggplot(plot_data, aes(x = time_h, y = mean_value, colour = transcript_id)) +
      geom_line(aes(group = interaction(transcript_id, condition_group), linetype = condition_group), linewidth = 1) +
      geom_point(aes(shape = condition_group), fill = "white", size = 3.5, alpha = 0.8) +
      scale_x_continuous(name = "Time since chilling (Hrs)", breaks = scales::breaks_width(4)) +
      geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 1.5, alpha = 0.5)
  } else {
    p <- ggplot(plot_data, aes(x = time_h, y = mean_value, colour = transcript_id)) +
      geom_line(aes(group = interaction(transcript_id, condition_group), linetype = condition_group), linewidth = 1) +
      geom_point(aes(shape = condition_group), fill = "white", size = 3.5, alpha = 0.8) +
      scale_x_continuous(name = "Time since dawn (Hrs)", breaks = scales::breaks_width(24)) +
      geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 1.5, alpha = 0.5)
  }

  # Add shared aesthetics and labels
  p +
   geom_line(aes(colour = transcript_id, linetype = condition_group), linewidth = 1) +
    geom_point(aes(colour = transcript_id, shape = condition_group), fill = "white", size = 3.5, alpha = 0.8) +
    scale_shape_manual(name = "Condition", values = aesthetics$shapes, labels = aesthetics$labels) +
    scale_linetype_manual(name = "Condition", values = aesthetics$linetypes, labels = aesthetics$labels) +
    labs(
      title = paste("Top variable transcripts of", gene_name),
      subtitle = paste("TAIR ID:", tair_code),
      y = "Mean Expression (TPM)", colour = "Transcript ID"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.1, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, colour = "gray40"),
      legend.position = "bottom", legend.box = "vertical"
    ) +
    scale_colour_viridis(discrete = TRUE, option = "viridis")
}

# Plot Transcript Structures and a Comparison of the Top Two
plot_gene_structure_and_comparison <- function(
    gene_id, gene_name,
    gtf_data, avg_df, manual_cutoff = NULL) {
  # Pivot data to long format to calculate variance and mean for each transcript
  transcript_metrics_df <- avg_df %>%
    dplyr::filter(str_starts(gene, gene_id)) %>%
    # Pivot all numeric condition columns to a long format
    pivot_longer(
      cols = where(is.numeric),
      names_to = "condition",
      values_to = "expression"
    ) %>%
    # Group by the full transcript ID
    group_by(gene) %>%
    # Calculate both variance and mean 
    summarise(
      variance = var(expression, na.rm = TRUE),
      overall_mean = mean(expression, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Remove transcripts where variance couldn't be calculated
    dplyr::filter(!is.na(variance))

  if (nrow(transcript_metrics_df) == 0) {
    warning(paste("No transcripts with calculable variance found for", gene_name))
    return(NULL)
  }

  # Transcript Selection Using L-BowR

  # Sort all transcripts by their variance
  sorted_transcripts <- transcript_metrics_df %>%
    arrange(desc(variance))

  # Automatically find the elbow point from the sorted variance values
  num_to_select <- find_elbow_point(sorted_transcripts$variance)
  if (!is.null(manual_cutoff)) {
    num_to_select <- manual_cutoff
    cat("Using manual cutoff of:", num_to_select, "transcripts.\n")
  }

  # Filter to keep only the most variable transcripts selected by L-BowR
  expressed_transcripts_df <- sorted_transcripts %>%
    slice_head(n = num_to_select)

  expressed_transcripts_list <- expressed_transcripts_df %>% pull(gene)


  expressed_gene_annotation <- gtf_data %>%
    dplyr::filter(!is.na(.data$gene_id) & .data$gene_id == !!gene_id) %>%
    dplyr::filter(transcript_id %in% expressed_transcripts_list)

  if (nrow(expressed_gene_annotation) == 0) {
    warning(paste(
      "Expressed transcripts for", gene_name,
      "not found in GTF file."
    ))
    return(NULL)
  }

  # Separate features for plotting
  gene_exons <- expressed_gene_annotation %>% dplyr::filter(type == "exon")
  gene_cds <- expressed_gene_annotation %>% dplyr::filter(type == "CDS")
  gene_myb <- expressed_gene_annotation %>% dplyr::filter(type == "mybdomain_part")
  gene_cct <- expressed_gene_annotation %>% dplyr::filter(type == "cctdomain_part")

  # PLotting gene structure
  p1 <- ggplot(gene_exons, aes(xstart = start, xend = end, y = transcript_id)) +
    geom_intron(
      data = to_intron(gene_exons, "transcript_id"),
      aes(strand = strand)
    ) +
    geom_range(aes(fill = "Exon"), height = 0.25) +
    geom_range(data = gene_cds, aes(fill = "CDS"), height = 0.5) + #CDS plot
    geom_range(
      data = gene_myb,
      aes(fill = "MYB Domain"), height = 0.4, colour = "black", linewidth = 0.001 #MYB plot
    ) +
    geom_range(
      data = gene_cct,
      aes(fill = "CCT Domain"), height = 0.4, colour = "black", linewidth = 0.001 #CCT plot
    ) +
    #Colour scale for different features
    scale_fill_manual(
     name = "Feature",
      values = c("Exon" = "white", "CDS" = "royalblue", "MYB Domain" = "coral", "CCT Domain" = "green"),
      breaks = c("Exon", "CDS", "MYB Domain", "CCT Domain")
    ) +
    labs(subtitle = paste0("Expressed transcript structures for ", gene_name), y = "Transcript ID") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  top_two_transcripts <- expressed_transcripts_df %>%
    slice_max(order_by = overall_mean, n = 2) %>%
    pull(gene)

  # Checks how many transcripts were found
  if (length(top_two_transcripts) < 2) {
    warning(paste(
      "Fewer than two expressed transcripts found for",
      gene_name, ". Cannot create comparison plot."
    ))
    return(
      p1 <- p1 +plot_annotation(
        title = paste("Transcript Structure for", gene_name),
        subtitle = paste("TAIR ID:", gene_id),
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, colour = "gray40", size = 12)
        )
      )
    )
  }

  transcript1_id <- top_two_transcripts[1]
  transcript2_id <- top_two_transcripts[2]

  # Filter annotations for the top two transcripts
  t1_exons <- gene_exons %>% dplyr::filter(transcript_id == transcript1_id)
  t1_cds <- gene_cds %>% dplyr::filter(transcript_id == transcript1_id)
  t1_myb <- gene_myb %>% dplyr::filter(transcript_id == transcript1_id)
  t1_cct <- gene_cct %>% dplyr::filter(transcript_id == transcript1_id)

  t2_exons <- gene_exons %>% dplyr::filter(transcript_id == transcript2_id)
  t2_cds <- gene_cds %>% dplyr::filter(transcript_id == transcript2_id)
  t2_myb <- gene_myb %>% dplyr::filter(transcript_id == transcript2_id)
  t2_cct <- gene_cct %>% dplyr::filter(transcript_id == transcript2_id)

  y_axis_label <- paste(
    transcript1_id, "\n(Top)", "/",
    transcript2_id, "\n(Bottom)"
  )

  p2 <- ggplot(mapping = aes(y = y_axis_label)) +
    geom_half_range(
      data = t2_exons,
      aes(xstart = start, xend = end),
      fill = "white", height = 0.2, colour = "black"
    ) +
    geom_intron(
      data = to_intron(t2_exons, "transcript_id"),
      aes(xstart = start, xend = end, strand = strand),
      arrow.min.intron.length = 300
    ) +
    geom_half_range(
      data = t2_cds,
      aes(xstart = start, xend = end),
      fill = "grey40", height = 0.2
    ) +
    geom_half_range(
      data = t2_myb,
      aes(xstart = start, xend = end),
      fill = "coral", height = 0.19,
      colour = "black", linewidth = 0.001
    ) +
    geom_half_range(
      data = t2_cct,
      aes(xstart = start, xend = end),
      fill = "green", height = 0.19,
      colour = "black", linewidth = 0.001
    ) +
    geom_half_range(
      data = t1_exons,
      aes(xstart = start, xend = end),
      range.orientation = "top",
      fill = "white", height = 0.2, colour = "black"
    ) +
    geom_intron(
      data = to_intron(t1_exons, "transcript_id"),
      aes(xstart = start, xend = end, strand = strand),
      arrow.min.intron.length = 300
    ) +
    geom_half_range(
      data = t1_cds,
      aes(xstart = start, xend = end),
      fill = "purple",
      range.orientation = "top", height = 0.2
    ) +
    geom_half_range(
      data = t1_myb,
      aes(xstart = start, xend = end),
      fill = "coral", range.orientation = "top",
      height = 0.19, colour = "black", linewidth = 0.001
    ) +
    geom_half_range(
      data = t1_cct,
      aes(xstart = start, xend = end),
      fill = "green", range.orientation = "top",
      height = 0.19, colour = "black", linewidth = 0.001
    ) +
    labs(
      subtitle = "Comparison of two most abundant transcripts",
      x = "Genomic Coordinates", y = ""
    ) +
    theme_minimal(base_size = 12)

  combined_plot <- p1/p2 +
    plot_layout(heights = c(length(unique(gene_exons$transcript_id)), 4)) +
    plot_annotation(
      title = paste("Transcript Structures for", gene_name),
      subtitle = paste("TAIR ID:", gene_id),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, colour = "gray40", size = 12)
      )
    )
  return(combined_plot)
}

# Create a Full Dashboard for a Gene
plot_expression_dashboard <- function(
    gene_id, gene_name,
    gtf_data, raw_df, avg_df, dataset, manual_cutoff = NULL) {
  p_expr <- plot_gene_expression(
    tair_code = gene_id, gene_name = gene_name,
    raw_df = raw_df, avg_df = avg_df, dataset = dataset
  )
  p_ind_transcripts <- plot_individual_transcripts(
    tair_code = gene_id, gene_name = gene_name,
    avg_df = avg_df, raw_df = raw_df, dataset = dataset, manual_cutoff = manual_cutoff
  )
  p_structure <- plot_gene_structure_and_comparison(
    gene_id = gene_id, gene_name = gene_name,
    gtf_data = gtf_data, avg_df = avg_df, manual_cutoff = manual_cutoff
  )

  combined_plot <- if (is.null(p_structure)) {
    p_expr + p_ind_transcripts
  } else {
    (p_expr + p_ind_transcripts) / p_structure +
      plot_layout(heights = c(1, 1.2))
  }

  final_dashboard <- combined_plot + plot_annotation(
    title = paste(
      "Expression and Structure Dashboard for",
      gene_name, "-", toupper(dataset)
    ),
    theme = theme(plot.title = element_text(
      hjust = 0.5,
      face = "bold", size = 20
    ))
  )

  output_filename <- file.path(figures_dir, paste0(
    "Dashboard_", gene_name,
    "_", toupper(dataset), ".svg"
  ))
  ggsave(
    plot = final_dashboard, filename = output_filename,
    width = 12, height = 12
  )
  cat(paste("Dashboard saved as", output_filename, "\n"))

  return(final_dashboard)
}

# Plot a Grid of Total Gene Expression for a List of Genes
# Generates a faceted grid plot showing total gene expression for each gene in the provided list.
plot_total_expression_grid <- function(raw_df, avg_df, gene_list, dataset) {
  if (!("gene" %in% colnames(raw_df))) {
    raw_df <- tibble::rownames_to_column(raw_df, var = "gene")
  }

  regex <- get_regex_patterns(dataset)
  aesthetics <- get_plot_aesthetics(dataset)
  time_map <- if (regex$pattern_type == "mapped") get_time_mapping() else NULL

  all_tair_codes <- gene_list$tair_code

  # Prepare data for plotting
  plot_data <- raw_df %>%
    mutate(tair_code = str_extract(gene, regex$tair_code)) %>%
    filter(tair_code %in% all_tair_codes) %>%
    pivot_longer(
      cols = where(is.numeric),
      names_to = "sample", values_to = "abundance"
    ) %>%
    mutate(
      condition_group = str_extract(sample, regex$condition_group),
      time_h = extract_time(sample, regex, time_map)
    ) %>%
    group_by(tair_code, sample, condition_group, time_h) %>%
    summarise(
      total_abundance_per_rep = sum(abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(tair_code, condition_group, time_h) %>%
    summarise(
      mean_abundance = mean(total_abundance_per_rep),
      sd_abundance = sd(total_abundance_per_rep), .groups = "drop"
    ) %>%
    left_join(gene_list, by = "tair_code") %>%
    filter(!is.na(time_h) & !is.na(condition_group))

  if (nrow(plot_data) == 0) {
    warning("No data available for plotting in the specified gene list.")
    return(NULL)
  }

  # Generate plot differently based on dataset
  if (dataset == "timeseries") {
    gap_start_time <- 48
    display_gap_width <- 4
    compression_amount <- (96 - gap_start_time) - display_gap_width
    plot_data <- plot_data %>%
      mutate(
        x_val = ifelse(time_h > gap_start_time, time_h - compression_amount, time_h),
        period = ifelse(time_h <= gap_start_time, "pre", "post")
      )
    original_breaks <- c(0, 48, 96, 120)
    compressed_breaks <- ifelse(original_breaks > gap_start_time, original_breaks - compression_amount, original_breaks)
    scale_layer <- scale_x_continuous(name = "Time (Hours)", breaks = compressed_breaks, labels = original_breaks)
    vline_layer <- geom_vline(xintercept = c(gap_start_time, gap_start_time + display_gap_width), linetype = "dashed", colour = "grey40")

    p <- ggplot(plot_data, aes(x = x_val, y = mean_abundance, colour = condition_group)) +
      geom_errorbar(aes(ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), width = 1.5, alpha = 0.5) +
      geom_line(aes(group = interaction(condition_group, period), linetype = condition_group), linewidth = 1) +
      geom_point(aes(shape = condition_group), fill = "white", size = 2.5) +
      facet_wrap(~gene_name, scales = "free_y", ncol = 3) +
      vline_layer +
      scale_layer
  } else {
    p <- ggplot(plot_data, aes(x = time_h, y = mean_abundance, colour = condition_group)) +
      geom_errorbar(aes(ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), width = 1.5, alpha = 0.5) +
      geom_line(aes(linetype = condition_group), linewidth = 1) +
      geom_point(aes(shape = condition_group), fill = "white", size = 2.5) +
      facet_wrap(~gene_name, scales = "free_y", ncol = 3) +
      scale_x_continuous(name = "Time (Hours)", breaks = scales::breaks_width(24))
  }

  # Add shared aesthetics and labels
  p +
    scale_colour_manual(name = "Condition", values = aesthetics$values, labels = aesthetics$labels) +
    scale_shape_manual(name = "Condition", values = aesthetics$shapes, labels = aesthetics$labels) +
    scale_linetype_manual(name = "Condition", values = aesthetics$linetypes, labels = aesthetics$labels) +
    labs(
      title = paste("Total Gene Expression Across All Genes -", toupper(dataset)),
      y = "Mean Expression (TPM)"
      
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 10)
    )
}


# Plot diversity comparison of relevant gene
plot_diversity_comparison <- function(
    tair_code, gene_name,
     raw_df, dataset,time_h= NULL,
    day_to_compare = NULL) {
  #  Input Validation 
  if (any(raw_df != round(raw_df), na.rm = TRUE)) {
    raw_df <- round(raw_df)
  }

  # Setup 
  regex <- get_regex_patterns(dataset)
  aesthetics <- get_plot_aesthetics(dataset)
  time_map <- if (regex$pattern_type == "mapped") get_time_mapping() else NULL

  # Data Preparation 
  # Create a base long-format data frame from raw data.
  base_long_df <- raw_df %>%
    tibble::rownames_to_column(var = "transcript_id") %>%
    dplyr::filter(str_starts(transcript_id, tair_code)) %>%
    pivot_longer(
      cols = -transcript_id,
      names_to = "sample",
      values_to = "abundance"
    ) %>%
    dplyr::rowwise() %>%
    mutate(
      condition_group = str_extract(sample, regex$condition_group),
      time_h_calc = extract_time(sample, regex, time_map)
    ) %>%
    ungroup()

  #  Conditional Logic for Data Filtering 
  if (dataset == "coolLL2") {
    if (is.null(day_to_compare) || !day_to_compare %in% c(1, 2, 3)) {
      stop("For 'coolLL2', 'day_to_compare' must be 1, 2, or 3 to define the start of the window.", call. = FALSE)
    }

    # Define the 36-hour rolling window starting from the beginning of the specified day
    start_hour <- (day_to_compare - 1) * 24
    potential_end_hour <- start_hour + 36
    max_data_hour <- max(base_long_df$time_h_calc, na.rm = TRUE)

    # Ensure the window does not exceed the available data
    end_hour <- min(potential_end_hour, max_data_hour)

    cat(paste("\ncoolLL2 Dataset: Finding Peak Abundance starting on Day ", day_to_compare, " (Using actual window: ", start_hour, "-", end_hour, "h) \n"))

    df_window <- base_long_df %>%
      filter(time_h_calc >= start_hour, time_h_calc <= end_hour)

    # Peak Matching Logic 
    cat("coolLL2: Finding temporally close high-abundance peaks\n")

    # Calculate total abundance for all time points in the window
    abundance_in_window <- df_window %>%
      group_by(condition_group, time_h_calc) %>%
      summarise(total_abundance_at_time = sum(abundance, na.rm = TRUE), .groups = "drop")

    # Find the absolute highest peak across both conditions
    anchor_peak <- abundance_in_window %>%
      slice_max(order_by = total_abundance_at_time, n = 1, with_ties = FALSE)

    anchor_condition <- anchor_peak$condition_group
    anchor_time <- anchor_peak$time_h_calc

    cat(paste("Anchor peak found in condition '", anchor_condition, "' at", anchor_time, "h\n"))

    # Identify the other condition.
    all_conditions <- unique(abundance_in_window$condition_group)
    search_condition <- setdiff(all_conditions, anchor_condition)

    # Define a search window around the anchor time to find the corresponding peak
    search_window_start <- max(start_hour, anchor_time - 8)
    search_window_end <- min(end_hour, anchor_time + 12)

    cat(paste("Searching for corresponding peak for '", search_condition, "' in window [", search_window_start, "h,", search_window_end, "h]\n"))

    # Find the highest peak for the search condition within the smaller window
    search_peak <- abundance_in_window %>%
      filter(
        condition_group == search_condition,
        time_h_calc >= search_window_start,
        time_h_calc <= search_window_end
      ) %>%
      slice_max(order_by = total_abundance_at_time, n = 1, with_ties = FALSE)

    search_time <- search_peak$time_h_calc
    # Combine the two selected peaks into the final peak_times data frame
    peak_times <- tibble(
      condition_group = c(anchor_condition, search_condition),
      peak_time_h = c(anchor_time, search_time)
    )

    print(peak_times)

    long_format_df <- df_window %>%
      inner_join(peak_times, by = "condition_group") %>%
      filter(time_h_calc == peak_time_h)

    peak_time_1 <- peak_times$peak_time_h[peak_times$condition_group == aesthetics$labels[1]]
    peak_time_2 <- peak_times$peak_time_h[peak_times$condition_group == aesthetics$labels[2]]

    plot_title <- paste("Shannon Diversity of", gene_name, "at Peak Abundance")
    plot_subtitle_extra <- paste0(aesthetics$labels[1], " at ", peak_times[[1, 2]], "h vs. ", aesthetics$labels[2], " at ", peak_times[[2, 2]], "h")
  } else if (dataset == "timeseries_df") {
    # ANCHOR POINT LOGIC FOR TIMESERIES_DF 
    time_window <- 24
    start_hour <- 0
    end_hour <- time_window
    cat(paste0("\n timeseries_df Dataset: Finding Peak Abundance using Anchor Logic (", start_hour, "-", end_hour, "h window) \n"))

    df_window <- base_long_df %>% filter(time_h_calc >= start_hour, time_h_calc <= end_hour)

    # Calculate total abundance across the window
    abundance_in_window <- df_window %>%
      group_by(condition_group, time_h_calc) %>%
      summarise(total_abundance_at_time = sum(abundance, na.rm = TRUE), .groups = "drop")

    # Find the absolute highest peak across both conditions to serve as the anchor
    anchor_peak <- abundance_in_window %>%
      slice_max(order_by = total_abundance_at_time, n = 1, with_ties = FALSE)

    anchor_condition <- anchor_peak$condition_group
    anchor_time <- anchor_peak$time_h_calc

    cat(paste("Anchor peak found in condition '", anchor_condition, "' at", anchor_time, "h\n"))

    # Identify the other condition to search within
    all_conditions <- unique(abundance_in_window$condition_group)
    search_condition <- setdiff(all_conditions, anchor_condition)

    # Define a search window around the anchor time
    search_window_start <- max(start_hour, anchor_time - 8)
    search_window_end <- min(end_hour, anchor_time + 8)

    cat(paste("Searching for corresponding peak for '", search_condition, "' in window [", search_window_start, "h,", search_window_end, "h]\n"))

    # Find the highest peak for the search condition within the search window
    search_peak <- abundance_in_window %>%
      filter(
        condition_group == search_condition,
        time_h_calc >= search_window_start,
        time_h_calc <= search_window_end
      ) %>%
      slice_max(order_by = total_abundance_at_time, n = 1, with_ties = FALSE)

    # Fallback: If no peak is found in the narrow search window, take the highest peak from the entire window
    if (nrow(search_peak) == 0) {
      warning(paste("No corresponding peak found for '", search_condition, "' within the search window. Using the highest peak from the entire 24h window instead."), call. = FALSE)
      search_peak <- abundance_in_window %>%
        filter(condition_group == search_condition) %>%
        slice_max(order_by = total_abundance_at_time, n = 1, with_ties = FALSE)
    }

    search_time <- search_peak$time_h_calc

    # Combine the anchor and search peaks
    peak_times <- tibble(
      condition_group = c(anchor_condition, search_condition),
      peak_time_h = c(anchor_time, search_time)
    )

    print(peak_times)

    long_format_df <- df_window %>%
      inner_join(peak_times, by = "condition_group") %>%
      filter(time_h_calc == peak_time_h)

    peak_time_1 <- peak_times$peak_time_h[peak_times$condition_group == names(aesthetics$labels)[1]]
    peak_time_2 <- peak_times$peak_time_h[peak_times$condition_group == names(aesthetics$labels)[2]]

    plot_title <- paste("Shannon Diversity of", gene_name, "at Peak Abundance")
    plot_subtitle_extra <- paste0(aesthetics$labels[1], " at ", peak_times[[1, 2]], "h vs. ", aesthetics$labels[2], " at ", peak_times[[2, 2]], "h")
  } else {
    #  Logic for 24h Dataset: Find Single Anchor Peak 
    # This logic finds the time of the absolute highest abundance across all conditions and uses that single time point for comparing diversity.
    cat("\n 24h Analysis: Finding single peak time (anchor) for comparison \n")

    # Calculate total abundance at each time point for each condition.
    abundance_by_time <- base_long_df %>%
      group_by(condition_group, time_h_calc) %>%
      summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

    # Find the single time point with the absolute highest abundance across all conditions.
    # If there's a tie, slice_max takes the first one.
    anchor_peak <- abundance_by_time %>%
      slice_max(order_by = total_abundance, n = 1, with_ties = FALSE)

    # Extract the single time point to use for comparison.
    peak_time_h <- anchor_peak$time_h_calc

    cat(paste("Overall peak abundance found at ", peak_time_h, "h (in condition '", anchor_peak$condition_group, "').\n"))
    cat(paste("Using ", peak_time_h, "h for comparison across all conditions.\n"))

    # Filter the original long-format data to only include rows at this single peak time.
    long_format_df <- base_long_df %>%
      dplyr::filter(time_h_calc == .env$peak_time_h)

    # Set plot titles to reflect the comparison at a single, data-driven time point.
    plot_title <- paste("Shannon Diversity of", gene_name, "at Peak Abundance")
    plot_subtitle_extra <- paste("Comparison at", peak_time_h, "h")
    # Ensure peak_times is null as this case compares at a single time point.
    peak_times <- NULL
  }

  #  Pooled Data Calculation 
  pooled_data <- long_format_df %>%
    group_by(condition_group, transcript_id) %>%
    summarise(pooled_abundance = sum(abundance, na.rm = TRUE), .groups = "drop")



  #  Main Diversity Calculation 
  diversity_summary_df <- pooled_data %>%
    filter(pooled_abundance > 0) %>%
    group_by(condition_group) %>%
    mutate(N = sum(pooled_abundance), P = pooled_abundance / N) %>%
    mutate(P_lnP = P * log(P), P_lnP_sq = P * (log(P))^2) %>%
    summarise(
      Richness = n(), N = mean(N, na.rm = TRUE), H = -sum(P_lnP),
      sum_P_lnP_sq = sum(P_lnP_sq), .groups = "drop"
    ) %>%
    mutate(Var_H = ((sum_P_lnP_sq - (H^2)) / N) + ((Richness - 1) / (2 * N^2))) %>%
    dplyr::select(condition_group, Richness, N, H, Var_H)
  print(" Main Diversity Summary ")
  print(diversity_summary_df)

  # Statistical Analysis 
  # Initialise p-value to NA to prevent errors if it is not calculated
  p_value_ttest <- NA

  if (nrow(diversity_summary_df) == 2) {
    #  Hutchinson's t-test 
    cat("\n=== Hutchinson's t-test ===\n")
    group1_H <- diversity_summary_df$H[1]
    group1_Var_H <- diversity_summary_df$Var_H[1]
    group1_N <- diversity_summary_df$N[1]
    group2_H <- diversity_summary_df$H[2]
    group2_Var_H <- diversity_summary_df$Var_H[2]
    group2_N <- diversity_summary_df$N[2]
    t_statistic <- (group1_H - group2_H) / sqrt(group1_Var_H + group2_Var_H)
    numerator_df <- (group1_Var_H + group2_Var_H)^2
    denominator_df <- ((group1_Var_H^2) / group1_N) + ((group2_Var_H^2) / group2_N)
    degrees_of_freedom <- numerator_df / denominator_df
    p_value_ttest <- 2 * pt(-abs(t_statistic), df = degrees_of_freedom)
    cat("p-value (t-test):", p_value_ttest, "\n\n")
  } else {
    warning("Skipping t-test because the number of groups is not equal to 2.")
  }

  # Plotting 
  if (!is.na(p_value_ttest)) {
    # Bar plot of Shannon Diversity 
    plot_data <- diversity_summary_df %>%
      mutate(CI_width = sqrt(Var_H) * 1.96, ymin = H - CI_width, ymax = H + CI_width)
    p_label_ttest <- ifelse(p_value_ttest < 0.001, "p < 0.001", paste("p =", round(p_value_ttest, 3)))
    subtitle_text <- paste(plot_subtitle_extra, paste("Hutchinson's t-test:", p_label_ttest), sep = "\n")
    y_position <- max(plot_data$ymax, na.rm = TRUE) * 1.15
    significance_bracket <- tibble(x = c(1, 1, 2, 2), y = c(y_position, y_position * 1.02, y_position * 1.02, y_position))

    p_bar <- ggplot(plot_data, aes(x = condition_group, y = H, fill = condition_group)) +
      geom_bar(stat = "identity", colour = "black", alpha = 0.6) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, linewidth = 0.7) +
      geom_line(data = significance_bracket, aes(x = x, y = y), inherit.aes = FALSE) +
      annotate("text", x = 1.5, y = y_position * 1.05, label = "p-value", size = 4.5) +
      scale_x_discrete(labels = aesthetics$labels) +
      scale_fill_manual(values = aesthetics$values, guide = "none") +
      labs(title = plot_title, subtitle = subtitle_text, x = "Condition", y = "Shannon Index (H)") +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 11, lineheight = 1.2))


    combined_plot <- p_bar

    # Save the plot 
    if (!dir.exists(figures_dir)) {
      dir.create(figures_dir, recursive = TRUE)
    }
    file_suffix <- if (dataset == "coolLL2") {
      paste0("PeakHours_Day", day_to_compare, "_36hWindow")
    } else if (dataset == "timeseries_df") {
      "PeakHours"
    } else {
      as.character(time_h)
    }
    output_filename <- file.path(figures_dir, paste0("SdComparison_", gene_name, "_", file_suffix, "_", toupper(dataset), ".svg"))
    ggsave(plot = combined_plot, filename = output_filename, width = 8, height = 10) 
    cat(paste("SdComparison saved as", output_filename, "\n"))
  }
  return(combined_plot)
}


# Finds the elbow point.
find_elbow_point <- function(values) {
  # Get the coordinates of the points on the curve
  n_points <- length(values)
  all_coords <- data.frame(x = 1:n_points, y = values)
  # Normalise the coordinates to a 0-1 scale.
  normalised_coords <- data.frame(
    x_norm = (all_coords$x - min(all_coords$x)) / (max(all_coords$x) - min(all_coords$x)),
    y_norm = (all_coords$y - min(all_coords$y)) / (max(all_coords$y) - min(all_coords$y))
  )
  # Define the line between the first and last points in normalised coordinates
  # The line is from P1(0, 1) to P2(1, 0) for a descending curve.
  # The equation of this line is x + y - 1 = 0.
  # Calculate the distance from each point (x_n, y_n) to the line.
  distances <- abs(normalised_coords$x_norm + normalised_coords$y_norm - 1) / sqrt(2)
  # The elbow point is the one with the maximum distance
  elbow_index <- which.max(distances)

  return(elbow_index)
}

plot_elbow_geometry <- function(values) {
  n_points <- length(values)
  x_vals <- 1:n_points
  # Normalise
  x_norm <- (x_vals - min(x_vals)) / (max(x_vals) - min(x_vals))
  y_norm <- (values - min(values)) / (max(values) - min(values))
  # Distance to line x + y = 1
  distances <- abs(x_norm + y_norm - 1) / sqrt(2)
  elbow_index <- which.max(distances)
  # Base plot
  plot(x_norm, y_norm,
    type = "b", main = "Elbow Detection Visualised",
    xlab = "Normalised Index", ylab = "Normalised Variance", pch = 16
  )
  # Reference line from (0,1) to (1,0)
  abline(a = 1, b = -1, col = "gray", lty = 2)

  for (i in 1:n_points) {
    # Original point
    x0 <- x_norm[i]
    y0 <- y_norm[i]
    # Calculate the x-coordinate of the projected point (xp)
    xp <- (x0 - y0 + 1) / 2
    # Calculate the y-coordinate of the projected point (yp)
    yp <- -xp + 1
    # Draw the original point
    points(x0, y0, pch = 19, col = "blue") # pch=19 is a solid circle
    # Draw the segment from the original point to its projection
    segments(x0, y0, xp, yp, col = "skyblue", lty = 3, lwd = 1.5)
    # Draw the projected point on the line
    points(xp, yp, pch = 19, col = "darkgreen")
  }
  # Highlight the elbow point
  points(x_norm[elbow_index], y_norm[elbow_index], col = "red", pch = 19, cex = 1.5)
  text(x_norm[elbow_index], y_norm[elbow_index], labels = "Elbow", pos = 3, col = "red")
}
