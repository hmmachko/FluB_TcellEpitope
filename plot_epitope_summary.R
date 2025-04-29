## Clear Global Environment
rm(list = ls())


#### Session prep ####
## Install packages and load libraries as required
if(!require(tidyverse)){
  install.packages("tidyverse",dependencies = T)
  library(tidyverse)
}
if(!require(dplyr)){
  install.packages("dplyr",dependencies = T)
  library(dplyr)
}

if(!require(ggplot2)){
  install.packages("ggplot2",dependencies = T)
  library(ggplot2)
}

if(!require(purrr)){
  install.packages("purrr")
  library(purrr)
}


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")

#### Directories ####
dir <-      paste("/Volumes/PrimateFS/TCF lab/Current Lab Members/machkovech_heather/bioinformatics/13197/processed_data/recall_fn_ann/", sep="")
dir_processed <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/machkovech_heather/bioinformatics/13197/processed_data/", sep="")
dir_save <-paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/machkovech_heather/bioinformatics/13197/figures/", sep="")
dir_resource <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/machkovech_heather/bioinformatics/13197/resources/", sep="")


## theme
axis_formatting <- theme(axis.text.x = element_text(size = 8),
                         axis.text.y = element_text(size = 8),
                         axis.title.x = element_text(size = 8, margin = margin(t = 4)),
                         axis.title.y = element_text(size = 8, margin = margin(r = 4)))

legend_formatting <- theme(legend.text = element_text(size = 8),
                           legend.key.height= unit(0.5, 'cm'),
                           legend.key.width= unit(0.5, 'cm'))

background_formatting <- theme(panel.border = element_rect(color = "grey", fill = NA, size = .5),
                               panel.grid = element_blank(),
                               strip.background = element_blank(),
                               panel.background = element_blank(),
                               legend.background = element_blank())

## plot_grid
Size_adjust = 12
LR_adjust = -0.5 # less = right
UD_adjust = 1.1 # less = down 

## palettes

# mutation types without stop
palette_muts_NS_S <- c("nonsynonymous" = "#FF7F20", "synonymous" = "#4F7899")


###plot mutations in epitope regions

#### read in variants table and wrangle to require variants in both replicates, filter to be at least 3% in both replicates

##read in processed file
final_result_filt= "processed_variants_table.csv"
final_result_filtered <- read.csv(file = file.path(dir_processed, final_result_filt))

# Transform data to wider format ensuring correct matching
final_result_wide <- final_result_filtered %>%
  mutate(
    sample = str_remove(Sample, "-R[12]$"),  # Remove -R1 or -R2 suffix
    Replicate = str_extract(Sample, "R[12]$")    # Extract R1 or R2
  ) %>%
  select(sample, segment, gene, POS, REF, ALT, substitution_type, hgvs_c, hgvs_p, Replicate, AF) %>%  # Keep relevant columns
  pivot_wider(
    names_from = Replicate, 
    values_from = AF,
    names_prefix = "AF_",
    values_fill = NA  # Fill missing values with NA
  ) %>% rename_with(tolower)

variants_combined_filter <- final_result_wide  %>%
  filter(af_r1 >= 0.03 & af_r2 >= 0.03)


final_result_f<- variants_combined_filter %>%
  as.data.frame() %>%
  # More granular priority system
  mutate(priority = case_when(
    #substitution_type == "nonsynonymous" & gene == "PA" ~ 1,  # Top priority: nonsynonymous PA
    substitution_type == "nonsynonymous" ~ 1,                 # Second: other nonsynonymous
    #gene == "PA" ~ 3,                                         # Third: synonymous PA
    TRUE ~ 2                                                 # Last: everything else
  )) %>%
  group_by(sample, segment, pos, ref, alt) %>%
  arrange(priority, .by_group = TRUE) %>%
  distinct(sample, segment, pos, ref, alt, .keep_all = TRUE) %>%
  select(-priority) %>%
  ungroup()

final_result_f  <- final_result_f  %>%
  mutate(
    segment = ifelse(is.na(segment) | segment == "NA", "neuraminidase", segment))

##############

final_result <- final_result_f 

# Convert 3-letter to 1-letter amino acid codes
# Create a lookup table for amino acid abbreviations
aa_map <- c("Ala" = "A", "Cys" = "C", "Asp" = "D", "Glu" = "E", "Phe" = "F",
            "Gly" = "G", "His" = "H", "Ile" = "I", "Lys" = "K", "Leu" = "L",
            "Met" = "M", "Asn" = "N", "Pro" = "P", "Gln" = "Q", "Arg" = "R",
            "Ser" = "S", "Thr" = "T", "Val" = "V", "Trp" = "W", "Tyr" = "Y")

# Function to modify hgvs_p entries
modify_aa <- function(hgvs_p) {
  # Remove the 'p.' prefix
  hgvs_p <- sub("^p\\.", "", hgvs_p)
  
  # Split the entry into parts
  matches <- regmatches(hgvs_p, gregexpr("[A-Za-z]+|\\d+", hgvs_p))
  
  # Ensure we have valid matches (3 parts: original AA, position, new AA)
  if (length(matches[[1]]) == 3) {
    aa_start <- matches[[1]][1]
    position <- matches[[1]][2]
    aa_end <- matches[[1]][3]
    
    # Convert to one-letter codes
    aa_start_one <- aa_map[aa_start]
    aa_end_one <- aa_map[aa_end]
    
    # Create the new substitution string
    if (!is.na(aa_start_one) && !is.na(aa_end_one)) {
      return(paste0(aa_start_one, position, aa_end_one))
    }
  }
  
  # Return NA if conversion fails
  return(NA)
}

# Apply the function to the hgvs_p column
final_result$aa_sub <- sapply(final_result$hgvs_p, modify_aa)

####################################

## extract epitope regions, 

# Load GTF file
gtf_path <- "B_Iowa_06_2017_VICTORIA.gtf"
setwd(dir_resource)
# Read the GTF file
gtf <- read_delim(gtf_path, delim = "\t", comment = "#", col_names = FALSE)

# Assign column names based on standard GTF format
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Filter only CDS regions
cds_regions <- gtf %>% filter(feature == "CDS")

# Extract gene names from the "attributes" column (assuming format: gene_id "GENE_NAME";)
cds_regions$gene <- str_extract(cds_regions$attributes, 'gene_id "[^"]+"') %>%
  str_remove_all('gene_id "|;"')  # Clean extracted names

# Define the genes of interest
genes_of_interest <- c("PB1",  "HA", "NP")

# Remove any double quotes (") and extra spaces from gene names
cds_regions$gene <- str_remove_all(cds_regions$gene, '"')

# Trim any accidental leading/trailing spaces
cds_regions$gene <- str_trim(cds_regions$gene)

fasta_path <- file.path(dir_resource, "B_Iowa_06_2017_VICTORIA_ref.fasta")
# Load reference FASTA

ref_fasta <- Biostrings::readDNAStringSet(fasta_path)

# Function to extract CDS sequences based on GTF coordinates
extract_cds_sequence <- function(seqname, start, end, strand, fasta) {
  if (!(seqname %in% names(fasta))) {
    warning(paste("Sequence", seqname, "not found in FASTA. Skipping."))
    return(NA)  # Return NA instead of error
  }
  
  seq <- fasta[[seqname]]  # Get the corresponding segment
  cds_seq <- Biostrings::subseq(seq, start = start, end = end)  # Extract CDS region
  
  # Reverse complement if on the negative strand
  if (strand == "-") {
    cds_seq <- reverseComplement(cds_seq)
  }
  
  return(as.character(cds_seq))
}
# Apply function to extract sequences for each row in cds_filtered
cds_regions$cds_seq <- pmap_chr(
  cds_regions[, c("seqname", "start", "end", "strand")], 
  extract_cds_sequence, fasta = ref_fasta
)

# Filter the CDS regions dataset to include only the selected genes
cds_filter <- cds_regions %>%
  filter(seqname %in% genes_of_interest)

# Translate CDS sequences to protein
cds_filter$protein_seq <- sapply(cds_filter$cds_seq, function(seq) {
  as.character(Biostrings::translate(Biostrings::DNAString(seq)))
})

# Define epitope regions
epitope_regions <- data.frame(
  gene = c("NP", "PB1", "HA", "PB1"),
  epitope_start = c(92, 413, 305,332),
  epitope_end = c(99, 421, 313,343),
  epitope_start_display = c(92, 413, 307,332),
  epitope_end_display  = c(99, 421, 315,343),
  stringsAsFactors = FALSE
)

# Initialize an empty list to store results
epitope_data <- list()

# Loop through each epitope region
for (i in seq_len(nrow(epitope_regions))) {
  gene_name <- epitope_regions$gene[i]
  start_pos <- epitope_regions$epitope_start[i]
  end_pos <- epitope_regions$epitope_end[i]
  
  # Find the corresponding gene in cds_regions
  gene_row <- cds_filter[cds_filter$seqname == gene_name, ]
  
  if (nrow(gene_row) > 0) {
    # Get the translated protein sequence
    protein_seq <- as.character(translate(DNAString(gene_row$cds_seq)))
    
    # Extract epitope sequence
    epitope_seq <- substr(protein_seq, start_pos, end_pos)
    
    # Extract pre-epitope and post-epitope sequences (handling boundaries)
    pre_epitope_start <- max(1, start_pos - 5)
    pre_epitope_seq <- substr(protein_seq, pre_epitope_start, start_pos - 1)
    
    post_epitope_end <- min(nchar(protein_seq), end_pos + 5)
    post_epitope_seq <- substr(protein_seq, end_pos + 1, post_epitope_end)
    
    # Store results in the list
    epitope_data[[i]] <- data.frame(
      gene = gene_name,
      epitope_start = start_pos,
      epitope_end = end_pos,
      epitope_seq = epitope_seq,
      pre_epitope_seq = pre_epitope_seq,
      post_epitope_seq = post_epitope_seq,
      stringsAsFactors = FALSE
    )
  }
}

# Combine list into a data frame
epitope_df <- do.call(rbind, epitope_data)

# Adding an order column, numbering from 1 to n
epitope_df <- epitope_df %>%
  mutate(epitope_start_display = c(92, 413, 307,332),
         epitope_end_display  = c(99, 421, 315,343))

# Create an empty list to store plots
epitope_plots <- list()

#updated loop 

for (i in seq_len(nrow(epitope_df))) {
  plot_num <- i
  gene_name <- epitope_df$gene[i]
  
  # Extract original numbering
  start_pos <- epitope_df$epitope_start[i]
  end_pos <- epitope_df$epitope_end[i]
  epitope_start_display <- epitope_df$epitope_start_display[i]
  epitope_end_display <- epitope_df$epitope_end_display[i]
  
  # Define extended region (extracted in YOUR numbering schema)
  plot_start <- start_pos - 5
  plot_end <- end_pos + 5
  
  # Extract the full sequence for plotting
  full_seq <- paste0(epitope_df$pre_epitope_seq[i], epitope_df$epitope_seq[i], epitope_df$post_epitope_seq[i])
  seq_positions <- seq(plot_start, plot_end)  # Positions in your schema
  
  # Compute display numbering (for collaborator's schema)
  display_positions <- seq(epitope_start_display - 5, epitope_end_display + 5)
  
  # Create a data frame for sequence text representation
  seq_df <- data.frame(
    position = seq_positions,  # YOUR schema
    display_position = display_positions,  # Collaborator's schema (for x-axis)
    amino_acid = strsplit(full_seq, "")[[1]],
    color = c(rep("gray", 5), rep("black", end_pos - start_pos + 1), rep("gray", 5))
  )
  
  mutation_df <- final_result %>%
    filter(segment == gene_name) %>%
    mutate(
      aa_position = as.numeric(str_extract(aa_sub, "[0-9]+")),  # Original numbering (collaborator's)
      
      # **Corrected shift calculation**
      aa_position_adjusted = aa_position + (epitope_start_display - start_pos),
      
      mutation_type = ifelse(substitution_type == "nonsynonymous", "nonsynonymous", "synonymous"),
      color_intensity = pmax(0.03, af_r1),
      color_value = case_when(
        substitution_type == "nonsynonymous" ~ scales::col_numeric(
          c("white", palette_muts_NS_S["nonsynonymous"]), c(0.03, 1)
        )(color_intensity),
        substitution_type == "synonymous" ~ scales::col_numeric(
          c("white", palette_muts_NS_S["synonymous"]), c(0.03, 1)
        )(color_intensity),
        TRUE ~ "gray"
      )
    ) %>%
    # Ensure mutations aren't filtered incorrectly
    filter(aa_position_adjusted >= (epitope_start_display - 5) & aa_position_adjusted<= (epitope_end_display + 5))
  
  
  
  # Generate the plot with the adjusted numbering
  p <- ggplot() +
    # Title uses adjusted numbering
    ggtitle(bquote(.(gene_name)[.(epitope_start_display) - .(epitope_end_display)])) +
    
    # Add mutations with correctly adjusted positions
    geom_point(data = mutation_df, aes(
      x = aa_position_adjusted, y = 0.08,  
      size = af_r1, 
      fill = color_value
    ), shape = 21, color = "black", alpha = 1) +
    
    # Add sequence text using collaborator's numbering
    geom_text(data = seq_df, aes(x = display_position, y = 0.03, label = amino_acid, color = color), 
              size = 3.5, family = "mono") +  
    
    # Ensure x-axis labels use collaborator's schema
    #scale_x_continuous(
    #breaks = seq_df$display_position, 
    #labels = seq_df$display_position
    # ) +
    
    # Adjust x-axis to label every 5th position
    scale_x_continuous(
      breaks = seq(min(seq_df$display_position), max(seq_df$display_position), by = 5), 
      labels = seq(min(seq_df$display_position), max(seq_df$display_position), by = 5)
    ) +
    
    # Set color scale for text  
    scale_color_manual(values = c("gray" = "gray", "black" = "black"), guide = "none") +
    scale_fill_identity() +  
    scale_size(range = c(2, 5), guide = "none") +  # No size legend
    
    # Very condensed y-axis
    ylim(0, 0.13) +
    
    
    
    # Formatting for compact display
    theme_minimal(base_family = "mono") +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, margin = margin(b = 2)),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.text.x = element_text(size = 8)
    ) +
    
    
    # X-axis label
    labs(x = NULL)
  
  # Store plot in list
  epitope_plots[[i]] <- p
}


#legend
data_syn <- data.frame(x = seq(0, 1, length.out = 100), y = 1)
syn_legend <- ggplot(data.frame(x = seq(0, 1, length.out = 100), y = 1)) +
  geom_tile(aes(x = x, y = y, fill = x)) +
  scale_fill_gradient(
    low = "white", 
    high = "#4F7899",
    name = "synonymous",
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.3, "cm")
  )
nonsyn_legend <- ggplot(data.frame(x = seq(0, 1, length.out = 100), y = 1)) +
  geom_tile(aes(x = x, y = y, fill = x)) +
  scale_fill_gradient(
    low = "white", 
    high = "#FF7F20",
    name = "nonsynonymous",
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.3, "cm")
  )
syn_components <- cowplot::get_plot_component(syn_legend, "guide-box", return_all = TRUE)
print(syn_components)
syn_grob <- syn_components[[3]]$grobs[[1]] 
nonsyn_components <- cowplot::get_plot_component(nonsyn_legend, "guide-box", return_all = TRUE)
print(nonsyn_components)
nonsyn_grob <- nonsyn_components[[3]]$grobs[[1]] 

legend_combined <- plot_grid(
  ggplot() + 
    annotate("text", x = 1, y = 1, label = "substitution type", size = 2.5) +
    theme_void(),
  plot_grid(syn_grob, nonsyn_grob, ncol = 1),
  ncol = 1,
  rel_heights = c(0.2, 0.8) #,

)


x_axis_title_plot <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Amino Acid Position", size = 3.5, fontface = "bold") +
  theme_void()  

final_plot <- wrap_plots(
  epitope_plots[[1]], epitope_plots[[2]], epitope_plots[[3]], epitope_plots[[4]], 
  nrow = 1
) /
  x_axis_title_plot / 
  legend_combined +
  plot_layout(nrow = 3, heights = c(2, 0.2, 2))  # You can tweak heights here

setwd(dir_save)
ggsave("epitope_isnv_summary_R1.png", final_plot , width = 7.08, height = 2.8, units = "in", dpi = 320)















