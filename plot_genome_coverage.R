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


############################################
### plot genome wide sequencing coverage for genes containing epitope regions with epitopes highlighted

## extract CDS and epitope regions

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


##########################################

###### READ IN COVERAGE 
coverage_file = "1920B007-R1_merged_position_coverage.tsv"
coverage <- read.csv(file = file.path(dir_processed, coverage_file))

# Read the TSV file (no header)
data <- read.table(file = file.path(dir_processed, coverage_file), header = FALSE, col.names = c("gene", "position", "depth"))

# Define facet order
gene_order <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS")

# Convert gene column to factor with the desired order
data$gene <- factor(data$gene, levels = gene_order)

cds_filter_df <- as.data.frame(cds_filter)


# Extract coding start sites for NP, PB1, and HA from cds_filter
coding_start_sites <- cds_filter %>%
  filter(seqname %in% c("NP", "PB1", "HA")) %>%
  select(seqname, start) 

# Define epitope amino acid positions
epitope_regions <- data.frame(
  seqname = c("NP", "PB1", "HA", "PB1"),
  aa_start = c(92, 413, 307,332),
  aa_end = c(99, 421, 315,343)
)

# Merge with coding start sites to calculate DNA positions
epitope_regions <- epitope_regions %>%
  left_join(coding_start_sites, by = "seqname") %>%
  mutate(
    dna_start = (aa_start ) * 3 + start,
    dna_end = (aa_end ) * 3 + start
  )

#####epitope_region  <- epitope_regions  %>%
#####rename(gene=seqname)


epitope_regions <- as.data.frame(epitope_regions)
epitope_regions$seqname <- as.character(epitope_regions$seqname)
# Rename 'seqname' to 'gene'
names(epitope_regions)[names(epitope_regions) == "seqname"] <- "gene"

epitope_region <- epitope_regions
# Print the updated dataframe
print(epitope_region)




# Compute the max depth value for ymax
max_depth <- max(data$depth, na.rm = TRUE)



# Subset data to include only NP, PB1, and HA
data_subset <- data %>%
  filter(gene %in% c("NP", "PB1", "HA")) %>%
  mutate(gene = factor(gene, levels = c("NP", "PB1", "HA")))  # Set the correct order

# Subset epitope_regions accordingly
epitope_regions_subset <- epitope_regions %>%
  filter(gene %in% c("NP", "PB1", "HA")) %>%
  mutate(gene = factor(gene, levels = c("NP", "PB1", "HA")))  # Ensure same order in epitope_regions


# Compute the max depth value for ymax
max_depth <- max(data_subset$depth, na.rm = TRUE)

# Generate the plot
p <- ggplot(data_subset, aes(x = position, y = depth)) +
  geom_line(color = "black", size = 0.2) +
  geom_hline(yintercept = 200, linetype = "dashed", color = "gray", size=0.2) + 
  facet_wrap(~ gene, scales = "free_x", ncol = 3) +  # Ensure 3 columns for NP, PB1, HA
  labs(title = " ",
       x = "Genome Position",
       y = "Read Depth") +
  axis_formatting + background_formatting # Ensure clean formatting
theme(strip.text = element_text(size = 12, face = "bold"))  # Format facet titles

# Add highlighted epitope regions with dynamic ymax
p <- p + geom_rect(data = epitope_regions_subset, inherit.aes = FALSE,
                   aes(xmin = dna_start, xmax = dna_end, ymin = 0, ymax = max_depth),
                   fill = "tan", alpha = 0.3)

# Display plot
print(p)

setwd(dir_save)
# Save the plot
ggsave("1920B007_genome_coverage.png", plot = p, width = 7.08, height = 3, dpi = 320)











