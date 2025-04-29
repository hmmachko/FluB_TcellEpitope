## Clear Global Environment
rm(list = ls())
allelefreq <- .03
depth <-200

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

if(!require(vcfR)){
  install.packages("vcfR")
  library(vcfR)
}


#### Directories ####
dir <-      paste("/Volumes/PrimateFS/TCF lab/Current Lab Members/machkovech_heather/bioinformatics/13197/processed_data/recall_fn_ann/", sep="")
dir_processed <- paste("/Volumes/PrimateFS/TCF lab/Current Lab Members/machkovech_heather/bioinformatics/13197/processed_data/", sep="")
dir_save <-paste("/Volumes/PrimateFS/TCF lab/Current Lab Members/machkovech_heather/bioinformatics/13197/figures/", sep="")
dir_resource <- paste("/Volumes/PrimateFS/TCF lab/Current Lab Members/machkovech_heather/bioinformatics/13197/resources/", sep="")
#dir <-      paste("/Volumes/primatefs/Lab-Friedrich/TCF lab/Current Lab Members/machkovech_heather/bioinformatics/13197/processed_data/recall_fn_ann/", sep="")
#dir_processed <-      paste("/Volumes/primatefs/Lab-Friedrich/TCF lab/Current Lab Members/machkovech_heather/bioinformatics/13197/processed_data/", sep="")
#dir_save <- paste("/Volumes/primatefs/Lab-Friedrich/TCF lab/Current Lab Members/machkovech_heather/bioinformatics/13197/figures/", sep="")
#dir_resource <- paste("/Volumes/primatefs/Lab-Friedrich/TCF lab/Current Lab Members/machkovech_heather/bioinformatics/13197/resources/", sep="")

#### iSNVs: Import, clean ####
## 
setwd(paste(dir)); getwd(); head(dir())

vcf <- dir(pattern=".vcf")
names_trunc <- gsub(".vcf","",vcf)
names_trunc <- gsub("_recall_fn_ann","",names_trunc)
n <- length(vcf)
list <- vector("list",n)

## Read all tables in vcf, apply to list, change columns
for (i in 1:n) {
  list[[i]] <- read.vcfR(vcf[i], verbose=T)
  list[[i]] <- cbind(as.data.frame(list[[i]]@fix), as.data.frame(list[[i]]@gt))
  colnames(list[[i]])[10] <- "AF"
  names(list) <- names_trunc}

## remove blank data frames
list <- Filter(function(x) dim(x)[1] > 0, list)

####################################

# Function to parse individual ANN entries (this can remain mostly the same)
parse_ann_entry <- function(ann_string) {
  fields <- unlist(strsplit(ann_string, "\\|"))
  
  result <- c(
    allele = fields[1],
    annotation = fields[2],
    impact = fields[3],      # Added impact field
    gene_name = fields[4],
    feature_type = fields[6],
    feature_ID = fields[7],
    hgvs_c = fields[10],
    hgvs_p = fields[11]
  )
  
  return(result)
}

# Modified FORMAT fields parser for the new format
parse_format_fields <- function(format_string, genotype_string) {
  format_fields <- unlist(strsplit(format_string, ":"))
  genotype_values <- unlist(strsplit(genotype_string, ":"))
  
  # Create named vector of values
  names(genotype_values) <- format_fields
  
  # Extract required values based on new format
  dp <- as.numeric(genotype_values["DP"])
  ad <- as.numeric(genotype_values["AD"])
  af <- as.numeric(genotype_values["AF"])
  
  return(c(DP = dp, AD = ad, AF = af))
}

# Modified INFO field parser
parse_vcf_annotations <- function(info_string, format_string, genotype_string) {
  # Extract ANN field
  ann_match <- regexpr("ANN=([^;]+)", info_string)
  if (ann_match == -1) return(NULL)
  
  ann_string <- substr(info_string, 
                       ann_match + 4, 
                       ann_match + attr(ann_match, "match.length") - 1)
  
  # Split multiple annotations
  ann_entries <- unlist(strsplit(ann_string, ","))
  
  # Process each annotation
  results <- lapply(ann_entries, parse_ann_entry)
  
  # Convert to data frame
  results_df <- as.data.frame(do.call(rbind, results))
  
  # Parse FORMAT and genotype fields
  format_values <- parse_format_fields(format_string, genotype_string)
  
  # Add FORMAT values to each row
  results_df$DP <- format_values["DP"]
  results_df$AD <- format_values["AD"]
  results_df$AF <- format_values["AF"]
  
  # Add additional INFO fields if needed
  info_fields <- unlist(strsplit(info_string, ";"))
  info_dict <- sapply(info_fields, function(x) {
    parts <- unlist(strsplit(x, "="))
    if(length(parts) == 2) {
      setNames(parts[2], parts[1])
    }
  })
  
  # Add relevant INFO fields to results
  results_df$SB <- as.numeric(info_dict["SB"])
  results_df$NVC <- as.numeric(info_dict["NVC"])
  
  return(results_df)
}

# The process_vcf function can remain largely the same
process_vcf <- function(vcf_df, sample_name) {
  parsed_list <- Map(parse_vcf_annotations, 
                     vcf_df$INFO, 
                     vcf_df$FORMAT, 
                     vcf_df[, ncol(vcf_df)])
  
  result <- data.frame()
  
  for(i in seq_along(parsed_list)) {
    if(!is.null(parsed_list[[i]]) && nrow(parsed_list[[i]]) > 0) {
      temp_df <- cbind(
        Sample = sample_name,
        CHROM = vcf_df$CHROM[i],
        POS = vcf_df$POS[i],
        REF = vcf_df$REF[i],
        ALT = vcf_df$ALT[i],
        QUAL = vcf_df$QUAL[i],
        FILTER = vcf_df$FILTER[i],
        parsed_list[[i]]
      )
      result <- rbind(result, temp_df)
    }
  }
  
  return(result)
}

##############################################
list_mod <- vector("list", length = length(list))

# Preserve names from the original list
names(list_mod) <- names(list)


for (i in seq_along(list)) {
  parsed_results <- process_vcf(list[[i]], names(list)[i])  # Pass the sample name
  list_mod[[i]] <- parsed_results
}

###################################################################

# join and reduce all samples into one df

df_VCF <- Reduce(full_join,list_mod)
df_VCF$DP <- as.integer(df_VCF$DP)
df_VCF$AF <- as.numeric(df_VCF$AF)
df_VCF$POS <- as.integer(df_VCF$POS)


#### filtering 
df_VCF <- df_VCF[df_VCF$AF >= allelefreq,]
df_VCF <- df_VCF[df_VCF$DP >= depth,]

df_VCF <- df_VCF %>%
  rename(substitution_type = annotation)

## Clarify mutation type, keep only NS and S and inframe deletion
df_VCF$substitution_type <- as.factor(df_VCF$substitution_type); levels(df_VCF$substitution_type)
df_VCF$substitution_type <- gsub("missense_variant","nonsynonymous",df_VCF$substitution_type)
df_VCF$substitution_type <- gsub("synonymous_variant","synonymous",df_VCF$substitution_type)
df_VCF$substitution_type <- gsub("conservative_inframe_deletion","deletion",df_VCF$substitution_type)
df_VCF_filtered <- df_VCF %>%
  filter(substitution_type %in% c("nonsynonymous", "synonymous","deletion" ))

#df_VCF_filtered <- df_VCF_filtered %>%
#rename(Sample = FILTER)

df_VCF_filtered <- df_VCF_filtered %>%
  rename(segment = CHROM)

df_VCF_filtered <- df_VCF_filtered %>%
  rename(gene = gene_name)

df_VCF_filtered <- df_VCF_filtered %>%
  select(Sample, segment, gene, POS, substitution_type, REF, ALT, DP, AD, AF,hgvs_c, hgvs_p,QUAL )


#save outfile for Supplemental table
write.csv(df_VCF_filtered, file = file.path(dir_processed, "processed_variants_table.csv"), row.names = FALSE)

#end for vcf processing 
##########################################################################

