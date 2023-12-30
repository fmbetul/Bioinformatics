# get_rsid_from_hg19_positions

# Connect to the Ensembl GRCh37 BioMart
ensembl37 = useMart(host = "https://grch37.ensembl.org", biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# Assuming GWAS is the data frame containing hg19 SNP positions
# Add 'chromosomal_region' column to GWAS (eg. 17:51393:51393)
GWAS[["chromosomal_region"]] <- paste(GWAS$CHR, GWAS$BP, GWAS$BP, sep = ":")

# Adjust batch size if needed
split_snps <- split(GWAS$chromosomal_region, ceiling(seq_along(GWAS$chromosomal_region)/10)) 

# Initialize the empty results data frame
results <- data.frame(
  refsnp_id = character(), 
  chr_name = integer(),
  chrom_start = integer(),
  chrom_end = integer(),
  stringsAsFactors = FALSE
)

# Iterate through each batch in split_snps with error handling
for(i in 1:length(split_snps)) {
  tryCatch({
    # Fetch results for the current batch
    batch_result <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"),
                          filters = "chromosomal_region",
                          values = split_snps[[i]],
                          mart = ensembl37)
    
    # Combine with the accumulated results
    results <- rbind(results, batch_result)
  }, error = function(e) {
    message(sprintf("Error in batch %d: %s", i, e$message))
  })
}
