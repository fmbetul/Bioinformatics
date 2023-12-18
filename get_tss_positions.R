# Get TSS Positions
library(EnsDb.Hsapiens.v86)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# get gene in the positive strand
pos <- subset(annotation, strand == "+")
pos <- pos[order(pos@ranges@start),] 
# remove all but the first exons per transcript
pos <- pos[!duplicated(pos@elementMetadata@listData$tx_id),] 
# convert to df
pos <- as.data.frame(pos, row.names = NULL)
# make a 1 base pair marker of the TSS
pos$end <- pos$start + 1 

# get genes on the negative strand
neg <- subset(annotation, strand == "-")
# convert to df
neg <- as.data.frame(neg, row.names = NULL)
neg <- neg[order(neg$start, decreasing = TRUE),] 
# remove all but the first exons per transcript
neg <- neg[!duplicated(neg$tx_id),] 
neg$start <- neg$end - 1

tss_positions <- rbind(pos, neg)

# rename columns
names(tss_positions)[names(tss_positions) == "seqnames"] <- "chromosome"
names(tss_positions)[names(tss_positions) == "gene_name"] <- "gene"

# Make a subset of the TSS annotation columns containing just the coordinates 
# and the gene name
tss_positions <- tss_positions[,c("chromosome", "start", "end", "gene")]
