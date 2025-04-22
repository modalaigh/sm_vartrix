# Read in required packages
library(Matrix)
suppressPackageStartupMessages(library(tidyverse))

# Integrate script with snakemake pipeline
vartrix <- as.data.frame(t(as.matrix(readMM(snakemake@input[[1]]))))
barcodes <- read.csv(snakemake@input[[2]], sep = "\t", header = FALSE)
sample_name <- snakemake@wildcards[[1]]
results <- snakemake@output[[1]]

# Recode the vartrix results into a useful format
#vartrix$barcode <- paste0(barcodes$V1, "-1")
vartrix$barcode <- barcodes$V1

df <- vartrix %>% mutate(status = case_when(
  V1 %in% c(0) ~ "no_coverage",
  V1 %in% c(1) ~ "non-mutant",
  V1 %in% c(2, 3) ~ "mutant"
)) %>%
  select(barcode, status) %>%
  mutate(sample = sample_name)

# Save results
write.table(df, results, row.names = FALSE, quote = FALSE, sep = "\t")
