library(tidyverse)

setwd("../output")
files <- list.files()

all_annotations <- data.frame()
for (file in files){
    all_annotations <- rbind(all_annotations, read.table(file, sep = "\t", header = T))
}

write.table(all_annotations, "all_annotations.tsv", sep = "\t", col.names = T, row.names = F, quote = F)