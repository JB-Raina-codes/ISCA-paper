options(stringsAsFactors = FALSE)
library(stringr)

## reformat data ##

#reformat metadata (to phenodata)
sample_metadata = read.csv("Samples_metadata.csv")
sample_metadata$Sample_ID = gsub("-","_",gsub(" ", "_", sample_metadata$Sample_ID))
sample_metadata$Description = gsub(" ", "_", sample_metadata$Description)

row.names(sample_metadata) = sample_metadata$Sample_ID
write.table(sample_metadata, file="samples_metadata.txt", sep="\t", quote=F, row.names = F)

# reformat counts
counts_csv = read.csv("KO_withBulk.csv")
KO = counts_csv$KO_ID
counts_csv = cbind(KO=1:nrow(counts_csv), counts_csv[,-1])
colnames(counts_csv) = gsub("[.]", "_", colnames(counts_csv))
colnames(counts_csv)[match(c("mockpos_S37_0194","neg_S38_0194","mockpos_S37","neg_S38"), colnames(counts_csv))] = c("Mock2","Neg5", "Mock1","Neg4")
write.table(counts_csv, file="KO.txt", sep="\t", quote=F, row.names = F)
