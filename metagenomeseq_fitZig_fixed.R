options(stringsAsFactors = FALSE)

library(metagenomeSeq)
library(biomformat)
library(stringr)

source("format_data.R")


# load in data
counts_data = loadMeta("KO.txt")

taxonomy_table = data.frame(KO=counts_data$taxa$taxa, KO)

phenotypes = loadPhenoData("samples_metadata.txt", tran=TRUE)
ord = match(colnames(counts_data$counts), rownames(phenotypes)) # reorder to match counts
phenotypes = phenotypes[ord, ]

phenotypeData = AnnotatedDataFrame(phenotypes)
KOdata = AnnotatedDataFrame(taxonomy_table)
obj = newMRexperiment(counts_data$counts, phenoData = phenotypeData, featureData = KOdata)


## Normalisation
# this doesn't work on this dataset as col 58 (Neg3) has one feature with a count > 0... 
p = cumNormStatFast(obj)
obj = cumNorm(obj, p)

# OR can use default 0.5
#obj = cumNorm(obj, 0.5)


# using wrench normalisation 
# see: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5160-5
# condition = phenotypes$Description
# library(Wrench)
# wrenchNorm <- function(obj, condition) {
#   count_data <- MRcounts(obj, norm = FALSE)
#   W <- wrench(count_data, condition = condition)
#   normFactors(obj) <- W$nf
#   return(obj)
# }
# 
# obj = wrenchNorm(obj, condition)

normalised_counts = MRcounts(obj, norm = TRUE, log = TRUE)
exportMat(normalised_counts, file = paste0("cumnorm_normalised_counts.txt"))

## Statistical Testing
# using fitZig (Zero-inflated Gaussian mixture model)
# fitFeatureModel doesn't work nicely with multiple conditions

# trim on controls
objTrim = obj
keepFeatures = which(rowSums(MRcounts(objTrim) > 0) >= 5) # keep features with >0 counts in >= 5 samples
objTrim = objTrim[keepFeatures, ]
objTrim = cumNorm(objTrim, p=cumNormStatFast(objTrim))

all_groups = factor(pData(objTrim)$Description)
mod = model.matrix(~0+all_groups)
colnames(mod) = levels(all_groups)

# fitZig model fitting
res = fitZig(obj = objTrim, mod = mod)
zigFit = res@fit
finalMod = res@fit$design

# run diff abundance testing on all combinations...
all_combinations = expand.grid(all_groups, all_groups)
all_combinations = all_combinations[!duplicated(with(all_combinations, paste(Var1,Var2))),]
all_combinations = all_combinations[all_combinations$Var1 != all_combinations$Var2,]
ordered_c = apply(all_combinations, 1, function(x) paste(sort(x), collapse = "-"))
all_combinations = all_combinations[which(!(duplicated(ordered_c))),]

for(i in 1:nrow(all_combinations)){
  group1 = as.character(all_combinations$Var1[i])
  group2 = as.character(all_combinations$Var2[i])

  myContrasts = paste0(group1,"-",group2)
  # contrasts (same as in limma)
  contrast.matrix = makeContrasts(myContrasts, levels = finalMod)
  fit = contrasts.fit(zigFit, contrast.matrix)
  fit = eBayes(fit)
  diffAbundant = (topTable(fit, number=9999999, p.value = 0.05)) # only write features with pvalue < 0.05
  diffAbundant = cbind(feature_id = taxonomy_table$KO.1[match(row.names(diffAbundant), (taxonomy_table$KO))], diffAbundant)
  write.table(diffAbundant, file = paste0("DiffAbundant_fitZig_", group1, "_v_", group2, "_", ".txt"), sep="\t", quote=F)
}
  
  



