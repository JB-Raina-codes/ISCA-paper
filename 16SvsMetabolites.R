library(pheatmap)
library(tidyverse)
library(psych)

# Loading the 16S data
data_16S <- read.csv(file='16S.csv')

#Loading the metabolome data
metabolome <- read.csv(file='metabolome.csv',colClasses = "character")

# get mean values for each condition...
conditions = unique(unlist(lapply(str_split(colnames(metabolome), "_"),"[[", 1)))
# you have metabolites in the 'X' column. It used to be called Name, which is why you weren't getting metabolite names in the heatmaps
conditions = conditions[!(conditions %in% c("Name", "X"))]
metabolome_numeric = mutate_all(metabolome[,-1], function(x) as.numeric(as.character(x)))
for(i in seq_along(conditions)){

  metabolome_numeric = mutate(metabolome_numeric, mean = rowMeans(select(metabolome_numeric, starts_with(conditions[i])), na.rm = TRUE)) 
  
  colnames(metabolome_numeric)[ncol(metabolome_numeric)] = paste0("mean_", conditions[i])
   
}
mean_metabolome = metabolome_numeric[,grep("mean_", colnames(metabolome_numeric))]
rownames(mean_metabolome) = metabolome[,1]
metabolome_ind = metabolome_numeric[,-grep("mean_", colnames(metabolome_numeric))]
rownames(metabolome_ind) = metabolome[,1]


# get mean values for each condition...
conditions = unique(unlist(lapply(str_split(colnames(data_16S), "_"),"[[", 1)))
conditions = conditions[conditions != "X"]
data_16S_numeric = mutate_all(data_16S[,-1], function(x) as.numeric(as.character(x)))
for(i in seq_along(conditions)){
  
  data_16S_numeric = mutate(data_16S_numeric, mean = rowMeans(select(data_16S_numeric, starts_with(conditions[i])), na.rm = TRUE)) 
  
  colnames(data_16S_numeric)[ncol(data_16S_numeric)] = paste0("mean_", conditions[i])
  
}
mean_data_16S = data_16S_numeric[,grep("mean_", colnames(data_16S_numeric))]
rownames(mean_data_16S) = data_16S$X
data_16S_ind = data_16S_numeric[,-grep("mean_", colnames(data_16S_numeric))]
rownames(data_16S_ind) = data_16S$X

# matrices need to have the same columns
conditions_in_both = colnames(mean_metabolome)[which(colnames(mean_metabolome) %in% colnames(mean_data_16S))]
mean_metabolome = mean_metabolome[,conditions_in_both]
mean_data_16S = mean_data_16S[,conditions_in_both]

# remove 0's
mean_data_16S = mean_data_16S[which(rowSums(mean_data_16S) > 0),]
mean_metabolome = mean_metabolome[which(rowSums(mean_metabolome) > 0),]


# loop through correlations of metabolites to 16s data
# no within-group correlations 
# (i.e. Root; k__Bacteria; p__Spirochaetes VS Root; k__Bacteria; p__FCPU426 is not computed)
t.cor_all = NULL
t.pvals.noadjust_all = NULL
for(i in 1:nrow(mean_metabolome)){
  
  t.cor = cor(t(mean_metabolome[i,]), t(mean_data_16S), method = 'pearson')
  t.pvals.noadjust = psych::corr.test(t(mean_metabolome[i,]), t(mean_data_16S), ci=FALSE, adjust = "none")
  
  t.cor_all = rbind(t.cor_all, t.cor)
  t.pvals.noadjust_all = rbind(t.pvals.noadjust_all, t.pvals.noadjust$p)
  
}


row_hclust  = hclust(dist(t.cor_all))
col_hclust  = hclust(dist(t(t.cor_all)))

# all correlations (no filtering)
pdf("pheatmap_all.pdf", height= 30, width=100, useDingbats = F)
pheatmap(t.cor_all, cluster_rows = row_hclust, cluster_cols = col_hclust)
dev.off()

# only p.val < 0.05 correlations
t.cor_all_signif= t.cor_all
t.cor_all_signif[(t.pvals.noadjust_all > 0.05)] = NA
pdf("pheatmap_p0.05.pdf", height= 30, width=100, useDingbats = F)
pheatmap(t.cor_all_signif, cluster_rows = row_hclust, cluster_cols = col_hclust)
dev.off()

# only ADJUSTED p.val < 0.05 correlations
adj = matrix(p.adjust(t.pvals.noadjust_all), nrow = nrow(t.pvals.noadjust_all), ncol=ncol(t.pvals.noadjust_all))
t.cor_all_signif[(adj > 0.05)] = NA
pdf("pheatmap_adjustedp0.05.pdf", height= 30, width=100, useDingbats = F)
pheatmap(t.cor_all_signif, cluster_rows = row_hclust, cluster_cols = col_hclust)
dev.off()



head(t.cor_all_signif)
signif_correlations = NULL
for(i in 1:ncol(t.cor_all_signif)){
  
  signif_rows = which(!is.na(t.cor_all_signif[,i]))
  
  if(length(signif_rows) > 0){
    signif_correlations = rbind(signif_correlations,
    data.frame(col = colnames(t.cor_all_signif)[i], col_n = i, 
               row = names(signif_rows),
               row_n = as.numeric(signif_rows), 
               correlation = t.cor_all_signif[,i][signif_rows],
               p_value_adj = adj[,i][signif_rows])
    )
  }

  
}

write_csv(signif_correlations, "significant_correlations.csv")

##################################################################
# using individual values (not averaged across each condition)
##################################################################

# matrices need to have the same columns
conditions_in_both = colnames(metabolome_ind)[which(colnames(metabolome_ind) %in% colnames(data_16S_ind))]
metabolome_ind = metabolome_ind[,conditions_in_both]
data_16S_ind = data_16S_ind[,conditions_in_both]

# remove 0's
data_16S_ind = data_16S_ind[which(rowSums(data_16S_ind) > 0),]
metabolome_ind = metabolome_ind[which(rowSums(metabolome_ind) > 0),]


# loop through correlations of metabolites to 16s data
# no within-group correlations 
# (i.e. Root; k__Bacteria; p__Spirochaetes VS Root; k__Bacteria; p__FCPU426 is not computed)
t.cor_all = NULL
t.pvals.noadjust_all = NULL
for(i in 1:nrow(metabolome_ind)){
  
  t.cor = cor(t(metabolome_ind[i,]), t(data_16S_ind), method = 'pearson')
  t.pvals.noadjust = psych::corr.test(t(metabolome_ind[i,]), t(data_16S_ind), ci=FALSE, adjust = "none")
  
  t.cor_all = rbind(t.cor_all, t.cor)
  t.pvals.noadjust_all = rbind(t.pvals.noadjust_all, t.pvals.noadjust$p)
  
}

# or pheatmaps
# these should be 'small' enough to make a pdf from and still be able to read col/rownames
pheatmap(t.cor_all)

row_hclust  = hclust(dist(t.cor_all))
col_hclust  = hclust(dist(t(t.cor_all)))

# all correlations (no filtering)
pdf("BS_pheatmap_all_ind.pdf", height= 30, width=100, useDingbats = F)
pheatmap(t.cor_all, cluster_rows = row_hclust, cluster_cols = col_hclust)
dev.off()

# only p.val < 0.05 correlations
t.cor_all_signif= t.cor_all
t.cor_all_signif[(t.pvals.noadjust_all > 0.05)] = NA
pdf("BS_pheatmap_p0.05_ind.pdf", height= 30, width=100, useDingbats = F)
pheatmap(t.cor_all_signif, cluster_rows = row_hclust, cluster_cols = col_hclust)
dev.off()

# only ADJUSTED p.val < 0.05 correlations
adj = matrix(p.adjust(t.pvals.noadjust_all), nrow = nrow(t.pvals.noadjust_all), ncol=ncol(t.pvals.noadjust_all))
t.cor_all_signif[(adj > 0.05)] = NA
pdf("BS_pheatmap_adjustedp0.05_ind.pdf", height= 30, width=100, useDingbats = F)
pheatmap(t.cor_all_signif, cluster_rows = row_hclust, cluster_cols = col_hclust)
dev.off()



head(t.cor_all_signif)
signif_correlations = NULL
for(i in 1:ncol(t.cor_all_signif)){
  
  signif_rows = which(!is.na(t.cor_all_signif[,i]))
  
  if(length(signif_rows) > 0){
    signif_correlations = rbind(signif_correlations,
                                data.frame(col = colnames(t.cor_all_signif)[i], col_n = i, 
                                           row = names(signif_rows),
                                           row_n = as.numeric(signif_rows), 
                                           correlation = t.cor_all_signif[,i][signif_rows],
                                           p_value_adj = adj[,i][signif_rows])
    )
  }
  
  
}

write_csv(signif_correlations, "significant_correlations_ind.csv")
