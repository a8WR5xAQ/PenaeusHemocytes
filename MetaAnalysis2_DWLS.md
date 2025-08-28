# Cell type deconvolution analysis (R)
## Library import
```{r}
library('DWLS')
library('zellkonverter')
library('zellkonverter')
library('reshape')
library('cowplot')
library('ggplot2')
library('reshape2')
library('nparcomp')
library('coin')
```

## Set working directory & read scRNA-seq data
```{r}
path = 'path/to/data'
setwd(path)

ad <- readH5AD(path+'/h5ad/VannameiWithSamapCluster.h5ad', X_name='counts')
smt <- as.matrix(ad@assays@data$counts)
rownames(smt) <- gsub('_', '-', rownames(smt))
id <- as.character(ad@colData$louvain)
```

## Read bulk data
```{r}
tmp2 <- read.csv(path+'/MetaAnalysis/MetaData2.csv', row.names = 1)
pdata2 <- AnnotatedDataFrame(tmp2)

tmp2 <- read.csv(path+'/MetaAnalysis/ExpMatrixdata', row.names = 1)
m2 <- as.matrix(tmp2)

eset <- new('ExpressionSet', exprs = m2, phenoData = pdata2)#, featureData = fdata)

rownames(tmp2) = gsub('_', '-', rownames(tmp2))
bulkData = tmp2
```
## DWLS
```{r}
Sig = DWLS::buildSignatureMatrixMAST(scdata = smt, id = id, 'results', diff.cutoff = 0.5, pval.cutoff = 0.01)
Genes<-intersect(rownames(Sig),rownames(bulkData))
B<-bulkData[Genes,]
S<-Sig[Genes,]
tr = list('sig'=S,'bulk'=B)

DWLSdf <- data.frame(matrix(NA, nrow=length(colnames(B)), ncol=7))
rownames(DWLSdf) <- colnames(B)
colnames(DWLSdf) <- as.character(0:6)

for (i in colnames(B)) {
  solDWLS<-DWLS::solveDampenedWLS(S,as.matrix(B[i]))
  DWLSdf[i,rownames(as.data.frame(solDWLS))]=as.data.frame(t(solDWLS))
}
```

## Correlation & p-value calculation
```{r}
d2 <- DGEList(m2)
d_UQ2 <- calcNormFactors(d2, method='upperquartile')

markerexp2<-as.data.frame(t(data.frame(d_UQ2$counts)[c('XM_027367792.1','XM_027363492.1','XM_027364731.1','XM_027376661.1','XM_027376190.1','XM_027376169.1','XM_027350837.1','XM_027350836.1','XM_027369772.1','XM_027352840.1','XM_027366024.1','XM_027380774.1'),]))
propcluster2<-DWLSdf

correl2 <- as.data.frame(cor(markerexp2,propcluster2))

p_values2 <- matrix(nrow = ncol(markerexp2), ncol = ncol(propcluster2))
rownames(p_values2) <- names(markerexp2)
colnames(p_values2) <- names(propcluster2)

for (i in seq_along(markerexp2)) {
  for (j in seq_along(propcluster2)) {
    test_result <- cor.test(markerexp2[,i], propcluster2[,j])
    p_values2[i, j] <- test_result$p.value
  }
}
p_values_df2 <- as.data.frame(p_values2)

correl2$gene <- rownames(correl2)
corr_s_melt2 <- reshape2::melt(correl2[nrow(correl2):1,], id.vars = 'gene', variable.name = 'cluster', value.name = 'Correlation')
p_values_df2$gene <- rownames(p_values_df2)
p_values_s_melt2 <- reshape2::melt(p_values_df2[nrow(p_values_df2):1,], id.vars = 'gene', variable.name = 'cluster', value.name = 'p_value')
```

## Dot plot for correlation & p-value
```{r}
combined_df2 <- cbind(corr_s_melt2, p_value=p_values_s_melt2[,'p_value'])#, by = c('gene', 'cluster'))

combined_df2$gene <- factor(combined_df2$gene, levels = unique(combined_df2$gene))
combined_df2$cluster <- factor(combined_df2$cluster, levels = unique(combined_df2$cluster))

correl_dotplot2 <- ggplot(combined_df2, aes(x = cluster, y = gene, color = Correlation, size = -log10(p_value))) +
  geom_point() +
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0, limits = c(-1, 1)) +
  labs(size = '-log10(p-value)', fill = 'Correlation') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  scale_size_area(max_size=10)
correl_dotplot2
```

## Save fig & matrix
```{r}
ggsave(filename=path+'/figure/correl_dotplot2.svg', plot=correl_dotplot2, width=5, height=3.5)

PropCluster <- DWLSdf[c('SRR24125251','SRR24125250','SRR24125249','SRR6917104','SRR11293531','SRR11293530','SRR11293529','SRR13480952','SRR13480963','SRR13480964','SRR24125254','SRR24125253','SRR24125252','SRR11293534','SRR11293533','SRR11293532','SRR6917103','SRR13480948','SRR13480949','SRR13480950','SRR24676009','SRR24676008','SRR24676007','SRR24676016','SRR24676015','SRR24676010','SRR24676014','SRR24676013','SRR24676012','SRR24676011','SRR24676006','SRR24676005','SRR24676004','SRR24676003','SRR13480942','SRR13480943','SRR13480944'),]
colnames(PropCluster) <- colnames(p_values2[,1:7])

metadataDIVAHP <- pData(eset)[c('SRR24125251','SRR24125250','SRR24125249','SRR6917104','SRR11293531','SRR11293530','SRR11293529','SRR13480952','SRR13480963','SRR13480964','SRR24125254','SRR24125253','SRR24125252','SRR11293534','SRR11293533','SRR11293532','SRR6917103','SRR13480948','SRR13480949','SRR13480950','SRR24676009','SRR24676008','SRR24676007','SRR24676016','SRR24676015','SRR24676010','SRR24676014','SRR24676013','SRR24676012','SRR24676011','SRR24676006','SRR24676005','SRR24676004','SRR24676003','SRR13480942','SRR13480943','SRR13480944'),]
MetadataCluster <- cbind(metadataDIVAHP,PropCluster)

write.table(MetadataCluster, path+'/MetaAnalysis/MetadataClusterPropSAMDWLS.csv', append = F,sep = ',', row.names = F, quote = F)
```