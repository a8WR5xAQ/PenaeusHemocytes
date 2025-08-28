# WGCNA (R)
## Library import
```{r}
library('edgeR')
library('sva')
library('limma')
library('tidyverse')
library('ggplot2')
library('WGCNA')
library('dendextend')
library('ggdendro')
library('gplots')
library('RColorBrewer')path='path/to/data'
library('svglite')
```

## Read data
```{r}
path='path/to/data'
ExpData <- read.csv(path+'/MetaAnalysis/ExpressionData.csv')
as.matrix(ExpData)

MetaData <- read.csv(path+'/MetaAnalysis/MetadataClusterPropSAMDWLS.csv')
as.matrix(MetaData)
```

## Normalize & Combat
```{r}
ExpData <- DGEList(counts=ExpData, group=MetaData$bioproject)
ExpData_UQ <- calcNormFactors(ExpData, method='upperquartile')
ExpUQ <- cpm(ExpData_UQ, log = FALSE, normalized.lib.sizes=TRUE)
rownames(ExpUQ)=ExpData_UQ$genes[, 1]
ExpUQ <- ExpUQ[apply(ExpUQ, 1, function(x) !any(x <= 0)), ]

ExpUQ_Combat = ComBat(dat=ExpUQ, batch=MetaData$bioproject, mod=NULL, par.prior=TRUE, prior.plots=TRUE)
```

## Outlier removal
```{r}
sampleTree = hclust(dist(t(ExpUQ_Combat)), method = 'average')
dend <- as.dendrogram(sampleTree)

dend_data <- dendro_data(dend, type = 'rectangle')

sledai_group_tmp <- tibble(array_ID=labels(dend)) %>% 
    left_join(MetaData, by=c('array_ID'='SRA') )

State = sledai_group_tmp$State
names(State) = sledai_group_tmp$array_ID
sledai_color_palette = c('#E34A33', '#FDBB84', 'grey', 'blue')
sledai_color = sledai_color_palette[unclass(as.factor(State))]

cutHeight = 100000
clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 1)
names(clust) = sampleTree$labels
clust = clust[labels(dend)]

cluster_color_palette = brewer.pal(length(table(clust)), 'Set3')
cluster_color = cluster_color_palette[clust]
names(cluster_color) = names(clust)

dend <- dend %>% 
  set('labels_col', cluster_color) %>% # change color
  set('labels_cex', 0) %>% 
  set('by_labels_branches_col', names(clust[clust == 1]),TF_values = cluster_color_palette[1])
plot(dend, ylim = c(-250, 200000))
#abline(h = cutHeight, col = 'red')
colored_bars(dend = dend, 
             colors = data.frame(sledai=sledai_color, cluster=cluster_color),
             rowLabels = c('State', 'cluster'),
             sort_by_labels_order = FALSE,
             y_shift=-80000, y_scale = 50000)
```
```{r}
keepSamples = (clust[colnames(ExpUQ_Combat)]==1)
notKeptSamples = colnames(ExpUQ_Combat)[!keepSamples]

datExpr = ExpUQ_Combat[,keepSamples]
nGenes = nrow(ExpUQ_Combat)
nSamples = ncol(ExpUQ_Combat)

MetaData_RM <- subset(MetaData, !(SRA %in% notKeptSamples))

datExpr=t(datExpr)
```

## WGCNA
```{r}
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = 'signed', verbose = 5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',type='n',
     main = paste('Scale independence'))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col='red')
abline(h=0.90,col='red')


par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = 'Soft Threshold (power)',
     ylab = 'Scale Free Topology Model Fit, signed R^2',
     main = paste('Scale independence')
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = 'red'
)
abline(h = 0.90, col = 'red')
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = 'Soft Threshold (power)',
     ylab = 'Mean Connectivity',
     type = 'n',
     main = paste('Mean connectivity')
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = 'red')
```
Select 'softPower.'
```{r}
softPower = 14
adjacency = adjacency(datExpr, power = softPower, type='signed')
TOM = TOMsimilarity(adjacency,TOMType = 'signed')
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = 'average')

deepSplit <- 0
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, 
                            distM = dissTOM,
                            deepSplit = deepSplit, 
                            pamStage = TRUE,
                            pamRespectsDendro = TRUE,
                            minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)

MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = 'average')
MEDissThres <- 0.25

merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
MEs <- merge$newMEs
```
## Dot plot for correlation & p-value
```{r}
MetadataCONTvs <- data.frame(DIV1 = ifelse(MetaData_RM$State == 'control', 0, ifelse(MetaData_RM$State == 'DIV1', 1, NA)),
                             VP_24 = ifelse(MetaData_RM$State == 'control', 0, ifelse(MetaData_RM$State == 'VP_24', 1, NA)),
                             VP_48 = ifelse(MetaData_RM$State == 'control', 0, ifelse(MetaData_RM$State == 'VP_48', 1, NA)))

MetadataProp <- cbind(MetadataCONTvs, MetaData_RM[,5:11])
colnames(MetadataProp)[1:3] <- c('DIV1', 'Vpara_24', 'Vpara_48')
moduleTraitCor = cor(MEs, MetadataProp, use = 'p')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitCor <- as.data.frame(moduleTraitCor)
moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)

moduleTraitCor$module <- rownames(moduleTraitCor)
moduleTraitCor_melt <- reshape2::melt(moduleTraitCor[nrow(moduleTraitCor):1,], id.vars = 'module', variable.name = 'trait', value.name = 'Correlation')
moduleTraitPvalue$module <- rownames(moduleTraitPvalue)
moduleTraitPvalue_melt <- reshape2::melt(moduleTraitPvalue[nrow(moduleTraitPvalue):1,], id.vars = 'module', variable.name = 'trait', value.name = 'p_value')

module_df <- cbind(moduleTraitCor_melt, p_value=moduleTraitPvalue_melt[,'p_value'])
module_df$module <- factor(module_df$module, levels = unique(module_df$module))
module_df$trait <- factor(module_df$trait, levels = unique(module_df$trait))

module_dotplot <- ggplot(module_df, aes(x = trait, y = module, color = Correlation, size = -log10(p_value))) +
  geom_point() +
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0, limits = c(-1, 1)) +
  labs(size = '-log10(p-value)', fill = 'Correlation') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_size_area(max_size=7)
module_dotplot
```
## Save fig & matrix
```{r}
ggsave(filename=path+'/figure/module_dotplot.svg', plot=module_dotplot, width=10, height=7)
module_gene_df <- data.frame(gene_id = colnames(datExpr), colors = merge$colors)
write.csv(module_gene_df, path+'/MetaAnalysis/genemodule.csv')
write.csv(MEs, path+'/MetaAnalysis/MEs.csv')
```