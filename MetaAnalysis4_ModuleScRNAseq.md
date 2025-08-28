# Gene module expression in scRNA-seq data (python)
## Library import
```python
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import scanpy as sc, anndata as ad, numpy as np, pandas as pd
from scipy import sparse
from anndata import AnnData
```
## Read gene module & scRNA-seq data
```{r}
path='path/to/data'
df = pd.read_csv(path+'/MetaAnalysis/genemodule.csv')
list(df[df['colors'] == 'red']['gene_id'])
colors = list(pd.read_csv(path+'/MetaAnalysis/MEs.csv', index_col = 0, header = 0).columns)
colors = [s.replace('ME', '') for s in colors]
v10_common = readH5AD(path+'/h5ad/VannameiWithSamapCluster.h5ad', X_name='counts')
```
## Dot plot
```{r}
colorlist = []
for color in colors:
    module = list(df[df['colors'] == color]['gene_id'])
    gene_list = [gene for gene in module if gene in list(v10_common.var_names)]
    v10_common.obs[color+' ('+str(len(gene_list))+')'] = v10_common[:, gene_list].X.sum(axis=1).A1/len(gene_list)
    colorlist.append(color+' ('+str(len(gene_list))+')')#Sig
sc.pl.dotplot(v10_common, colorlist, groupby='louvain', vmax=3, save = path+'/figure/genemodule.svg')
```
