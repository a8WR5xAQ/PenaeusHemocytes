# Cross-species analysis of scRNA-seq data (python)
## Library import
```python
from samap import q, ut, pd, sp, np, warnings, sc
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles, find_cluster_markers)
from samap.utils import (save_samap, load_samap,
                            prepend_var_prefix, df_to_dict, to_vn, 
                            to_vo, substr,
                            sparse_knn)
from samalg import SAM
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
import anndata as ad
from scipy import sparse
from anndata import AnnData
import networkx as nx
import itertools
import sklearn.utils.sparsefuncs as sf
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec
import matplotlib as mpl
from scipy.sparse import csr_matrix
```

## Read data & Data normalization (For *P. vannamei*)
```python
path='path/to/data'
sample_metadata = [
    {'sample_id': 'shrimp_1', 'path': path+'/matrixdata/shrimp1', 'method': 'Dropseq', 'state': 'Live'},
    {'sample_id': 'shrimp_2', 'path': path+'/matrixdata/shrimp2', 'method': 'Dropseq', 'state': 'Live'},
    {'sample_id': 'shrimp_3', 'path': path+'/matrixdata/shrimp3', 'method': 'Dropseq', 'state': 'Live'},
    {'sample_id': 'shrimp_4', 'path': path+'/matrixdata/shrimp4', 'method': 'Dropseq', 'state': 'Live'},
    {'sample_id': 'shrimp_5', 'path': path+'/matrixdata/shrimp5', 'method': 'Dropseq', 'state': 'Live'},
    {'sample_id': 'shrimp10X', 'path': path+'/matrixdata/shrimp10X', 'method': '10X', 'state': 'Live'},
]

adatas = []
for sample_info in sample_metadata:
    sample_id = sample_info['sample_id']
    data_path = sample_info['path']
    temp_adata = sc.read_mtx(data_path+'/matrix.mtx.gz').T
    barcodes_df = pd.read_csv(data_path + '/barcodes.tsv.gz', header=None, sep='\t')
    features_df = pd.read_csv(data_path + '/features.tsv.gz', header=None, sep='\t')
    temp_adata.obs_names = barcodes_df[0].values
    temp_adata.var_names = features_df[1].values
    temp_adata.var['gene_ids'] = features_df[0].values
    temp_adata.var_names_make_unique()
    sc.experimental.pp.normalize_pearson_residuals(temp_adata)
    temp_adata.X = np.nan_to_num(temp_adata.X, nan=0.0)
    temp_adata.X = csr_matrix(temp_adata.X)

    for key, value in sample_info.items():
        if key != 'path':
            temp_adata.obs[key] = value

    adatas.append(temp_adata)

combined_adata = ad.AnnData.concatenate(*adatas)

combined_adata.write(path+'/h5ad/Vannamei.h5ad')
```
Do the same for *P. monodon*, *P. japonicus*, & *D. melanogaster*.
## Read data & SAM
```python
path='path/to/data'
vannamei = sc.read(path+'/h5ad/Vannamei.h5ad')
kuruma = sc.read(path+'/h5ad/Kuruma.h5ad')
monodon = sc.read(path+'/h5ad/Monodon.h5ad')
drosophila = sc.read(path+'/h5ad/Drosophila.h5ad')

sam1=SAM(counts=vannamei)
sam2=SAM(counts=kuruma)
sam3=SAM(counts=monodon)
sam4=SAM(counts=drosophila)

sam1.run(batch_key='batch')
sam2.run(batch_key='batch')
sam3.run(batch_key='batch')
sam4.run()
```

## SAMap
```python
sams = {'va':sam1,'ku':sam2,'mo':sam3,'do':sam4}

sm = SAMAP(
        sams,
        f_maps=path+'/blast_maps/',
    )

sm.run(pairwise=True)
samap = sm.samap
```

## Louvain clustering
```python
sc.tl.louvain(adata=sm.samap.adata,key_added='samap_louvain', resolution=0.7)

sm.sams['va'].adata.obs['samap_louvain']=sm.samap.adata.obs.loc[sm.samap.adata.obs['species'] == 'va']['samap_louvain']
sm.sams['va'].adata.obsm['X_umap_old']=sm.sams['va'].adata.obsm['X_umap']
sm.sams['va'].adata.obsm['X_umap']=sm.sams['va'].adata.obsm['X_umap_samap']
fig, ax = plt.subplots()
sc.pl.umap(sm.sams['va'].adata, size=10, color='samap_louvain', ax=ax, return_fig=fig)
plt.savefig(path+'/figures/umap_louvain_va.svg')
plt.show()

sm.sams['ku'].adata.obs['samap_louvain']=sm.samap.adata.obs.loc[sm.samap.adata.obs['species'] == 'ku']['samap_louvain']
sm.sams['ku'].adata.obsm['X_umap_old']=sm.sams['ku'].adata.obsm['X_umap']
sm.sams['ku'].adata.obsm['X_umap']=sm.sams['ku'].adata.obsm['X_umap_samap']
sc.pl.umap(sm.sams['ku'].adata, size=10, color='samap_louvain', ax=ax, return_fig=fig)
plt.savefig(path+'/figures/umap_louvain_ku.svg')
plt.show()

sm.sams['mo'].adata.obs['samap_louvain']=sm.samap.adata.obs.loc[sm.samap.adata.obs['species'] == 'mo']['samap_louvain']
sm.sams['mo'].adata.obsm['X_umap_old']=sm.sams['mo'].adata.obsm['X_umap']
sm.sams['mo'].adata.obsm['X_umap']=sm.sams['mo'].adata.obsm['X_umap_samap']
sc.pl.umap(sm.sams['mo'].adata, size=10, color='samap_louvain', ax=ax, return_fig=fig)
plt.savefig(path+'/figures/umap_louvain_mo.svg')
plt.show()

sm.sams['do'].adata.obs['samap_louvain']=sm.samap.adata.obs.loc[sm.samap.adata.obs['species'] == 'do']['samap_louvain']
sm.sams['do'].adata.obsm['X_umap_old']=sm.sams['do'].adata.obsm['X_umap']
sm.sams['do'].adata.obsm['X_umap']=sm.sams['do'].adata.obsm['X_umap_samap']
sc.pl.umap(sm.sams['do'].adata, size=10, color='samap_louvain', ax=ax, return_fig=fig)
plt.savefig(path+'/figures/umap_louvain_do.svg')
plt.show()
```

## Mapping score
Modify 'CellTypeTriangle()' function for analysis with 4 species.
```python
def CellTypeSquare(sm,keys, align_thr=0.1): #modification of 'CellTypeTriangle()'
    '''Outputs a table of cell type triangles.
    
    Parameters
    ----------
    sm: SAMAP object - assumed to contain at least three species.
       
    keys: dictionary of annotation keys (`.adata.obs[key]`) keyed by  species.

    align_thr: float, optional, default, 0.1
        Only keep triangles with minimum `align_thr` alignment score.        
    '''
    
    D,A = get_mapping_scores(sm,keys=keys)
    x,y = A.values.nonzero()
    all_pairsf = np.array([A.index[x],A.columns[y]]).T.astype('str')
    alignmentf = A.values[x,y].flatten()

    alignment = alignmentf.copy()
    all_pairs = all_pairsf.copy()
    all_pairs = all_pairs[alignment > align_thr]
    alignment = alignment[alignment > align_thr]
    all_pairs = to_vn(np.sort(all_pairs, axis=1))

    x, y = substr(all_pairs, ';')
    ctu = np.unique(np.concatenate((x, y)))
    Z = pd.DataFrame(data=np.arange(ctu.size)[None, :], columns=ctu)
    nnm = sparse.lil_matrix((ctu.size,) * 2)
    nnm[Z[x].values.flatten(), Z[y].values.flatten()] = alignment
    nnm[Z[y].values.flatten(), Z[x].values.flatten()] = alignment
    nnm = nnm.tocsr()

    import networkx as nx

    G = nx.Graph()
    gps=ctu[np.vstack(nnm.nonzero()).T]
    G.add_edges_from(gps)
    alignment = pd.Series(index=to_vn(gps),data=nnm.data)
    all_triangles = [c for c in nx.enumerate_all_cliques(G) if len(c)==4]
    Z = np.sort(np.vstack(all_triangles), axis=1)
    DF = pd.DataFrame(data=Z, columns=[x.split('_')[0] for x in Z[0]])
    for i,sid1 in enumerate(sm.ids):
        for sid2 in sm.ids[i:]:
            if sid1!=sid2:
                DF[sid1+';'+sid2] = [alignment[x] for x in DF[sid1].values.astype('str').astype('object')+';'+DF[sid2].values.astype('str').astype('object')]
    DF = DF[sm.ids]
    return DF, G, alignment


keys={'va':'samap_louvain','ku':'samap_louvain','mo':'samap_louvain','do':'samap_louvain'}
DF, G, alignment = CellTypeSquare(sm,keys,align_thr=0.05)

data=alignment
edges_with_weights = []
for edge, weight in data.items():
    source, target = edge.split(';')
    edges_with_weights.append((source, target, {'weight': weight}))

G2 = nx.Graph()
G2.add_edges_from(edges_with_weights)

node_names = list(G2.nodes())
sorted_node_names = sorted(node_names)
sorted_G = nx.Graph()
for node in sorted_node_names:
    sorted_G.add_node(node)
for u, v, data in G2.edges(data=True):
    if u in sorted_node_names and v in sorted_node_names:
        sorted_G.add_edge(u, v, weight=data['weight'])


plt.figure(figsize=(13, 10)) 
pos = nx.circular_layout(sorted_G)
edge_colors = [d['weight'] for u, v, d in sorted_G.edges(data=True)]
cmap = plt.get_cmap('Blues')
node_labels=[0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3]
nx.draw(sorted_G, pos, with_labels=True, edge_color=edge_colors, edge_cmap=cmap, width=2, node_size=600, node_color=node_labels, cmap=mpl.colors.ListedColormap(['cyan','limegreen','dodgerblue','tomato']), font_size=10, font_color='black')
colorbar = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
colorbar.set_array([])
plt.colorbar(colorbar)
plt.rcParams['svg.fonttype'] = 'none'
plt.savefig(path+'/figures/allignmentcom.svg')
plt.show()
```

## Conserved marker gene identification
Modify 'GeneTriangles()' function for analysis with 4 species.
```python
ortholog_pairs = np.load(path'/InParanoidDiamond/ortholog_pairs.npy')

# modification of 'GeneTriangles()'
'''Outputs a table of gene triangles.

Parameters
----------
sm: SAMAP object which contains at least three species

orths: (n x 2) ortholog pairs

keys: dict of strings corresponding to each species annotation column keyed by species, optional, default None
    If you'd like to include information about where each gene is differentially expressed, you can specify the
    annotation column to compute differential expressivity from for each species.

compute_markers: bool, optional, default True
    Set this to False if you already precomputed differential expression for the input keys.

corr_thr: float, optional, default, 0.3
    Only keep triangles with minimum `corr_thr` correlation.

pval_thr: float, optional, defaul, 1e-10
    Consider cell types as differentially expressed if their p-values are less than `pval_thr`.
'''
sm=sm
orth=ortholog_pairs
keys=keys
compute_markers=True
corr_thr=0.3
psub_thr = 0.3
pval_thr=1e-10

FINALS = []

orth = np.sort(orth,axis=1)
orthsp = np.vstack([q([x.split('_')[0] for x in xx]) for xx in orth])

RES = ParalogSubstitutions(sm, orth, psub_thr = psub_thr)
if RES.shape[0] > 0:
    op = to_vo(q(RES['ortholog pairs']))
    pp = to_vo(q(RES['paralog pairs']))
    ops = np.vstack([q([x.split('_')[0] for x in xx]) for xx in op])
    pps = np.vstack([q([x.split('_')[0] for x in xx]) for xx in pp])
    doPsubsAll=True
else:
    doPsubsAll=False
gnnm = sm.samap.adata.varp['homology_graph_reweighted']
gn = q(sm.samap.adata.var_names)
gnsp = q([x.split('_')[0] for x in gn])

import itertools
combs = list(itertools.combinations(sm.ids,3))
for comb in combs:
    A,B,C = comb
    smp1 = SAM(counts=sm.samap.adata[np.logical_or(sm.samap.adata.obs['species']==A,sm.samap.adata.obs['species']==B)])
    smp2 = SAM(counts=sm.samap.adata[np.logical_or(sm.samap.adata.obs['species']==A,sm.samap.adata.obs['species']==C)])
    smp3 = SAM(counts=sm.samap.adata[np.logical_or(sm.samap.adata.obs['species']==B,sm.samap.adata.obs['species']==C)])

    sam1=sm.sams[A]
    sam2=sm.sams[B]
    sam3=sm.sams[C]
    A1,A2=A,B
    B1,B2=A,C
    C1,C2=B,C

    f1 = ((orthsp[:,0]==A1) * (orthsp[:,1]==A2) + (orthsp[:,0]==A2) * (orthsp[:,1]==A1)) > 0
    f2 = ((orthsp[:,0]==B1) * (orthsp[:,1]==B2) + (orthsp[:,0]==B2) * (orthsp[:,1]==B1)) > 0
    f3 = ((orthsp[:,0]==C1) * (orthsp[:,1]==C2) + (orthsp[:,0]==C2) * (orthsp[:,1]==C1)) > 0
    orth1 = orth[f1]
    orth2 = orth[f2]
    orth3 = orth[f3]

    gnnm1 = sp.sparse.vstack((
                                sp.sparse.hstack((sp.sparse.csr_matrix(((gnsp==A1).sum(),)*2),gnnm[gnsp==A1,:][:,gnsp==A2])),
                                sp.sparse.hstack((gnnm[gnsp==A2,:][:,gnsp==A1],sp.sparse.csr_matrix(((gnsp==A2).sum(),)*2)))
                                )).tocsr()
    gnnm2 = sp.sparse.vstack((
                                sp.sparse.hstack((sp.sparse.csr_matrix(((gnsp==B1).sum(),)*2),gnnm[gnsp==B1,:][:,gnsp==B2])),
                                sp.sparse.hstack((gnnm[gnsp==B2,:][:,gnsp==B1],sp.sparse.csr_matrix(((gnsp==B2).sum(),)*2)))
                            )).tocsr()
    gnnm3 = sp.sparse.vstack((
                                sp.sparse.hstack((sp.sparse.csr_matrix(((gnsp==C1).sum(),)*2),gnnm[gnsp==C1,:][:,gnsp==C2])),
                                sp.sparse.hstack((gnnm[gnsp==C2,:][:,gnsp==C1],sp.sparse.csr_matrix(((gnsp==C2).sum(),)*2)))
                            )).tocsr()
    gn1 = np.append(gn[gnsp==A1],gn[gnsp==A2])
    gn2 = np.append(gn[gnsp==B1],gn[gnsp==B2])
    gn3 = np.append(gn[gnsp==C1],gn[gnsp==C2])

    if doPsubsAll:
        f1 = np.logical_and(((ops[:,0]==A1) * (ops[:,1]==A2) + (ops[:,0]==A2) * (ops[:,1]==A1)) > 0,
                            ((pps[:,0]==A1) * (pps[:,1]==A2) + (pps[:,0]==A2) * (pps[:,1]==A1)) > 0)
        f2 = np.logical_and(((ops[:,0]==B1) * (ops[:,1]==B2) + (ops[:,0]==B2) * (ops[:,1]==B1)) > 0,
                            ((pps[:,0]==B1) * (pps[:,1]==B2) + (pps[:,0]==B2) * (pps[:,1]==B1)) > 0)
        f3 = np.logical_and(((ops[:,0]==C1) * (ops[:,1]==C2) + (ops[:,0]==C2) * (ops[:,1]==C1)) > 0,
                            ((pps[:,0]==C1) * (pps[:,1]==C2) + (pps[:,0]==C2) * (pps[:,1]==C1)) > 0)
        doPsubs = f1.sum() > 0 and f2.sum() > 0 and f3.sum() > 0
    else:
        doPsubs = False

    if doPsubs:
        RES1=RES[f1]
        RES2=RES[f2]
        RES3=RES[f3]

        op1 = to_vo(q(RES1['ortholog pairs']))
        op2 = to_vo(q(RES2['ortholog pairs']))
        op3 = to_vo(q(RES3['ortholog pairs']))
        pp1 = to_vo(q(RES1['paralog pairs']))
        pp2 = to_vo(q(RES2['paralog pairs']))
        pp3 = to_vo(q(RES3['paralog pairs']))

        # suppress warning
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            T1 = pd.DataFrame(data=np.arange(gn1.size)[None, :], columns=gn1)
            x, y = T1[op1[:, 0]].values.flatten(), T1[op1[:, 1]].values.flatten()
            gnnm1[x, y] = gnnm1[x, y]
            gnnm1[y, x] = gnnm1[y, x]

            T1 = pd.DataFrame(data=np.arange(gn2.size)[None, :], columns=gn2)
            x, y = T1[op2[:, 0]].values.flatten(), T1[op2[:, 1]].values.flatten()
            gnnm2[x, y] = gnnm2[x, y]
            gnnm2[y, x] = gnnm2[y, x]

            T1 = pd.DataFrame(data=np.arange(gn3.size)[None, :], columns=gn3)
            x, y = T1[op3[:, 0]].values.flatten(), T1[op3[:, 1]].values.flatten()
            gnnm3[x, y] = gnnm3[x, y]
            gnnm3[y, x] = gnnm3[y, x]

        gnnm1.data[gnnm1.data==0]=1e-4
        gnnm2.data[gnnm2.data==0]=1e-4
        gnnm3.data[gnnm3.data==0]=1e-4

    pairs1 = gn1[np.vstack(gnnm1.nonzero()).T]
    pairs2 = gn2[np.vstack(gnnm2.nonzero()).T]
    pairs3 = gn3[np.vstack(gnnm3.nonzero()).T]
    data = np.concatenate((gnnm1.data, gnnm2.data, gnnm3.data))

    CORR1 = pd.DataFrame(data=gnnm1.data[None, :], columns=to_vn(pairs1))
    CORR2 = pd.DataFrame(data=gnnm2.data[None, :], columns=to_vn(pairs2))
    CORR3 = pd.DataFrame(data=gnnm3.data[None, :], columns=to_vn(pairs3))

    pairs = np.vstack((pairs1, pairs2, pairs3))
    all_genes = np.unique(pairs.flatten())
    Z = pd.DataFrame(data=np.arange(all_genes.size)[None, :], columns=all_genes)
    x, y = Z[pairs[:, 0]].values.flatten(), Z[pairs[:, 1]].values.flatten()
    GNNM = sp.sparse.lil_matrix((all_genes.size,) * 2)
    GNNM[x, y] = data
    GNNM=GNNM.tocsr()
    GNNM.data[GNNM.data<corr_thr]=0
    GNNM.eliminate_zeros()

    G = nx.from_scipy_sparse_array(GNNM, create_using=nx.Graph)
    all_triangles = [c for c in nx.enumerate_all_cliques(G) if len(c)==3]
    Z = all_genes[np.sort(np.vstack(all_triangles), axis=1)] 

    DF = pd.DataFrame(data=Z, columns=[x.split('_')[0] for x in Z[0]])
    DF = DF[[A, B, C]]

    orth1DF = pd.DataFrame(data=orth1, columns=[x.split('_')[0] for x in orth1[0]])[
        [A, B]
    ]
    orth2DF = pd.DataFrame(data=orth2, columns=[x.split('_')[0] for x in orth2[0]])[
        [A, C]
    ]
    orth3DF = pd.DataFrame(data=orth3, columns=[x.split('_')[0] for x in orth3[0]])[
        [B, C]
    ]

    if doPsubs:
        ps1DF = pd.DataFrame(
            data=np.sort(pp1, axis=1),
            columns=[x.split('_')[0] for x in np.sort(pp1, axis=1)[0]],
        )[[A, B]]
        ps2DF = pd.DataFrame(
            data=np.sort(pp2, axis=1),
            columns=[x.split('_')[0] for x in np.sort(pp2, axis=1)[0]],
        )[[A, C]]
        ps3DF = pd.DataFrame(
            data=np.sort(pp3, axis=1),
            columns=[x.split('_')[0] for x in np.sort(pp3, axis=1)[0]],
        )[[B, C]]

        A_AB = pd.DataFrame(data=to_vn(op1)[None, :], columns=to_vn(ps1DF.values))
        A_AC = pd.DataFrame(data=to_vn(op2)[None, :], columns=to_vn(ps2DF.values))
        A_BC = pd.DataFrame(data=to_vn(op3)[None, :], columns=to_vn(ps3DF.values))
    else:
        ps1DF,ps2DF,ps3DF = None,None,None
        A_AB,A_AC,A_BC = None,None,None

    AB = to_vn(DF[[A, B]].values)
    AC = to_vn(DF[[A, C]].values)
    BC = to_vn(DF[[B, C]].values)

    AVs = []
    CATs = []
    CORRs = []
    for i, X, O, P, Z, R in zip(
        [0, 1, 2],
        [AB, AC, BC],
        [orth1DF, orth2DF, orth3DF],
        [ps1DF, ps2DF, ps3DF],
        [A_AB, A_AC, A_BC],
        [CORR1, CORR2, CORR3],
    ):
        cat = q(['homolog'] * X.size).astype('object')
        cat[np.in1d(X, to_vn(O.values))] = 'ortholog'
        AV = np.zeros(X.size, dtype='object')
        if doPsubs:
            ff = np.in1d(X, to_vn(P.values))
            cat[ff] = 'substitution'

            z = Z[X[ff]] #problem line here
            x = X[ff]
            av = np.zeros(x.size, dtype='object')
            for ai in range(x.size):
                v=pd.DataFrame(z[x[ai]]) #get ortholog pairs - paralog pairs dataframe
                vd=v.values.flatten() #get ortholog pairs
                vc=q(';'.join(v.columns).split(';')) # get paralogous genes
                temp = np.unique(q(';'.join(vd).split(';'))) #get orthologous genes
                av[ai] = ';'.join(temp[np.in1d(temp,vc,invert=True)]) #get orthologous genes not present in paralogous genes
            AV[ff] = av
        corr = R[X].values.flatten()
        AVs.append(AV)
        CATs.append(cat)
        CORRs.append(corr)

    tri_pairs = np.vstack((AB, AC, BC)).T
    cat_pairs = np.vstack(CATs).T
    corr_pairs = np.vstack(CORRs).T
    homology_triangles = DF.values
    substituted_genes = np.vstack(AVs).T
    substituted_genes[substituted_genes == 0] = 'N.S.'
    data = np.hstack(
        (
            homology_triangles.astype('object'),
            substituted_genes.astype('object'),
            tri_pairs.astype('object'),
            corr_pairs.astype('object'),
            cat_pairs.astype('object'),
        )
    )

    FINAL = pd.DataFrame(data = data, columns = [f'{A} gene',f'{B} gene',f'{C} gene',
                                                f'{A}/{B} subbed',f'{A}/{C} subbed',f'{B}/{C} subbed',
                                                f'{A}/{B}',f'{A}/{C}',f'{B}/{C}',
                                                f'{A}/{B} corr',f'{A}/{C} corr',f'{B}/{C} corr',
                                                f'{A}/{B} type',f'{A}/{C} type',f'{B}/{C} type'])
    FINAL['#orthologs'] = (cat_pairs=='ortholog').sum(1)
    FINAL['#substitutions'] = (cat_pairs=='substitution').sum(1)    
    FINAL = FINAL[(FINAL['#orthologs']+FINAL['#substitutions'])==3]
    x = FINAL[[f'{A}/{B} corr',f'{A}/{C} corr',f'{B}/{C} corr']].min(1)
    FINAL['min_corr'] = x
    FINAL = FINAL[x>corr_thr]
    if keys is not None:
        keys3 = [keys[A],keys[B],keys[C]] # 変更
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if keys3 is not None:
            for i,sam,n in zip([0,1,2],[sam1,sam2,sam3],[A,B,C]):
                if compute_markers:
                    find_cluster_markers(sam,keys3[i])
                a = sam.adata.varm[keys3[i]+'_scores'].T[q(FINAL[n+' gene'])].T
                p = sam.adata.varm[keys3[i]+'_pvals'].T[q(FINAL[n+' gene'])].T.values
                p[p>pval_thr]=1
                p[p<1]=0
                p=1-p
                f = a.columns[a.values.argmax(1)]
                res=[]
                for i in range(p.shape[0]):
                    res.append(';'.join(np.unique(np.append(f[i],a.columns[p[i,:]==1]))))
                FINAL[n+' cell type'] = res
    FINAL = FINAL.sort_values('min_corr',ascending=False)
    FINALS.append(FINAL)
FINAL = pd.concat(FINALS,axis=0)
#return FINAL

FINAL.to_csv(path+'/out/GeneTrianglesINPARANOIDDIAMOND.csv')
```
Select conserved marker genes manualy, then saved as 'markergene.csv.'
```python
metadata1 = 'samap_louvain'
def mean_and_fract_cells_by_gene(df, species):
    fraction_of_cells = (df.drop(metadata1, axis=1)>0).sum()/len(df)
    mean = df.mean(numeric_only=True)
    return pd.concat([fraction_of_cells.rename('fraction_of_cells')*400, mean.rename(species)], axis=1).reset_index()

fig, ax = plt.subplots(figsize = (25, 10))

sps=['mo','ku','va']
spl=['P.monodon', 'P.japonicus', 'P.vannamei']
col=['green', 'blue','orangered']

for i in range(3):
    gene = list(markergene[sps[i]+' gene'])
    df = sc.get.obs_df(sm.sams[sps[i]].adata, keys=gene+[metadata1])
    df = df.groupby([metadata1]).apply(lambda df: mean_and_fract_cells_by_gene(df, species=spl[i])).reset_index()
    df['pos']=df2['samap_louvain'].astype(int)*3 + 2 - i

    df.plot.scatter(x='pos', y='level_1', s='fraction_of_cells', 
                    c=spl[i], cmap=mpl.colors.LinearSegmentedColormap.from_list('colormap_name', ['lightgrey',col[i]]), colorbar=True, edgecolors=None,
                    norm=Normalize(vmin=-2, vmax=20), ax=ax);

ax.invert_yaxis()
ax.set_title('Dotplot')
fig.set_tight_layout(True)

plt.tick_params(length=0)

plt.xticks(list(range(21)), ['' if i % 3 !=1 else str(i // 3) for i in list(range(21))])
plt.yticks(list(range(26)), list(markergene['va gene'].str.strip('va_')+'  '+markergene['ku gene'].str.strip('ku_')+'  '+markergene['mo gene'].str.strip('mo_')))

ax.set_xlabel('clusters')
ax.set_ylabel('gene symbols')

plt.rcParams['svg.fonttype'] = 'none'

plt.savefig(path+'/figures/Dotplot.svg')
```