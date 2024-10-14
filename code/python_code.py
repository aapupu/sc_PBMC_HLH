import numpy as np
import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import pygam
import seaborn as sns
import anndata
from scipy import io
#################### PAGA
adata = sc.read_h5ad('CD8_counts.h5ad')
adata.obs.Celltype3 = pd.Categorical(adata.obs.Celltype3, 
                                    categories=['CD8_Tn_CCR7','CD8_Tn_BACH2','CD8_Tcm_RPS10',
                                                 'CD8_Tcm_FOS','CD8_Tem_IL7R','CD8_Tem_DUSP2','CD8_CTL',
                                                 'CD8_Tex'],
                                    ordered=True)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=10, use_rep='X_pca', method='gauss')
sc.tl.paga(adata, groups='Celltype3')
with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.paga(adata, color=['Celltype3'], threshold=0.2, show=False, fontsize=14, frameon=True, node_size_scale=0.5, edge_width_scale=0.5,save=True)

############### pyscenic
import loompy as lp
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import glob

# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = "/CD8_scanpy/pyscenic_out/CD8_scenic.loom"

# path to anndata object, which will be updated to store Scanpy results as they are generated below
f_anndata_path = "/CD8_scanpy/pyscenic_out/CD8.h5ad"

# path to pyscenic output
f_pyscenic_output = "/CD8_scanpy/pyscenic_out/CD8_pyscenic_output.loom"

# loom output, generated from a combination of Scanpy and pySCENIC results:
f_final_loom = '/CD8_scanpy/pyscenic_out/CD8_scenic_integrated-output.loom'
row_attrs = { 
    "Gene": np.array(adata.var.index) ,
}
col_attrs = { 
    "CellID":  np.array(adata.obs.index) ,
    "nGene": np.array( adata.obs['nFeature_RNA']) ,
    "nUMI": np.array(  adata.obs['nCount_RNA']) ,
}

lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs )

# transcription factors list
f_tfs = "/pyscenic_file/allTFs_hg38.txt" # human
# tf_names = load_tf_names( f_tfs )
!pyscenic grn {f_loom_path_scenic} {f_tfs} -o adj.csv --num_workers 30

adjacencies = pd.read_csv("adj.csv", index_col=False, sep='\,')

# ranking databases
f_db_glob = "/data4/laiwp/pyscenic_file/*feather"
f_db_names = ' '.join( glob.glob(f_db_glob) )

# motif databases
f_motif_path = "/pyscenic_file/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

!pyscenic ctx adj.csv \
    {f_db_names} \
    --annotations_fname {f_motif_path} \
    --expression_mtx_fname {f_loom_path_scenic} \
    --output reg.csv \
    --mask_dropouts \
    --num_workers 30

!pyscenic aucell \
    {f_loom_path_scenic} \
    reg.csv \
    --output {f_pyscenic_output} \
    --num_workers 30

lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = lf.ra.Regulons

rss_cellType = regulon_specificity_scores( auc_mtx, adata.obs['Celltype3'] )

########################## scirpy
import scirpy as ir
import muon as mu

def load_scTCR(path, batch, type='10X'):
    assert type in ['10X', 'tracer', 'BD', 'h5ad'], 'type of scTCR-seq data must in 10X, tracer, BD, or h5ad.'
    adatas_tcr = {}
    
    for sample, sample_meta in zip(path, batch):
        if type == '10X':
            adata_tcr = ir.io.read_10x_vdj(sample)
        elif type == 'tracer':
            adata_tcr = ir.io.read_tracer(sample)
        elif type == 'BD':
            adata_tcr = ir.io.read_bd_rhapsody(sample)
        elif type == 'h5ad':
            adata_tcr = ir.io.read_h5ad(sample)
        adata_tcr.obs['batch'] = sample_meta
        adata_tcr.obs.index = sample_meta+'_'+adata_tcr.obs.index
        adatas_tcr[sample_meta] = adata_tcr

    adata_tcr = anndata.concat(adatas_tcr, index_unique=None)
    return adata_tcr

adata = sc.read_h5ad('TCR_counts.h5ad')
adata = adata[adata.obs.Group3.isin(['5R'])]
adata_tcr = load_scTCR(path=path_list, batch=['F012','F013','F014','P4','P5','P6'])
intersect_index = list(set(adata_tcr.obs.index).intersection(set(adata.obs.index)))
adata_tcr = adata_tcr[intersect_index,]
adata = adata[intersect_index,]
mdata = mu.MuData({"gex": adata, "airr": adata_tcr})

ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)

_ = ir.pl.group_abundance(mdata, groupby="airr:receptor_subtype", target_col="gex:Batch")

_ = ir.pl.group_abundance(mdata, groupby="airr:chain_pairing", target_col="gex:Batch")

mu.pp.filter_obs(mdata, "airr:chain_pairing", lambda x: ~np.isin(x, ["orphan VDJ", "orphan VJ"]))

ir.pp.ir_dist(mdata)
ir.tl.define_clonotypes(mdata, receptor_arms="all", dual_ir="primary_only")

ir.tl.clonotype_network(mdata, min_cells=2)
_ = ir.pl.clonotype_network(mdata, color="gex:Batch", base_size=5, label_fontsize=6, panel_size=(7, 7))

_ = ir.pl.clonal_expansion(mdata, target_col="clone_id", groupby="gex:Celltype3", breakpoints=(1, 2, 5), normalize=False,fig_kws = {'figsize': (8, 3)})
_ = ir.pl.clonal_expansion(mdata, target_col="clone_id", groupby="gex:Celltype3",breakpoints=(1, 2, 5),fig_kws = {'figsize': (8, 3)})
_ = ir.pl.alpha_diversity(mdata, metric="normalized_shannon_entropy", groupby="gex:Celltype3",fig_kws = {'figsize': (8, 3)})
_ = ir.pl.alpha_diversity(mdata, metric="normalized_shannon_entropy", groupby="gex:Batch",fig_kws = {'figsize': (3, 3)})
_ = ir.pl.group_abundance(
    mdata,
    groupby="airr:clone_id",
    target_col="gex:Celltype3",
    max_cols=10,
    normalize='gex:Batch',
    linewidth=0
)

_ = ir.pl.vdj_usage(
    mdata,
    full_combination=False,
    max_segments=None,
    max_ribbons=15,
    fig_kws={"figsize": (8, 5)},
)

df, dst, lk = ir.tl.repertoire_overlap(mdata, "gex:Celltype3", inplace=False)
ir.pl.repertoire_overlap(
    mdata,
    "gex:Celltype3",
    yticklabels=True,
    xticklabels=True,
)

vdjdb = ir.datasets.vdjdb()
ir.pp.ir_dist(mdata, vdjdb, metric="identity", sequence="aa")
ir.tl.ir_query(
    mdata,
    vdjdb,
    metric="identity",
    sequence="aa",
    receptor_arms="any",
    dual_ir="any",
)
ir.tl.ir_query_annotate(
    mdata,
    vdjdb,
    metric="identity",
    sequence="aa",
    include_ref_cols=["antigen.species"],
    strategy="most-frequent",
)






