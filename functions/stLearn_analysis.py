import stlearn as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

# Loading raw data #
p0928s1_data_dir = "./20221129WF3"
p0928s1_data = st.Read10X(p0928s1_data_dir)
p0928s1_data.var_names_make_unique()
HE_path="./spatial/tissue_hires_image.png"
st.add.image(adata=p0928s1_data,
             imgpath=HE_path,
             library_id="20211129WF3", visium=True)

# Basic normalisation #
st.pp.filter_genes(p0928s1_data, min_cells=3)
st.pp.normalize_total(p0928s1_data)
p0928s1_spot_mixtures = pd.read_csv("p0928s1_label_transfer_bc.csv", 
                                    index_col=0)
p0928s1_labels = p0928s1_spot_mixtures.loc[:,'spotlight.cluster'].values.astype(str)
p0928s1_spot_mixtures = p0928s1_spot_mixtures.drop(['spotlight.cluster'],
                                                   axis=1)

# NOTE: using the same key in p0928s1_data.obs & p0928s1_data.uns
p0928s1_data.obs['cell_type'] = p0928s1_labels # Adding the dominant cell type labels per spot
p0928s1_data.obs['cell_type'] = p0928s1_data.obs['cell_type'].astype('category')
p0928s1_data.uns['cell_type'] = p0928s1_spot_mixtures # Adding the cell type scores
st.pl.cluster_plot(p0928s1_data, use_label='cell_type',figsize=(12,8),size=35)
plt.savefig("stlearn_p0928s1_celltype.png")
plt.show()
# Loading the LR databases available within stlearn (from NATMI)
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
# Running the analysis #
st.tl.cci.run(p0928s1_data, lrs,
              min_spots = 20, #Filter out any LR pairs with no scores for less than min_spots
              distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
              n_pairs=1000, # Number of random pairs to generate; low as example, recommend ~10,000
              n_cpus=128, # Number of CPUs for parallel. If None, detects & use all available.
)
st.tl.cci.adj_pvals(p0928s1_data, correct_axis='spot', pval_adj_cutoff=0.05, adj_method='fdr_bh')
st.tl.cci.run_cci(p0928s1_data, 'cell_type', min_spots=3, spot_mixtures=True, cell_prop_cutoff=0.2, sig_spots=True, n_perms=100)

#CCI 
# Visualising the no. of interactions between cell types across all LR pairs #
pos_1 = st.pl.ccinet_plot(p0928s1_data, 'cell_type', return_pos=True)
plt.savefig('./p0928s1.ccinet_plot.png')
plt.close()

lrs = p0928s1_data.uns['lr_summary'].index.values[0:20]
st.pl.lr_cci_map(p0928s1_data, 'cell_type', lrs=lrs, min_total=100, figsize=(20,12),square_scaler=5000)
plt.savefig('./p0928s1.lr_cci_map.png')
plt.close()

p0928s1_data.filename = './p0928s1_LR.h5ad'
print(p0928s1_data.isbacked) # True
