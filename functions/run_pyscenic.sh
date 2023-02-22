## pyscenic shell ----
# pyscenic grn
pyscenic grn --num_workers 32 --output adj.fibroblast.tsv --method grnboost2 fibroblast.loom hs_hgnc_tfs.txt

#pyscenic cistarget
pyscenic ctx adj.fibroblast.tsv hg19-tss-centered-10kb-7species.mc9nr.feather --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname fibroblast.loom \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 32 \
--mask_dropouts

#pyscenic AUCell
pyscenic aucell fibroblast.loom reg.csv --output fibroblast_SCENIC.loom --num_workers 32
