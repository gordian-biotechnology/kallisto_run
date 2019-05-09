import subprocess
import sys
import os
from .kallisto_parser import get_parser
from .build_kallisto_idx import build_kallisto
from .kallisto_bus import *
from .genes_from_ens_kallisto import *
from .scvi_scanpy_anlysis import *
import scanpy as sc
import numpy as np
import pickle

def main():
    args = get_parser()
    kallisto_idx = ''
    genome = []
    annotation = []
    if os.path.exists(os.path.abspath(args.ref_dir)):
        for f in os.listdir(args.ref_dir):
            if f[-4:] == '.idx':
                print("Kallisto index file found:")
                print(f)
                kallisto_idx = os.path.join(os.path.abspath(args.ref_dir), f)


    if not os.path.exists(os.path.abspath(args.ref_dir)):
        os.mkdir(os.path.abspath(args.ref_dir))
    path_to_ref, path_to_anno, tr2g_path = make_kallisto_refs(os.path.abspath(args.ref_dir), args)

    if not kallisto_idx:
        build_kallisto(path_to_ref)


    if os.path.exists(args.out):
        if args.rerun:
            print('Overwriting output folder.')
        else:
            outpath = os.path.abspath(args.out)
    else:
        outpath = os.path.abspath(args.out)
        os.mkdir(outpath)
    if not os.path.exists(os.path.join(outpath, 'output.bus')):
        print('Running kallisto bus')
        run_kallisto_bus(args.input, kallisto_idx, outpath, args)
    elif args.rerun:
            print('Rerun: Running kallisto bus')
            run_kallisto_bus(args.input, kallisto_idx, outpath, args)
    if not os.path.exists(os.path.join(outpath, 'output.sort.txt')):
        print('Running bustools')
        run_bustools(outpath, args)
    elif args.rerun:
            print('Rerun: Running bustools')
            run_bustools(outpath, args)

    readENS_ids(path_to_anno, tr2g_path, outpath)

    if not all([os.path.exists(os.path.join(outpath, 'matrix.mtx')), os.path.exists(os.path.join(outpath, 'genes.tsv')),os.path.exists(os.path.join(outpath, 'barcodes.tsv'))]):
        matrix_from_kallisto_bus(outpath, tr2g_path)
    if not args.name:
        name = 'scanpy'
    else:
        name= args.name
    results_file = os.path.join(outpath, name+'_raw.h5ad')
    if not os.path.exists(results_file):
        adata = anndata_from_mtx(outpath, args.name)
    else:
        adata = sc.read(results_file)
    adata_filtered_path = os.path.join(outpath,'filtered_adata.h5ad')
    if not os.path.exists(adata_filtered_path):
        adata_filtered = run_analysis(adata, outpath, name, args.species, tr2g_path)
    else:
        adata_filtered = sc.read(adata_filtered_path)
        sc.settings.set_figure_params(dpi=200)
        figdir = os.path.join(outpath, 'figures')
        if not os.path.exists(figdir):
            os.mkdir(figdir)
        sc.settings.figdir= figdir
    scvi_posterior, scvi_latent = compute_scvi_latent(os.path.join(outpath,'original_for_scvi.h5ad'), n_epochs=50, n_latent=6)
    adata_filtered.obsm["X_scvi"] = scvi_latent
    sc.tl.pca(adata_filtered, svd_solver="arpack")
    sc.pp.neighbors(adata_filtered, n_neighbors=20, n_pcs=40, use_rep="X_scvi")
    sc.tl.umap(adata_filtered)
    sc.tl.louvain(adata_filtered, key_added="louvain_scvi", resolution=0.7)
    sc.pl.umap(adata_filtered, color="louvain_scvi", save='_louvain.png')
    n_genes = 20
    rank_genes_groups_bayes(adata_filtered, scvi_posterior, label_name="louvain_scvi", n_genes=n_genes)
    sc.pl.rank_genes_groups(adata_filtered, key="rank_genes_groups_scvi", sharey=False, n_genes=n_genes, save='_bayes_scvi.png')
    all_genes = len(adata_filtered.var_names)
    adata_filtered.write_h5ad(os.path.join(outpath,'filtered_adata_post_scvi.h5ad'))
    sc.tl.rank_genes_groups(adata_filtered, 'louvain_scvi', method='t-test',   use_raw=False, key_added='rank_genes_groups_ttest',  n_genes=all_genes)
    sc.tl.rank_genes_groups(adata_filtered, 'louvain_scvi', method='wilcoxon', use_raw=False, key_added='rank_genes_groups_wilcox', n_genes=all_genes)
    differential_expression = rank_genes_groups_bayes(adata_filtered, scvi_posterior, label_name='louvain_scvi', n_genes=all_genes)
    de_table = store_de_scores(adata_filtered, differential_expression, save_path=os.path.join(outpath, 'de_genes_bayes.csv'))

if __name__ == '__main__':
    main()
