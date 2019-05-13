# Seed for reproducability
import torch
import numpy as np
import pandas as pd
import scanpy as sc
import os
from scipy import sparse

def run_analysis(adata, outpath, name, species, tr2g_path):
    sc.settings.verbosity = 0  # verbosity: errors (0), warnings (1), info (2), hints (3)
    torch.manual_seed(0)
    np.random.seed(0)
    sc.settings.set_figure_params(dpi=200)
    figdir = os.path.join(outpath, 'figures')
    if not os.path.exists(figdir):
        os.mkdir(figdir)
    sc.settings.figdir =figdir
    min_genes = 200
    min_cells = 3
    if adata.var_names[0][0:3] == 'ENS':
        gene_names = {}
        with open(tr2g_path) as f:
            f.readline()
            for line in f:
                t,g,gn = line.split()
                gene_names[g] = gn
        gnames = [gene_names[a] for a in adata.var_names]
        adata.var['gene_ids'] = adata.var.index
        adata.var.index = gnames
        adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_genes=1)
    adata = adata[adata.obs['doublet_score']<adata.uns['doublet_threshold']]
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    if species == 'mouse':
        mito_genes = adata.var_names.str.startswith("mt-")
    else:
        mito_genes = adata.var_names.str.startswith("MT-")
    num_mito_found = len(adata.var[mito_genes])
    if num_mito_found> 10:
        adata.obs["percent_mito"] = (np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1)
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save='_prefilter.png')
        adata = adata[adata.obs["percent_mito"] < np.percentile(np.array(adata.obs["percent_mito"]), 99), :]
    else:
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True, save='_prefilter.png')
    adata = adata[adata.obs["n_genes_by_counts"] < np.percentile(np.array(adata.obs["n_genes_by_counts"]), 99), :]

    if num_mito_found > 10:
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save='_postfilter.png')
    else:
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True, save='_postfilter.png')
    adata_original = adata.copy()

    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)


    min_mean = 0.0125
    max_mean = 3
    min_disp = 0.5
    max_disp = np.inf

    sc.pp.highly_variable_genes(
        adata,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        max_disp=max_disp
        # n_top_genes=500
    )

    adata.raw = adata

    highly_variable_genes = adata.var["highly_variable"]
    adata = adata[:, highly_variable_genes]
    print(highly_variable_genes.sum(), 'genes most variable.')
    if num_mito_found> 10:
        sc.pp.regress_out(adata, ["n_counts", "percent_mito"])
    else:
        sc.pp.regress_out(adata, ["n_counts"])
    sc.pp.scale(adata, max_value=10)

    # Also filter the original adata genes
    adata_original = adata_original[:, highly_variable_genes]

    # We also store adata_original into adata.raw
    # (which was designed for this purpose but actually has limited functionnalities)
    adata.raw = adata_original
    adata_original.write_h5ad(os.path.join(outpath,'original_for_scvi.h5ad'))
    adata.write_h5ad(os.path.join(outpath,'filtered_adata.h5ad'))
    return adata

import scvi
from scvi.dataset.anndata import AnnDataset
from scvi.inference import UnsupervisedTrainer
from scvi.models.vae import VAE
from typing import Tuple


def compute_scvi_latent(
    adata: sc.AnnData,
    n_latent: int = 5,
    n_epochs: int = 100,
    lr: float = 1e-3,
    use_batches: bool = False,
    use_cuda: bool = False,
) -> Tuple[scvi.inference.Posterior, np.ndarray]:
    """Train and return a scVI model and sample a latent space

    :param adata: sc.AnnData object non-normalized
    :param n_latent: dimension of the latent space
    :param n_epochs: number of training epochs
    :param lr: learning rate
    :param use_batches
    :param use_cuda
    :return: (scvi.Posterior, latent_space)
    """
    # Convert easily to scvi dataset
    scviDataset = AnnDataset(adata)

    # Train a model
    vae = VAE(
        scviDataset.nb_genes,
        n_batch=scviDataset.n_batches * use_batches,
        n_latent=n_latent,
    )
    trainer = UnsupervisedTrainer(vae, scviDataset, train_size=1.0, use_cuda=use_cuda)
    trainer.train(n_epochs=n_epochs, lr=lr)
    ####

    # Extract latent space
    posterior = trainer.create_posterior(
        trainer.model, scviDataset, indices=np.arange(len(scviDataset))
    ).sequential()

    latent, _, _ = posterior.get_latent()

    return posterior, latent

def rank_genes_groups_bayes(
    adata: sc.AnnData,
    scvi_posterior: scvi.inference.Posterior,
    n_samples: int = None,
    M_permutation: int = None,
    n_genes: int = 25,
    label_name: str = "louvain_scvi",
) -> pd.DataFrame:
    """
    Rank genes for characterizing groups.
    Computes Bayes factor for each cluster against the others to test for differential expression.
    See Nature article (https://rdcu.be/bdHYQ)

    :param adata: sc.AnnData object non-normalized
    :param scvi_posterior:
    :param n_samples:
    :param M_permutation:
    :param n_genes:
    :param label_name: The groups tested are taken from adata.obs[label_name] which can be computed
                       using clustering like Louvain (Ex: sc.tl.louvain(adata, key_added=label_name) )
    :return: Summary of Bayes factor per gene, per cluster
    """

    # Call scvi function
    per_cluster_de, cluster_id = scvi_posterior.one_vs_all_degenes(
        cell_labels=np.asarray(adata.obs[label_name].values).astype(int).ravel(),
        min_cells=1,
        n_samples=n_samples,
        M_permutation=M_permutation,
    )

    # convert to ScanPy format -- this is just about feeding scvi results into a format readable by ScanPy
    markers = []
    scores = []
    names = []
    for i, x in enumerate(per_cluster_de):
        subset_de = x[:n_genes]
        markers.append(subset_de)
        scores.append(tuple(subset_de["bayes1"].values))
        names.append(tuple(subset_de.index.values))

    markers = pd.concat(markers)
    dtypes_scores = [(str(i), "<f4") for i in range(len(scores))]
    dtypes_names = [(str(i), "<U50") for i in range(len(names))]
    scores = np.array([tuple(row) for row in np.array(scores).T], dtype=dtypes_scores)
    scores = scores.view(np.recarray)
    names = np.array([tuple(row) for row in np.array(names).T], dtype=dtypes_names)
    names = names.view(np.recarray)

    adata.uns["rank_genes_groups_scvi"] = {
        "params": {
            "groupby": "",
            "reference": "rest",
            "method": "",
            "use_raw": True,
            "corr_method": "",
        },
        "scores": scores,
        "names": names,
    }
    return markers

def store_de_scores(
    adata: sc.AnnData, differential_expression: pd.DataFrame, save_path: str = None
):
    """Creates, returns and writes a DataFrame with all the differential scores used in this notebook.

    Args:
        adata: scRNAseq dataset
        differential_expression: Pandas Dataframe containing the bayes factor for all genes and clusters
        save_path: file path for writing the resulting table

    Returns:
        pandas.DataFrame containing the scores of each differential expression test.

    """
    # get shapes for array initialisation
    n_genes_de = differential_expression[
        differential_expression["clusters"] == 0
    ].shape[0]
    all_genes = adata.shape[1]
    # check that all genes have been used
    if n_genes_de != all_genes:
        raise ValueError(
            "scvi differential expression has to have been run with n_genes=all_genes"
        )
    # get tests results from AnnData unstructured annotations
    rec_scores = []
    rec_names = []
    test_types = ["ttest", "wilcox"]
    for test_type in test_types:
        res = adata.uns["rank_genes_groups_" + test_type]
        rec_scores.append(res["scores"])
        rec_names.append(res["names"])
    # restrict scvi table to bayes factor
    res = differential_expression[["bayes1", "clusters"]]
    # for each cluster join then append all
    dfs_cluster = []
    groups = res.groupby("clusters")
    for cluster, df in groups:
        for rec_score, rec_name, test_type in zip(rec_scores, rec_names, test_types):
            temp = pd.DataFrame(
                rec_score[str(cluster)],
                index=rec_name[str(cluster)],
                columns=[test_type],
            )
            df = df.join(temp)
        dfs_cluster.append(df)
    res = pd.concat(dfs_cluster)
    if save_path:
        res.to_csv(save_path)
    return res
