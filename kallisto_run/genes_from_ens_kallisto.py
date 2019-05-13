import numpy as np
import pickle
import mygene
import gzip
import subprocess
import os


def create_transcript_list(input, use_name = True, use_version = True):
    r = {}
    for line in input:
        if len(line) == 0 or line[0] == '#':
            continue
        l = line.strip().split('\t')
        if l[2] == 'transcript':
            info = l[8]
            d = {}
            for x in info.split('; '):
                x = x.strip()
                p = x.find(' ')
                if p == -1:
                    continue
                k = x[:p]
                p = x.find('"',p)
                p2 = x.find('"',p+1)
                v = x[p+1:p2]
                d[k] = v


            if 'transcript_id' not in d or 'gene_id' not in d:
                continue

            tid = d['transcript_id']
            gid = d['gene_id']
            if use_version:
                if 'transcript_version' not in d or 'gene_version' not in d:
                    continue

                tid += '.' + d['transcript_version']
                gid += '.' + d['gene_version']
            gname = None
            if use_name:
                if 'gene_name' not in d:
                    continue
                gname = d['gene_name']

            if tid in r:
                continue

            r[tid] = (gid, gname)
    return r

def ec2g(ec, tr2g, ecs):
    if ec in ecs:
        return list(set(tr2g[trlist[t]] for t in ecs[ec]))
    else:
        return []

def print_output(output, r, use_name = True):
    for tid in r:
        if use_name:
            output.write("%s\t%s\t%s\n"%(tid, r[tid][0], r[tid][1]))
        else:
            output.write("%s\t%s\n"%(tid, r[tid][0]))

def readENS_ids(path_to_gtf, t2g_path, outpath):
    if path_to_gtf[-3:] == '.gz':
        if os.path.exists(path_to_gtf):
            subprocess.call('gunzip -v %s' %(path_to_gtf), shell=True)
        path_to_gtf = path_to_gtf[0:-3]
    with open(path_to_gtf) as file:
        r = create_transcript_list(file, use_name = True, use_version = True)
    with open(t2g_path, "w+") as output:
        print_output(output, r, use_name = True)



def make_kallisto_refs(ref_path, args):
    species = args.species
    run_wget = False
    if species.lower() == 'human':
        ref='GRCh38'
        ens_='ENS'
        t2g_path = os.path.join(ref_path, 'human_transcript_to_gene.tsv')
        path_to_ref = os.path.join(ref_path, 'Homo_sapiens.GRCh38.cdna.all.fa.gz')
        path_to_anno = os.path.join(ref_path, 'Homo_sapiens.GRCh38.94.gtf.gz')
        if not os.path.exists(path_to_ref) and not os.path.exists(path_to_ref[0:-3]):
            wget_cmd = "wget -O %s ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz" %(path_to_ref)
            run_wget = True
        if not os.path.exists(path_to_anno) and not os.path.exists(path_to_anno[0:-3]):
            wget_cmd += '&& wget -O %s ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz' %(path_to_anno)
    elif species.lower() == 'mouse':
        ref='Mus_musculus.GRCm38'
        ens_='ENSMUS'
        t2g_path = os.path.join(ref_path, 'mouse_transcript_to_gene.tsv')
        path_to_ref = os.path.join(ref_path, 'Mus_musculus.GRCm38.cdna.all.fa.gz')
        path_to_anno = os.path.join(ref_path, 'Mus_musculus.GRCm38.94.gtf.gz')
        if not os.path.exists(path_to_ref) and not os.path.exists(path_to_ref[0:-3]):
            wget_cmd = "wget -O %s ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz" %(path_to_ref)
            run_wget = True
        if not os.path.exists(path_to_anno) and not os.path.exists(path_to_anno[0:-3]):
            wget_cmd += '&& wget -O %s ftp://ftp.ensembl.org/pub/release-94/gtf/mus_musculus/Mus_musculus.GRCm38.94.gtf.gz' %(path_to_anno)
    elif species.lower() == 'mouse_human':
        ref='Mus_human.GRCm38'
        ens_='ENSMUS'
        t2g_path = os.path.join(ref_path, 'human_mouse_transcript_to_gene.tsv')
        path_to_ref1 = os.path.join(ref_path, 'Mus_musculus.GRCm38.cdna.all.fa.gz')
        path_to_ref2 = os.path.join(ref_path, 'Homo_sapiens.GRCh38.cdna.all.fa.gz')
        path_to_anno1 = os.path.join(ref_path, 'Mus_musculus.GRCm38.cdna.all.fa.gz')
        path_to_anno2 = os.path.join(ref_path, 'Homo_sapiens.GRCh38.94.gtf.gz')
        path_to_ref  = os.path.join(ref_path, 'human_mouse_contatenated_transcriptome.fa')
        path_to_anno = os.path.join(ref_path, 'human_mouse_contatenated_GTF.gtf')
        wget_cmd = ''
        if not os.path.exists(path_to_ref1) and not os.path.exists(path_to_ref1[0:-3]):
            wget_cmd = "wget -O %s ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz" %(path_to_ref)
            run_wget = True
        if not os.path.exists(path_to_ref2) and not os.path.exists(path_to_ref2[0:-3]):
            if wget_cmd:
                wget_cmd += "&& wget -O %s ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz" %(path_to_ref)
            else:
                wget_cmd = "wget -O %s ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz" %(path_to_ref)
            run_wget = True
        if not os.path.exists(path_to_anno1) and not os.path.exists(path_to_anno1[0:-3]):
            wget_cmd += '&& wget -O %s ftp://ftp.ensembl.org/pub/release-94/gtf/mus_musculus/Mus_musculus.GRCm38.94.gtf.gz' %(path_to_anno1)
        if not os.path.exists(path_to_anno2) and not os.path.exists(path_to_anno2[0:-3]):
            wget_cmd += '&& wget -O %s ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz' %(path_to_anno2)
        if not os.path.exists(path_to_ref):
            wget_cmd += 'gunzip -v %s && gunzip -v %s && cat %s %s > %s' %(path_to_ref1,  path_to_ref2, path_to_ref1.strip('.gz'), path_to_ref2.strip('.gz'), os.path.join(ref_path, 'human_mouse_contatenated_transcriptome.fa'))
        if not os.path.exists(path_to_anno):
            wget_cmd += 'gunzip -v %s && gunzip -v %s && cat %s %s > %s' %(path_to_anno1,  path_to_anno2, path_to_anno1.strip('.gz'), path_to_anno2.strip('.gz'), os.path.join(ref_path, 'human_mouse_contatenated_GTF.gtf'))
    if not args.dry_run:
        if run_wget:
            process = subprocess.Popen(wget_cmd, stdout=subprocess.PIPE, shell=True)
            out, err = process.communicate()
            print(out)

    return path_to_ref, path_to_anno, t2g_path

def matrix_from_kallisto_bus(outpath, tr2g_path):
    import csv
    from collections import defaultdict
    import collections
    #load transcript to gene file
    tr2g = {}
    trlist = []
    genes = []
    gene_names ={}
    with open(tr2g_path) as f:
        f.readline()
        for line in f:
            t,g,gn = line.split()
            gene_names[g] = gn
            tr2g[t] = g
            trlist.append(t)
            genes.append(g)

    genes = list(set(genes))

    # load equivalence classes
    ecs = {}
    with open(os.path.join(outpath,'matrix.ec')) as f:
        for line in f:
            l = line.split()
            ec = int(l[0])
            trs = [int(x) for x in l[1].split(',')]
            ecs[ec] = trs

    def ec2g(ec):
        if ec in ecs:
            return list(set(tr2g[trlist[t]] for t in ecs[ec]))
        else:
            return []

    cell_gene = collections.defaultdict(lambda: collections.defaultdict(float))
    pbar=None
    pumi=None
    with open(os.path.join(outpath,'output.sort.txt')) as f:
        gs = set()
        for line in f:
            l = line.split()
            barcode,umi,ec,count = line.split()
            ec = int(ec)

            if barcode == pbar:
                # same barcode
                if umi == pumi:
                    # same UMI, let's update with intersection of genelist
                    gl = ec2g(ec)
                    gs.intersection_update(gl)
                else:
                    # new UMI, process the previous gene set
                    for g in gs:
                        cell_gene[barcode][g] += 1.0/len(gs)
                    # record new umi, reset gene set
                    pumi = umi
                    gs = set(ec2g(ec))
            else:
                # work with previous gene list
                for g in gs:
                    cell_gene[pbar][g] += 1.0/len(gs)

                if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
                    del cell_gene[pbar]

                pbar = barcode
                pumi = umi

                gs = set(ec2g(ec))
        #remember the last gene
        for g in gs:
            cell_gene[pbar][g] += 1.0/len(gs)

        if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
            del cell_gene[pbar]

    barcode_hist = collections.defaultdict(int)
    for barcode in cell_gene:
        cg = cell_gene[barcode]
        s = len([cg[g] for g in cg])
        barcode_hist[barcode] += s
    bad_barcode = [x for x in barcode_hist if  barcode_hist[x] <= 20]
    print("Percent barcodes not used: "+str(len(bad_barcode)/len(cell_gene)))

    s = 0
    bad_s = 0
    bad_barcode_set = set(bad_barcode)
    for barcode in cell_gene:
        cg = cell_gene[barcode]
        cgs =  sum(cg[g] for g in cg)
        s += cgs
        if barcode in bad_barcode_set:
            bad_s += cgs

    print('Percent cells bad: '+str(bad_s/s))

    outfile = os.path.join(outpath,'matrix.mtx')

    gene_to_id = dict((g,i+1) for i,g in enumerate(genes))
    barcodes_to_use = [b for b,x in barcode_hist.items() if x > 500 and x < 12000]

    num_entries = 0
    for barcode in barcodes_to_use:
        num_entries += len([x for x in cell_gene[barcode].values() if round(x)>0])
    with open(outfile, 'w') as of:
        of.write('%%MatrixMarket matrix coordinate real general\n%\n')
        #number of genes
        of.write("%d %d %d\n"%(len(genes), len(barcodes_to_use), num_entries))
        bcid = 0
        for barcode in barcodes_to_use:
            bcid += 1
            cg = cell_gene[barcode]
            gl = [(gene_to_id[g],round(cg[g])) for g in cg if round(cg[g]) > 0]
            gl.sort()
            for x in gl:
                of.write("%d %d %d\n"%(x[0],bcid,x[1]))

    gl = []
    for g in genes:
        gl.append((g,gene_names[g]))

    with open(os.path.join(outpath, 'genes.tsv'),'w') as of:
        for g,gn in gl:
            of.write("%s\t%s\n"%(g,gn))

    with open(os.path.join(outpath, 'barcodes.tsv'),'w') as of:
        of.write('\n'.join(x + '-1' for x in barcodes_to_use))
        of.write('\n')

def anndata_from_mtx(outpath, name):
    import numpy as np
    import pandas as pd
    import scanpy.api as sc
    import scrublet as scr
    from scipy import sparse
    from skimage.filters import threshold_minimum

    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_versions()
    if not name:
        name = 'scanpy'
    results_file = os.path.join(outpath, name+'_raw.h5ad')
    sc.settings.set_figure_params(dpi=80)

    adata = sc.read(os.path.join(outpath,'matrix.mtx'), cache=False).T  # transpose the data
    adata.var_names = pd.read_csv(os.path.join(outpath,'genes.tsv'), header=None, sep='\t')[0]
    adata.obs_names = pd.read_csv(os.path.join(outpath,'barcodes.tsv'), header=None, sep='\t')[0]
    adata.var_names_make_unique()
    counts_matrix = sparse.csc_matrix(adata.X)

    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=round(counts_matrix.shape[0]/125000, 4))
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                          min_cells=3,
                                                          min_gene_variability_pctl=85,
                                                          n_prin_comps=30)
    threshold = threshold_minimum(scrub.doublet_scores_sim_)
    adata.obs['doublet_score'] = scrub.doublet_scores_obs_
    adata.uns['doublet_threshold'] = threshold
    adata.write_h5ad(results_file)
    return adata
