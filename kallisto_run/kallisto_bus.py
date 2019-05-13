import os
import sys
import subprocess
import gzip
from Bio import SeqIO
from collections import Counter
import numpy as np

def check_r1_polyT(input, r1_files, r2_files):
    poly_t_r1 = []
    poly_t_r2 = []
    for r1 in r1_files:
        i= 0
        with gzip.open(os.path.join(input, r1), "rt") as handle1:
            for record1 in SeqIO.parse(handle1, "fastq"):
                if i <300:
                    i+=1
                    seq = record1.seq
                    poly_t_r1.append(Counter(seq)['T']/len(seq))
                else:
                    break
    for r2 in r2_files:
        i= 0
        with gzip.open(os.path.join(input, r2), "rt") as handle2:
            for record2 in SeqIO.parse(handle2, "fastq"):
                if i <300:
                    i+=1
                    seq = record2.seq
                    poly_t_r2.append(Counter(seq)['T']/len(seq))
                else:
                    break
    r2_t_mean = np.array(poly_t_r2).mean()
    r1_t_mean = np.array(poly_t_r1).mean()
    if  r2_t_mean > r1_t_mean:
        print('R2 files appear polyT enriched, switching order.')
        print('R2 mean: ', r2_t_mean)
        print('R1 mean: ', r1_t_mean)
        fq_zipped = zip(r2_files, r1_files)
        return fq_zipped
    else:
        print('R2 mean: ', r2_t_mean)
        print('R1 mean: ', r1_t_mean)
        fq_zipped = zip(r1_files, r2_files)
        return fq_zipped



def run_kallisto_bus(input, kallisto_idx, out, args):
    r1_files = []
    r2_files = []
    if os.path.exists(input):
        for fq in os.listdir(input):
            if '.fastq.gz' in fq:
                fq_split = fq.split('_')
                if fq_split[-2] =='R1' or fq_split[-1].split('.')[0] =='1':
                    r1_files.append(fq)
                elif fq_split[-2] =='R2' or fq_split[-1].split('.')[0] =='2':
                    r2_files.append(fq)
        r1_files.sort()
        r2_files.sort()
        fq_zipped  = check_r1_polyT(input, r1_files, r2_files)

        ordered_fastqs = []
        for (r1, r2) in fq_zipped:
            ordered_fastqs.append(os.path.join(input,r1))
            ordered_fastqs.append(os.path.join(input,r2))
        fastq_cmd_str = ' '.join(ordered_fastqs)
        kallisto_bs_cmd = 'kallisto bus -i %s -x %s -t 8 -o %s %s' %(kallisto_idx,args.version, out, fastq_cmd_str)
        if args.dry_run:
            print(kallisto_bs_cmd)
        else:
            process = subprocess.Popen(kallisto_bs_cmd, stdout=subprocess.PIPE, shell=True)
            out, err = process.communicate()
            print(out)

def run_bustools(input_dir, args):
    bus_in = os.path.join(input_dir,'output.bus')

    if not os.path.exists(bus_in) and not args.dry_run:
        sys.exit("Kallisto bus output doesn't exist at "+input_dir)
    out_sort = os.path.join(input_dir, 'output.sort.bus')
    out_txt = os.path.join(input_dir, 'output.sort.txt')
    bus_cmd = 'bustools sort -t 8 -o %s %s' %(out_sort, bus_in)
    txt_cmd = 'bustools text -o %s %s' %(out_txt, out_sort)
    if not args.dry_run:
        process = subprocess.Popen(bus_cmd, stdout=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        print(out)
    else:
        print(bus_cmd)
        print(txt_cmd)
    if os.path.exists(out_sort):
        if not args.dry_run:
            process = subprocess.Popen(txt_cmd, stdout=subprocess.PIPE, shell=True)
            out, err = process.communicate()
            print(out)
    elif not args.dry_run:
        sys.exit('bustolls sort failed to produce proper output.')
