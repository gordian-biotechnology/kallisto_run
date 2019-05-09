import os
import sys
import subprocess

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
        print(r1_files, r2_files)
        fq_zipped = zip(r1_files, r2_files)
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
