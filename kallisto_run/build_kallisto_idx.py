import os
import subprocess

def build_kallisto(genome):
  genome_name = genome.split('.fa')[0]
  idx_name = os.path.join(os.path.dirname(genome),genome_name +'.idx')
  cmd = "kallisto index -i %s %s" %(idx_name, genome)
  process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
  out, err = process.communicate()
  print(out)
