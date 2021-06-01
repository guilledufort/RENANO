import subprocess
import operator as op
import numpy as np
import csv
import operator
import sys
import os
from IPython import embed

import argparse
parser = argparse.ArgumentParser()

#-d DATASET -r REPORT_NAME -c MINCOVERAGE -m MAXSPECIES
parser.add_argument("-d", "--dataset", help="Dataset name")
parser.add_argument("-r", "--report", help="Report file name")
parser.add_argument("-c", "--coverage", help="Min coverage", type=float)
parser.add_argument("-m", "--maxspecies", help="Max num. species", type=int)

args = parser.parse_args()

db_name = args.dataset
print("Dataset: ", db_name)

f_name = args.report
if f_name == None:
    f_name = "datasets/{}/{}.report".format(db_name, db_name)

print("Report file: ", f_name)

min_cov = args.coverage
if min_cov==None:
    min_cov = 0.3

print("Minimum species coverage: ", str(min_cov))

max_species = args.maxspecies
if max_species==None:
    max_species = 20

print("Maximum number of species to add to the ref: ", str(max_species))

total_comp = 0
total = 0

species = []

with open(f_name) as csvfile:
    readCSV = csv.reader(csvfile, delimiter='\t')
    for row in readCSV:
        if row[3] == 'S':
            species.append(row)

species = sorted(species, key = lambda x: float(x[0]))[::-1]

max_refs_per_species = 1
max_fails = 5
    
dir_name = "aux_refs"
os.system("mkdir -p " + dir_name)

cov = 0
j = 0

print("Downloading reference genomes from NCBI database:")
for s in species:
    s_cov = float(s[0])
    if j == max_species or s_cov < min_cov:
        break;
    accs_file = dir_name + "/accs.txt"
    spec = s[5].lstrip()
    print(spec)
    cmd = "esearch -db assembly -query \'\"{}\"[Organism] AND latest[filter] AND (all[filter] NOT anomalous[filter] AND all[filter] NOT \"derived from surveillance project\"[filter])\' | esummary | xtract -pattern DocumentSummary -element AssemblyAccession > {}".format(spec, accs_file)
    os.system(cmd)
    with open(accs_file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter='\t')
        i = 0
        fails = 0
        for row in readCSV:
            if i < max_refs_per_species and fails < max_fails:
                cmd = "esearch -db assembly -query {} | elink -target nucleotide -name assembly_nuccore_refseq | efetch -format fasta > {}/{}.fa".format(row[0], dir_name, spec.replace(" ", "") + "_" + str(i))
                output = subprocess.getoutput(cmd)
                if (output == ''):
                    fails = 0
                    i += 1
                    j += 1
                    cov += s_cov
                else:
                    fails += 1
            else:
                if (fails == max_fails):
                    print("Failed to download NCBI reference in {} tries".format(max_fails))
                break;
    
    
    os.system("rm -rf " + accs_file)

print("The total perecentage of reads covered by these species is {}\%".format(cov))

os.system("cat {}/*.fa > datasets/{}/{}_genome.fna".format(dir_name, db_name, db_name))

os.system("rm -rf " + dir_name)