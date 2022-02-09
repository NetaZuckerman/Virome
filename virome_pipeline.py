import os
import multiprocessing as mp
import sys
import time
import subprocess
import io
import logging
import argparse
from pathlib import Path
import csv
import pandas as pd
import gzip
from functools import partial
import json

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import (pyplot as plt)
from scipy.stats import entropy
import matplotlib
matplotlib.use('Agg')

def define_parser():
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '-i', '--input',
        dest='input',
        type=Path,
        default='None',
        help="Input path"
        )
    
    parser.add_argument(
        '-o', '--output',
        dest='output',
        type=Path,
        default='./',
        help="Output path"
        )
    
    parser.add_argument(
        '-t', '--threads',
        dest='threads',
        type=int,
        default=1,
        help="Number of threads"
        )
    
    parser.add_argument(
        '-c', '--ref_clean',
        dest='ref_clean',
        type=str,
        default='',
        help="references_to_clean"
        )
    
    parser.add_argument(
        '-r', '--viral_ref',
        dest='viral_ref',
        type=str,
        default='',
        help="viral reference"
        )

    return parser



def verify_path(args, arg_name):
    if arg_name not in args:
        raise ValueError(f"{arg_name} not provided.")
    _path = args.__dict__[arg_name]
    if not _path.exists():
        raise FileNotFoundError(f"{arg_name}: {_path} doesn't exist")
    return _path

#%% GLOBAL    



class MP():
    
    def __init__(self, num_threads):
        self.num_threads = num_threads
        
    def mltiprocess_command(self, func,list_paths,arguments):
        sema = mp.Semaphore(self.num_threads)
        pss = []
        for path in list_paths:
            p = mp.Process(target=func, args=(arguments,path,sema))
            p.start()
            pss.append(p)
      
        for p in pss:
            p.join()
        
    def set_threads(self,t):
        self.num_threads = t
        

        
#%% FASTQ FILTERING 

FILTER_FASTQ = 'bwa mem -v 1%(p_param)s -a -t 20 %(genome)s %(fastq)s %(fastq2)s | samtools view -b -F 2 - | samtools fastq - > %(path)s/%(name)s.fastq'

def filter_fastq_command(arguments,fastq,sema):
    former_key,path,genome = arguments
    sema.acquire()
    name = os.path.basename(fastq)
    if former_key:
        name = name.replace('.fastq','')
        fastq2 = ''
        p_param = ' -p'
    else:
        name = name.replace('_R1_001.fastq.gz','')
        fastq2 = fastq.replace('_R1','_R2')
        p_param = ''
        
    print(name)
    print('bwa mem -v 1%(p_param)s -a -t 20 %(genome)s %(fastq)s %(fastq2)s | samtools view -b -F 2 - | samtools fastq - > %(path)s/%(name)s.fastq' % dict(name=name,p_param=p_param,path=path,genome=genome,fastq=fastq,fastq2=fastq2))
    subprocess.call(FILTER_FASTQ % dict(name=name,p_param=p_param,path=path,genome=genome,fastq=fastq,fastq2=fastq2), shell=True)
    sema.release()

    
def filter_fastq(args,MP_obj, offset=30):
    path_list = [os.path.join(r, f) for r, d, fs in os.walk(args.input) for f in fs if f.endswith('_R1_001.fastq.gz') and 'singelton' not in f]
    path_new_dir = None
    
    references_to_clean = json.loads(args.ref_clean)
    former_key = None
    for key in list(references_to_clean.keys()): 
        dir_name = 'fastq/no' + key.split('_')[0]
        path_new_dir = os.path.join(args.output,dir_name)
        os.makedirs(path_new_dir,exist_ok=True)
            
        if former_key:
            path_list = [os.path.join(former_key, file) for file in os.listdir(former_key)]
        
        if 'bacteria' in references_to_clean.get(key):
            MP_obj.set_threads(1)
        MP_obj.mltiprocess_command(filter_fastq_command,path_list,[former_key,path_new_dir,references_to_clean.get(key)])
        if MP_obj.num_threads != args.threads:
            MP_obj.set_threads(args.threads)
        
        former_key = path_new_dir
    return(path_new_dir)
            

#%% FASTQ METASPADES SCAFFOLDS

SCAFFOLD = 'metaspades.py --12 %(fastq)s --threads 10 -o %(assembly_dir)s/%(name)s'

def scaffold_command(assembly_dir,fastq,sema):
    sema.acquire()
    name = os.path.basename(fastq)
    name = name.replace('.fastq','')  
    print(name)
    subprocess.call(SCAFFOLD % dict(name=name,fastq=fastq,assembly_dir=assembly_dir), shell=True)
    sema.release()

def scaffold(filtered_fastq,out,MP_obj):
    filtered_fastq = json.loads(filtered_fastq)
    filtered_fastq_dir = list(filtered_fastq)[-1]
    dir_name = 'fastq/no' + filtered_fastq_dir.split('_')[0]
    filtered_fastq_dir = os.path.join(out,dir_name)
    path_list = [os.path.join(r, f) for r, d, fs in os.walk(filtered_fastq_dir) for f in fs]
    assembly_path = os.path.join(out,'assembly')
    os.makedirs(assembly_path, exist_ok=True)
    MP_obj.mltiprocess_command(scaffold_command,path_list,assembly_path)
    return (assembly_path)
    
#%% MAPPING AND REPORT the contigs

# MAP_AND_REPORT = '''bwa mem -v 1 -t 10 -a %(genome)s %(assembly_dirs)s/scaffolds.fasta \
#  | samtools view -b - | tee %(contig_path)s/%(name)s.bam | samtools view -b -F 4 - \
#      | tee %(contig_path)s/%(name)s.mapped.bam | samtools view - | cut -f3 | sort | uniq -c \
#          | awk -v i=$i '{print i$1"\t"$2}' > %(report_file)s '''
MAP_AND_REPORT = '''bwa mem -v 1 -t 10 -a %(genome)s %(assembly_dirs)s/scaffolds.fasta \
  | samtools view -b - | tee %(contig_path)s/%(name)s.bam | samtools view -b -F 4 - \
      | tee %(contig_path)s/%(name)s.mapped.bam | samtools view -q 1 -F 260 - | cut -f3 | sort | uniq -c \
          | awk -v i=$i '{print i$1"\t"$2}' > %(report_file)s'''


def map_and_report_command(arguments,assembly_dirs,sema):
    sema.acquire()
    report_path,contig_path,genome = arguments
    name = os.path.basename(assembly_dirs)
    # IN REPORT FILE: 1. how many alignments does the BAM file contain. 2. sequence of virus
    report_file = name + '.txt'
    report_file = open(os.path.join(report_path,report_file),"w+")
    report_file = report_path + '/' + name + '.txt'
    subprocess.call(MAP_AND_REPORT % dict(name=name,genome=genome,report_file=report_file,contig_path=contig_path,assembly_dirs=assembly_dirs), shell=True)
    sema.release()

def map_and_report(out,viral_ref,MP_obj):
    assembly_dir = os.path.join(out,'assembly')
    path_list = [os.path.join(assembly_dir,dir) for dir in os.listdir(assembly_dir)]
    print(path_list)
    contig_path = os.path.join(out,'BAM/virome_contigs')
    os.makedirs(contig_path, exist_ok=True)
    report_path = os.path.join(out,'reports')
    os.makedirs(report_path, exist_ok=True)
    MP_obj.mltiprocess_command(map_and_report_command,path_list,[report_path,contig_path,viral_ref])

#%% SORT BAM from contigs

SORT = 'samtools sort %(bam)s -o %(output_bam)s'

def sort_bam_command(noparam,bam,sema):
    sema.acquire()
    output_bam = bam.replace('.bam','.sorted.bam')  
    print(output_bam)
    subprocess.call(SORT % dict(bam=bam,output_bam=output_bam), shell=True)
    sema.release()

def sort_bam(out,MP_obj):
    BAM_DIR = os.path.join(out,'BAM/virome_contigs')
    path_list = [os.path.join(r, f) for r, d, fs in os.walk(BAM_DIR) for f in fs if f.endswith('.mapped.bam')]
    MP_obj.mltiprocess_command(sort_bam_command,path_list,'')
    

#%% DEPTH

DEPTH = 'cut -f2 %(reports)s > %(refs_dir)s/%(sample)s.txt; \
    samtools faidx %(genome)s -r %(refs_dir)s/%(sample)s.txt -o %(refs_dir)s/%(sample)s.fasta; \
        bwa index %(refs_dir)s/%(sample)s.fasta; \
            bwa mem -a -t 10 -v 1 -a %(refs_dir)s/%(sample)s.fasta %(fastq)s %(fastq2)s \
                | samtools view -b - | samtools sort - | tee %(bam_dir)s/%(sample)s.mapped.sorted.bam \
                    | samtools depth -a - > %(depth_dir)s/%(sample)s.depth;'# \
                        #python /mnt/data3/code/bard/automation/depth_allviruses.py %(depth_dir)s/%(sample)s.depth %(depth_dir)s/%(sample)s.pdf'
"""
depth graph of all viruses refseq (all viruses identified from contigs)
vs. original raw sample
"""
def depth_allviruses(depth_dir,sample,genome):
    genome = os.path.dirname(genome)
    genome = os.path.join(genome,'sequences.csv')
    annotations = pd.read_csv(genome)
    depth_file = depth_dir + '/' + sample + '.depth'
    pdf_file = depth_file.replace('.depth', '.pdf')
    file = pd.read_csv(depth_file, header=None, delimiter='\t', names=['acc', 'pos', 'depth'])
    file = file.groupby('acc')
    viruses = [file.get_group(x) for x in file.groups]
    all_accs = sorted(viruses, key=lambda x: entropy(x['depth']), reverse=True)
    with PdfPages(pdf_file) as pdf:
        for df in all_accs:
            depth_col = df['depth']
            enrpy = entropy(depth_col)

            size = len(depth_col)
            coverage = float(sum([x for x in depth_col if x > 0]))/size
            # add annotations! Host, species!

            species = annotations.loc[annotations['Accession'] == df['acc'].iloc[0].split('.')[0], 'Species'].item()
            host = annotations.loc[annotations['Accession'] == df['acc'].iloc[0].split('.')[0], 'Host'].item()

            plt.plot(range(size), depth_col, 'ok', markersize=0.5, rasterized=True)
            plt.xlabel('position')
            plt.ylabel('depth')
            plt.title(f"Species:{species}   Host:{host}")
            plt.suptitle(f"Entropy = {enrpy}    Acc:{df['acc'].iloc[0]}")
            pdf.savefig()
            plt.close()
    
    
def depth_from_source_command(arguments,reports,sema):
    sema.acquire()
    fastq_dir,refs_dir,depth_dir,bam_dir,genome = arguments
    sample = os.path.basename(reports)
    sample = sample.replace('.txt','')  
    fastq = sample + '_R1_001.fastq.gz'
    fastq = os.path.join(fastq_dir,fastq)
    fastq2 = fastq.replace('_R1','_R2')
    #genome = references_to_clean.get('viral_ref')
    print(sample)
    subprocess.call(DEPTH % dict(sample=sample,genome=genome,refs_dir=refs_dir,depth_dir=depth_dir,fastq=fastq,fastq2=fastq2,reports=reports,bam_dir=bam_dir), shell=True)
    depth_allviruses(depth_dir,sample,genome)
    sema.release()

def depth_from_source(args,MP_obj):
    refs_path = os.path.join(args.output,'refs/scaffoldsVirome')
    os.makedirs(refs_path, exist_ok=True)
    depth_path = os.path.join(args.output,'depth')
    os.makedirs(depth_path, exist_ok=True)
    bam_path = os.path.join(args.output,'BAM/virome_fastq')
    os.makedirs(bam_path, exist_ok=True)
    
    report_dir = os.path.join(args.output,'reports')
    path_list = [os.path.join(r, f) for r, d, fs in os.walk(report_dir) for f in fs if f.endswith('.txt')]
    MP_obj.mltiprocess_command(depth_from_source_command,path_list,[args.input,refs_path,depth_path,bam_path,args.viral_ref])
    

#%% INDEX BAM OF FASTQ (and contig)

INDEX = '''samtools index %(bam)s'''
def index_bam_command(no,bam_fastq,sema):
    sema.acquire()  
    print(bam_fastq)
    bam_contig = bam_fastq.replace('virome_fastq','virome_contigs')
    subprocess.call(INDEX % dict(bam=bam_fastq), shell=True)
    subprocess.call(INDEX % dict(bam=bam_contig), shell=True)
    sema.release()

def index_bam(out,MP_obj):
    
    BAM_DIR = os.path.join(out,'BAM/virome_fastq')
    path_list = [os.path.join(r, f) for r, d, fs in os.walk(BAM_DIR) for f in fs if f.endswith('.mapped.sorted.bam')]
    MP_obj.mltiprocess_command(index_bam_command,path_list,'')
    
#%% REPORT FROM FASTQ BAM 

COVERAGE_STATS = '''samtools coverage -H -r %(virus)s %(out)s/BAM/virome_fastq/%(sample)s.mapped.sorted.bam;'''
CONTIG_NUM = '''samtools view -q 1 -c -F 260 %(out)s/BAM/virome_contigs/%(sample)s.mapped.sorted.bam %(virus)s'''
LARGEST_CONTIG = '''samtools stats %(out)s/BAM/virome_contigs/%(sample)s.mapped.sorted.bam %(virus)s | grep ^RL | cut -f 2-'''

def depthEntropy(out_file):
    with open(out_file, 'r') as depth_file:
        depth_dist = [int(x.split()[2]) for x in depth_file]
        entrpy = entropy(depth_dist)
        return entrpy

def report_command(arguments,depth,sema):
    sema.acquire()
    out,breadth_dir,genome = arguments
    sample = os.path.basename(depth)
    sample = sample.replace('.depth','')   
    
    df = pd.DataFrame(columns = ["virusAcc","startpos","length","numreadsmapped","covbases","coverage","meandepth","meanbaseq","meanmapq","numcontigsmapped","longestcontig","entropy"])

    depth_table = pd.read_table(depth,header=None)

    viruses = depth_table.iloc[:,0]
    viruses = sorted(set(viruses))
    for virus in viruses:
        out_file = os.path.join(breadth_dir,sample+'tmp.depth')
        depth_table[depth_table[0] == virus].to_csv(out_file,sep='\t', index=False, header=False)
        entropy = depthEntropy(out_file)
        coverage_stats = subprocess.check_output(COVERAGE_STATS % dict(virus=virus,out=out,sample=sample), shell=True).strip().decode("utf-8")
        contig_num = subprocess.check_output(CONTIG_NUM % dict(virus=virus,out=out,sample=sample), shell=True).strip().decode("utf-8")
        largest_contig = subprocess.check_output(LARGEST_CONTIG % dict(virus=virus,out=out,sample=sample), shell=True).strip().decode("utf-8")
        largest_contig = largest_contig.split('\n')[-1].split('\t')[0]
        new_row = coverage_stats + '\t' + str(contig_num) + '\t' + str(largest_contig) + '\t' + str(entropy)
        new_row = new_row.split('\t')
        df.loc[len(df)] = new_row
    
    df['Accession'] = df['virusAcc'].str.split('.').str.get(0)
    genome = os.path.dirname(genome)
    genome = os.path.join(genome,'sequences.csv')
    annotations = pd.read_csv(genome)
    annotations = annotations[['Accession','Species','Genus','Family','Host']]
    df = pd.merge(df, annotations, how="left", on='Accession')
    df = df.drop(['startpos','Accession'],axis=1)
    os.remove(os.path.join(breadth_dir,sample+'tmp.depth'))
    out_file = os.path.join(breadth_dir,sample+'.csv')
    df.to_csv(out_file,index=False)
    sema.release()

def report(out,viral_ref,MP_obj):
    breadth_path = os.path.join(out,'breadths')
    os.makedirs(breadth_path, exist_ok=True)
    
    depth_DIR = os.path.join(out,'depth')
    path_list = [os.path.join(r, f) for r, d, fs in os.walk(depth_DIR) for f in fs if f.endswith('.depth')]
    MP_obj.mltiprocess_command(report_command,path_list,[out,breadth_path,viral_ref])
    #report_command([out,breadth_path,viral_ref],path_list[0],1)
    

#%% COUNT READS IN FILTERED FASTQ


def count_reads_command(fastq,out,references_to_clean):
    
    #fastq,out = arguments
    print(out)
    fastq2 = fastq.replace('_R1','_R2')
    name = os.path.basename(fastq).replace('_R1_001.fastq.gz','')
    col_df = ['sample','totalReads']
    # TOTAL NUM OF READS
    tot_reads1 = sum(1 for line in gzip.open(fastq)) / 4
    tot_reads2 = sum(1 for line in gzip.open(fastq2)) / 4
    former_reads = int(tot_reads1)+int(tot_reads2)
    tot_reads = [name,former_reads]
    for key in list(references_to_clean.keys()):
        dir_name = 'no' + key.replace('_ref','')
        fastq_path = os.path.join(out,'fastq',dir_name)
        fastq = [os.path.join(r, f) for r, d, fs in os.walk(fastq_path) for f in fs if f.startswith(name)]
        no_reads = sum(1 for line in open(fastq[0])) / 4
        tot_reads.append(former_reads - no_reads)
        former_reads = no_reads
        col_df.append(key.replace('_ref','')+'Reads')
       
    df = pd.DataFrame(columns = col_df)

    
    df.loc[len(df)] = tot_reads
   
    return(df)


def count_reads(args):
    path_list = [os.path.join(r, f) for r, d, fs in os.walk(args.input) for f in fs if f.endswith('_R1_001.fastq.gz') and 'singelton' not in f]
    
    references_to_clean = json.loads(args.ref_clean)
    with mp.Pool(args.threads) as pool:
        dfs = pool.map(partial(count_reads_command, out=args.output, references_to_clean=references_to_clean), path_list)
    #dfs = map(partial(count_reads_command, out=args.output), path_list)
    
    merged_samples = pd.concat(dfs, axis=0) \
        .reset_index()
    merged_samples = merged_samples.drop(columns='index')
    out_path = args.output / ('read_count.csv')
    merged_samples.to_csv(out_path,index=False)
        
#%% Main
if __name__ == "__main__":
    DEBUG = False
    parser = define_parser()
    if DEBUG:
        inline = ['-i','/home/orz/virom_new_ref/raw_fastq/', '-o', '/home/orz/virom_new_ref/', '-c', '{"HG":"/mnt/data3/code/bard/references/Human_Hg38/hg38.fa","Bacteria":"/mnt/data3/code/bard/references/bacteria/bacteria_fullGenomes.fna.gz"}', '-r',"/mnt/data3/code/bard/references/virome/refseq_new_bwa/viral.genomic.all.fasta" ]
        args = parser.parse_args(inline)
        
    else:
        args = parser.parse_args()

    MP_obj = MP(args.threads)
    filter_fastq(args,MP_obj)
    scaffold(args.ref_clean,args.output,MP_obj)
    map_and_report(args.output,args.viral_ref,MP_obj)
    sort_bam(args.output,MP_obj)
    depth_from_source(args,MP_obj)
    index_bam(args.output,MP_obj)
    report(args.output,args.viral_ref,MP_obj)
    count_reads(args)

    
        
    