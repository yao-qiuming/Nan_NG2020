import pandas as pd
import multiprocessing
import gzip
from collections import Counter
from functools import partial
import logging
import sys
import argparse
import os

FORMAT = '%(levelname)s @ %(asctime)s: %(message)s'
DATEFMT = '%m/%d/%Y %I:%M:%S %p'

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG,
                    format=FORMAT,
                    datefmt=DATEFMT)

parser = argparse.ArgumentParser()
parser.add_argument('library', type=str,
                    help="single column tab delimited file (w/ header) " +\
                    "containing sgRNAs")
parser.add_argument('fastq', nargs='+', type=str,
                    help="Space seperated list of fastq files.\n" +\
                    "For umi samples, seperate sgRNA and UMI files by " +\
                    "a \',\' (i.e. /path/to/sgRNAs.fastq,/path/to/UMIs.fastq)")
parser.add_argument('-l', '--labels', nargs='+', type=str,
                    help="Space seperated list of labels in same order as " +\
                    "fastq files.", default="")
parser.add_argument('-n','--name', type=str, help="prefix for output files")
parser.add_argument('-o', '--outdir', type=str,
                    help="/path/to/output.csv")
parser.add_argument('-p','--processes', type=int, default=1,
                    help="Number of processes")

args = parser.parse_args()

fastqs = args.fastq
labels = args.labels
df_lib = pd.read_table(args.library)
lib = set(df_lib.iloc[:,0].values)

if "," in fastqs[0]:
    fastqs = [tuple(fastq.split(",") + [label]) for fastq, label in 
              zip(fastqs, labels)]
else:
    fastqs = zip(fastqs, labels)
    
def gini(counts):
    """
    Adapted from MAGeCK from the Shirley Liu Lab
    """
    counts = pd.np.array(counts)
    xs = sorted(pd.np.log(counts + 1))
    n = len(xs)
    gssum = sum([(i + 1.0) * xs[i] for i in range(n)])
    ysum = sum(xs)
    if not ysum:
        ysum = 1.0
    return 1 - 2.0*(n - gssum/ysum)/(n-1)

def countFastq(lib, files):
    sgrna_handle, umi_handle, sample = files
    with gzip.open(sgrna_handle, 'rt') as sgrnas,\
    gzip.open(umi_handle, 'rt') as umis:
        logging.info("Processing {}".format(sample))
        count = 0
        counter = Counter()
        observed = set()
        for sgrna, umi in zip(sgrnas, umis):
            if count % 4 == 1:
                #sgrna = sgrna.strip()[:20]
                sgrna = sgrna.strip()[1:21]
                if sgrna in lib:
                    umi = umi.strip()[:8]
                    counter[(sgrna, umi)] += 1
                    observed.add(sgrna)
            count += 1
            if count % (10 ** 6) == 0:
                logging.info("Processed {}M Lines of {}".format(count / (10**6),
                      sample))
        zero_counts = len(lib.difference(observed))
        total_reads = int(count/4)
        sorted_keys = sorted(counter)
        df_counter = pd.DataFrame([counter[key] for key in sorted_keys],
                      index=pd.MultiIndex.from_tuples(
                              [key for key in sorted_keys],
                              names=['sgRNA','umi']))
        df_trad = df_counter.groupby('sgRNA').sum()
        gini_idx = gini(df_trad[0].tolist() + [0] * zero_counts)
        mapped_reads = int(df_counter[0].sum())
        percentiles = [0.1, 0.5, 0.9]
        sgrna_read_stats = df_trad[0].describe(percentiles=percentiles).tolist()
        umi_read_stats = df_counter[0].describe(percentiles=percentiles).tolist()
        df_umi = df_counter[~df_counter.index.get_level_values(1).str.contains('N')].groupby(
                'sgRNA').count()
        umi_stats = df_umi[0].describe(percentiles=percentiles).tolist()
        
        stats = [sample, total_reads, mapped_reads, mapped_reads/total_reads,
                 gini_idx, zero_counts] + sgrna_read_stats[1:] + umi_read_stats[1:] +\
                 umi_stats[1:]
        df_counter.columns = [sample]
        df_trad.columns = [sample]
        return df_trad, stats
    
countFastq_mp = partial(countFastq, lib)
p = multiprocessing.Pool(processes=args.processes)
counters = p.map(countFastq_mp, fastqs)
p.close()
p.join()
dfs = [df for df,_ in counters]

vals = ["mean","std","min","10%","50%","90%","max"]
stat_cols = ["sample","total_reads","mapped_reads","percent_reads_mapped",
             "gini_index","zero_counts"] + ["{}_reads/sgrna".format(x) for x in vals] +\
             ["{}_reads/sgrna_umi".format(x) for x in vals] +\
             ["{}_umis/sgrna".format(x) for x in vals]
df_stat = pd.DataFrame([stat for _,stat in counters],
                       columns=stat_cols)
df_stat.to_csv(os.path.join(args.outdir,"{}_stats.csv".format(args.name)),
               index=False)
df_counts = pd.concat(dfs, axis=1)
df_counts.fillna(0, inplace=True)
df_counts.to_csv(os.path.join(args.outdir,"{}_count.csv".format(args.name)),
                 index_label="sgRNA")
