#!/usr/bin/env python

import regex
import gzip
import argparse
import os
import sys
import warnings
import math
import multiprocessing as mp
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord,_RestrictedDict
from sklearn.neighbors import KDTree
from itertools import product,combinations
#from scipy.spatial.distance import hamming
import pybktree

if not sys.warnoptions:
    warnings.simplefilter("ignore")

def warn(message):
    warnings.warn(message)
    print(message)

parser = argparse.ArgumentParser(description="""
Author: Jianhong Ou @ duke, July, 2023
This source code is licensed under the MIT license
This code depends on regex, gzip, argparse and Biopython.
This code will extract read by R2 reads if R2 reads fit the regrex string
eg:
    pattern: "(.{17,20})GGATTCGAGGAGCGTGTGCGAACTCAGACCA{i<=2,d<=2,e<=3}(.{6})ATCCACGTGCTTGAGAGGCCAGAGCATTCG(((?P<TYPE>AG).{3})|((?P<TYPE>TC).{3}))"
    (.{17,20}) will match the first 17~20 nucleotide until reached GGATT...
    GGATT...GCGTGTGC...ACCA{i<=2,d<=2,e<=3} will match the sequence with max 2 insertion, 2 deletion, and total 3 mismatch
        The GGATT...GCGT is from the 3' end of third round primer
        The GTGC...ACCA is from the 5' end of the second round primer
    (.{6}) will match the third round barcode (6 nucleotide)
    ATC...GAGAGG...CG will match the sequence, you can also allow max 2 mismatch, 1 insertion, 1 deletion by change it to
        ATCCACGTGCTTGAGAGGCCAGAGCATTCG{i<=1,d<=1,e<=2}
    (((?P<TYPE>AG).{3})|((?P<TYPE>TC).{3})) will match the first round insertion barcode. Please note that (?P<TYPE>AG) and (?P<TYPE>TC) must named with TYPE
The code will rewrite R1 and R2 reads:
    the id will with the matched barcode
    the R2 reads will be the concatenated barcode
If you want to run cellRanger, please swith the output reads name (R2->R1) and make sure the pattern follow 16(BAR)+10(UMI) rule
 eg:
     pattern: pattern ="(?P<UMI>.{10})(?P<BAR>.{7})(.{0,3})GGATTCGAGGAGCGTGTGCGAACTCAGACCA{i<=2,d<=2,e<=3}(?P<BAR>.{6})ATCCACGTGCTTGAGAGGCCAGAGCATTCG(((?P<TYPE>AG)(?P<BAR>.{3}))|((?P<TYPE>TC)(?P<BAR>.{3})))"

example:
  python trim_R2.py -i1 input_R1_001.fastq.gz -i2 input_R2_001.fastq.gz \\
        -d1 dna_R1_001.fastq.gz -d2 dna_R2_001.fastq.gz \\
        -r1 rna_R1_001.fastq.gz -r2 rna_R2_001.fastq.gz \\
        -n 500 -rb TC -db AG \\
        -x -w 3M-february-2018.txt.gz \\
        -p "(?P<UMI>.{10})(?P<BAR>.{7})(.{0,3})GGATTCGAGGAGCGTGTGCGAACTCAGACCA{i<=2,d<=2,e<=3}(?P<BAR>.{6})ATCCACGTGCTTGAGAGGCCAGAGCATTCG(((?P<TYPE>AG)(?P<BAR>.{3}))|((?P<TYPE>TC)(?P<BAR>.{3})))"
""", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i1', '--read1', type=str, required=True, help="input R1 reads")
parser.add_argument('-i2', '--read2', type=str, required=True, help="input R2 reads")
parser.add_argument('-r1', '--outRNA1', type=str, required=True, help="output RNA R1 reads")
parser.add_argument('-r2', '--outRNA2', type=str, required=True, help="output RNA R2 reads")
parser.add_argument('-d1', '--outDNA1', type=str, required=True, help="output DNA R1 reads")
parser.add_argument('-d2', '--outDNA2', type=str, required=True, help="output DNA R2 reads")
parser.add_argument('-rb', '--rnaBarcode', type=str, required=True, help="rna barcode")
parser.add_argument('-db', '--dnaBarcode', type=str, required=True, help="dna barcode")
parser.add_argument('-p', '--pattern', type=str, required=True, help="pattern for R2 reads to match")
parser.add_argument('-n', '--nrecord', type=int, default=0, help="number of record to be read")
parser.add_argument('-x', '--cellranger', action='store_true', help="output for cellranger")
parser.add_argument('-w', '--barcodeWhitelist', type=str, help="filename of barcode whitelist for cellranger")
parser.add_argument('-m', '--barcodeMap', type=str, default="barcodeMap.tsv.gz", help="output filename for barcode file. The tsv file contains two columns originalBarcode and cellRangerBarcode")
parser.add_argument('-b', '--originBarcodePrefix', type=str, default="ori_", help="prefix of output filename for original barcode R2 file.")
parser.add_argument('-c', '--cores', dest="cores", default=1, type=int, required=False, help="number of cores for multiprocessing")

args = parser.parse_args()

in_r1=args.read1 #in_r1="LY27_R1_001.fastq.gz"
in_r2=args.read2 #in_r2="LY27_R2_001.fastq.gz"
out_rna_r1=args.outRNA1 #out_rna_r1="r_R1_001.fastq.gz"
out_rna_r2=args.outRNA2 #out_rna_r2="r_R2_001.fastq.gz"
out_dna_r1=args.outDNA1 #out_dna_r1="d_R1_001.fastq.gz"
out_dna_r2=args.outDNA2 #out_dna_r2="d_R2_001.fastq.gz"
pattern =args.pattern #pattern ="(?P<UMI>.{10})(?P<BAR>.{7})(.{0,3})GGATTCGAGGAGCGTGTGCGAACTCAGACCA{i<=2,d<=2,e<=3}(?P<BAR>.{6})ATCCACGTGCTTGAGAGGCCAGAGCATTCG(((?P<TYPE>AG)(?P<BAR>.{3}))|((?P<TYPE>TC)(?P<BAR>.{3})))"
N=args.nrecord
cellranger=args.cellranger
rna_barcode=args.rnaBarcode
dna_barcode=args.dnaBarcode
cores=args.cores

if cellranger:
    handle_barcodeWhitelist=gzip.open(args.barcodeWhitelist, "rt")
    barcodeWhitelist={line.rstrip('\n'):True for line in handle_barcodeWhitelist}
    handle_barcodeWhitelist.close()
    barcode_map = {}
    ori_barcode = []
    handle_barcode_map = gzip.open(args.barcodeMap, "wt")
    originBarcodePrefix = args.originBarcodePrefix
    rna_r2_realpath = os.path.abspath(out_rna_r2)
    out_ori_rna_r2 = os.path.join(os.path.dirname(rna_r2_realpath), args.originBarcodePrefix+os.path.basename(rna_r2_realpath))
    handle2_rna_out = gzip.open(out_ori_rna_r2, "wt")
    dna_r2_realpath = os.path.abspath(out_dna_r2)
    out_ori_dna_r2 = os.path.join(os.path.dirname(dna_r2_realpath), args.originBarcodePrefix+os.path.basename(dna_r2_realpath))
    handle2_dna_out = gzip.open(out_ori_dna_r2, "wt")
    handle_new_rna_r2 = gzip.open(out_rna_r2, "wt")
    handle_new_dna_r2 = gzip.open(out_dna_r2, "wt")
else:
    handle2_rna_out = gzip.open(out_rna_r2, "wt")
    handle2_dna_out = gzip.open(out_dna_r2, "wt")

handle1_in = gzip.open(in_r1, "rt")
handle2_in = gzip.open(in_r2, "rt")
handle1_rna_out = gzip.open(out_rna_r1, "wt")
handle1_dna_out = gzip.open(out_dna_r1, "wt")

itter1 = SeqIO.parse(handle1_in, "fastq")
itter2 = SeqIO.parse(handle2_in, "fastq")
nr = 0
for record_r2 in itter2:
    record_r1=next(itter1)
    nr = nr + 1
    if nr>N & N!=0:
      break
    if 0 == nr % 10000000:
        print("handle reads %s\n" % nr)
    m = regex.match(pattern, str(record_r2.seq))
    if m:
        seq0 = str(record_r2.seq)
        quality0 =record_r2.letter_annotations
        qk = list(quality0)[0]
        seq = ""
        nid = []
        quality = []
        if cellranger:
            bar = ""
            q_bar = []
            for i in m.spans("BAR"):
                id = slice(i[0], i[1])
                bar +=seq0[id]
                nid.append(seq0[id])
                q_bar +=quality0[qk][id]
            seq = bar[:16]
            if not seq in ori_barcode:
                ori_barcode.append(seq)
            quality = q_bar[:16]
            umi = ""
            q_bar = []
            for i in m.spans("UMI"):
                id = slice(i[0], i[1])
                umi +=seq0[id]
                nid.append(seq0[id])
                q_bar +=quality0[qk][id]
            seq += umi[:10]
            quality += q_bar[:10]
        else:
            idx = [i[0] for i in m.allspans()[1:] if len(i) ]
            id_arr = []
            for i in idx:
                id = slice(i[0], i[1])
                nid.append(seq0[id])
                id_arr += range(i[0], i[1])
            id_arr = set(sorted(id_arr))
            for i in id_arr:
                seq += seq0[i]
                quality.append(quality0[qk][i])
        id=m.spans("TYPE")[0]
        TYPE = seq0[slice(id[0], id[1])]
        rna=rna_barcode==TYPE
        dna=dna_barcode==TYPE
        nid.insert(0, TYPE)
        nid = ':'.join(nid)
        nq = _RestrictedDict(length=len(seq))
        nq[qk] = quality
        recordN = SeqRecord(Seq(seq), id=record_r2.id, description=nid, letter_annotations=nq)
        record_r1.description=nid
        if rna:
            handle2_rna_out.write(recordN.format("fastq"))
            handle1_rna_out.write(record_r1.format("fastq"))
        else:
            if dna:
                handle2_dna_out.write(recordN.format("fastq"))
                handle1_dna_out.write(record_r1.format("fastq"))

handle1_in.close()
handle2_in.close()
handle1_rna_out.close()
handle2_rna_out.close()
handle1_dna_out.close()
handle2_dna_out.close()

print("reads barcode split done\n")

def checkKey(key, barcodeWhitelist):
    try:
        return barcodeWhitelist[key]
    except:
        return False

def hamming_distance(a, b):
    """
    hamming distance for fixed width string a and b (both are 16 bp)
    """
    if len(a) != len(b):
        return len(a) if len(a) > len(b) else len(b)
    return sum([not x == y for x, y in zip(a, b)])

# Function to generate k-mers from DNA sequence
def generate_kmers(dna_sequence, k):
    return [dna_sequence[i:i+k] for i in range(len(dna_sequence) - k + 1)]

# Function to generate the unique k-mers for a given K
def generate_unique_kmers(k):
    return [''.join(c) for c in product('ACGT', repeat=k)]

# Function to convert DNA sequence into k-mer counts vector
def dna_to_kmer_counts(dna_sequence, k, unique_kmers):
    kmers = generate_kmers(dna_sequence, k)
    kmer_counts = [kmers.count(kmer) for kmer in unique_kmers]
    return np.array(kmer_counts)

# Function get one barcode by nearest neighbor, not promising the best match
def getKNN(query, subject, kdtree, k, unique_kmers):
        # Convert query sequence to k-mer counts vector
        query_vector = dna_to_kmer_counts(query, k, unique_kmers)
        # Find the nearest neighbor and its distance
        distances, indices = kdtree.query([query_vector], k=1)
        # Get the index of the nearest neighbor
        nearest_neighbor_index = indices[0][0]
        # Get the DNA sequence of the nearest neighbor
        nearest_neighbor_dna = subject[nearest_neighbor_index]
        return nearest_neighbor_dna

# Function get first barcode that not used
def getFirstBar(dic):
    for k,v in dic.items():
        if v:
            return k
    return ''

# Function the update the barcode_map, whitelist and unique_barcode
def updateBarcodesAndWhitelist(barcode_map, barcodeWhitelist, unique_barcode, bar, new_bar):
    barcode_map[bar] = new_bar
    barcodeWhitelist[new_bar] = False
    tmp = unique_barcode.remove(bar)

# multiple process workers
def worker1(ori_barcode_srt, barcodeWhitelist, notuse):
    barcode_map = {}
    for bar in ori_barcode_srt:
        if checkKey(bar, barcodeWhitelist):
           barcode_map[bar] = bar
    return barcode_map

def worker2(ori_barcode_srt, barcodeWhitelist, bk_tree):
    barcode_map = {}
    for bar in ori_barcode_srt:
        candidates = bk_tree.find(bar, 1)
        for d,b in candidates:
            barcode_map[bar] = b
            break
    return barcode_map

# callback
def collect_barcode_map(b_map):
    for bar,new_bar in b_map.items():
        barcode_map[bar] = new_bar
        barcodeWhitelist[new_bar] = False
        tmp = unique_barcode.remove(bar)

# multiple process listener, stop listen utill get keywords kill
def listener(q):
    while 1:
        m = q.get()
        if m == "kill":
            break

def find_barcode(cores, workerFUN, ori_barcode_srt, barcode_map, barcodeWhitelist, unique_barcode, additionalParam):
    manager = mp.Manager()
    q = manager.Queue()
    if cores > mp.cpu_count():
        cores = mp.cpu_count()
    pool = mp.Pool(cores)
    watcher = pool.apply_async(listener, q)
    jobs = []
    N = len(ori_barcode_srt)
    step = math.ceil(N/cores)
    for i in range(0, N, step):
        f = 0+i
        t = step+i
        if t > N:
            t = N
        job = pool.apply_async(workerFUN, (ori_barcode_srt[f:t], barcodeWhitelist, additionalParam), callback=collect_barcode_map)
        jobs.append(job)
    for job in jobs:
        job.get()
    q.put("kill")
    pool.close()
    pool.join()

if cellranger:
    print("Remap barcode to whitelist\n")
    barcode_map = {}
    # get unique barcodes
    unique_barcode = list(set(ori_barcode))
    if len(unique_barcode) > len(barcodeWhitelist):
        message = '''\
            The length of barcode ({blen}}) is longer than that of whitelist({wlen}).
            If you do not want to ran cellranger, you can stop at anytime.\
        '''.format(blen=len(unique_barcode), wlen=len(barcodeWhitelist))
        warn(message)

    print("Remapping barcode by hamming distance 1\n")
    ## for distance 0 and 1, one whitelist barcode may map to multiple barcodes
    # distance 0
    ori_barcode_srt = sorted(unique_barcode)
    if __name__ == '__main__':
        fb = find_barcode(cores, worker1, ori_barcode_srt, barcode_map, barcodeWhitelist, unique_barcode, [])

    # distance 1
    ori_barcode_srt = sorted(unique_barcode)
    bk_tree = pybktree.BKTree(hamming_distance, barcodeWhitelist.keys())
    if __name__ == '__main__':
        fb = find_barcode(cores, worker2, ori_barcode_srt, barcode_map, barcodeWhitelist, unique_barcode, bk_tree)
    print("hamming distance 1 done! Continue remapping.\n")

    # use nearest neighbor to speed up, this may be not the closest word because
    # the K=3 is not good for alignment. K=7 may work but ask too much memory.
    if len(unique_barcode):
        k = 2 # 3 mer ask at least 16G memory, k=2 is two-level, the distance is about 8
        unique_kmers = generate_unique_kmers(k)
        available_dna = [dna_sequence for dna_sequence in barcodeWhitelist.keys() if barcodeWhitelist[dna_sequence]]
        data_vectors = [dna_to_kmer_counts(dna_sequence, k, unique_kmers) for dna_sequence in available_dna]
        # Create a KDTree with the data vectors
        kdtree = KDTree(data_vectors)
        ori_barcode_srt = sorted(unique_barcode)
        for bar in ori_barcode_srt:
            new_bar = getKNN(bar, available_dna, kdtree, k, unique_kmers)
            if new_bar:
                if not checkKey(new_bar, barcodeWhitelist):
                    new_bar = getFirstBar(barcodeWhitelist)
                    if not new_bar:
                        # no barcode left in barcodeWhitelist
                        break
            else:
                new_bar = getFirstBar(barcodeWhitelist)
                if not new_bar:
                    break
            barcode_map[bar] = new_bar
            barcodeWhitelist[bar] = False
            unique_barcode.remove(bar)

    if len(unique_barcode):
        warn("There are %0.2f%%(%d/%d) barcode can not map to given whitelist!" % (100*len(unique_barcode)/len(set(ori_barcode)), len(unique_barcode), len(set(ori_barcode))))
        for bar in unique_barcode:
            bacode_map[bar] = bar

    print("Remap barcode done\nWriting new barcode files\n")
    handle_barcode_map.write("originalBarcode\tcellRangerBarcode\n")
    for key,value in barcode_map.items():
        handle_barcode_map.write("{}\t{}\n".format(key, value))
    handle_barcode_map.close()

    handle2_rna_out=gzip.open(out_ori_rna_r2, "rt")
    itter2 = SeqIO.parse(handle2_rna_out, "fastq")
    for record_r2 in itter2:
        seq0 = str(record_r2.seq) # seq0 = barcode[0:16] + UNI[16:]
        record_r2.seq = Seq(barcode_map[seq0[:16]] + seq0[16:])
        handle_new_rna_r2.write(record_r2.format("fastq"))
    handle_new_rna_r2.close()
    handle2_rna_out.close()

    handle2_dna_out=gzip.open(out_ori_dna_r2, "rt")
    itter2 = SeqIO.parse(handle2_dna_out, "fastq")
    for record_r2 in itter2:
        seq0 = str(record_r2.seq)
        record_r2.seq = Seq(barcode_map[seq0[:16]] + seq0[16:])
        handle_new_dna_r2.write(record_r2.format("fastq"))
    handle_new_dna_r2.close()
    handle2_dna_out.close()
