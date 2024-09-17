import pysam
import argparse
import sys
import pandas as pd
import numpy as np
from dataclasses import dataclass


def parse_args(args):
    # Instantiate the parser
    parser = argparse.ArgumentParser(description='Summarize BAM files')
    # Required positional argument
    parser.add_argument('chromosome_position_file', type=str, help='List of chromosome and position')
    parser.add_argument('bamlist', type=str, help='List of bam files to summarize')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Output file')
    # Parse the arguments
    return parser.parse_args(args)


@dataclass
class Alleles:
    allele_list: list

    def __post_init__(self):
        self.allele_tally()
        
    def allele_tally(self):
        allele_array = np.array([a.upper() for a in self.allele_list])
        allele_counts = np.unique(allele_array, return_counts=True)
        self.allele_counts = allele_counts
        
    def allele_balance(self):
        balance = np.min(self.allele_counts[1])/np.sum(self.allele_counts[1])
        return balance

    def allele_count(self):
        return len(self.allele_counts[1])
    
    def depth(self):
        return np.sum(self.allele_counts[1])


def read_chromosome_position_file(chromosome_position_file):
    with open(chromosome_position_file, "r") as f:
        position_list = []
        for line in f:
            position_list.append(line.strip())
    return position_list
            
def read_bam_list(bamlist):
    with open(bamlist, "r") as f:
        bam_list = []
        for line in f:
            bam_list.append(line.strip())
    return bam_list

def summarize_bam_sites(infile, chromosome_position_list):
    samfile = pysam.AlignmentFile(infile, "rb" )
    data = []
    for chrom_pos in chromosome_position_list:
        chromosome, position, *extra = chrom_pos.split()
        position = int(position)
        id = infile.split("/")[-1].replace(".bam", "").replace(".cram", "").replace(".sam", "").replace(".sorted", "")
        data_row = make_bam_data(id, samfile, chromosome, position)
        data.append(data_row)
    samfile.close()
    return pd.concat(data)

def make_bam_data(id, samfile, chromosome, position):
    iter = samfile.pileup(chromosome, position, position+1, truncate=True)
    base_data = {
        "position": [],
        "base": [], 
        "base_qual": [], 
        "mapping_qual": [], 
        "duplicate": [], 
        "paired": [], 
        "proper": [], 
        "rc": [], 
        "secondary": [], 
        "supplementary": [], 
        "unmapped": [], 
        "qc_fail": [], 
        "read1": [], 
        "read2": []}
    for pileupcolumn in iter:
        alleles_list = []
        position = pileupcolumn.pos # 0-based to 1-based
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                alleles_list.append(base)
                base_data["position"].append(position)
                base_data["base"].append(base)
                base_data["base_qual"].append(pileupread.alignment.query_qualities[pileupread.query_position])
                base_data["mapping_qual"].append(pileupread.alignment.mapping_quality)
                base_data["duplicate"].append(pileupread.alignment.is_duplicate)
                base_data["paired"].append(pileupread.alignment.is_paired)
                base_data["proper"].append(pileupread.alignment.is_proper_pair)
                base_data["rc"].append(pileupread.alignment.is_reverse)
                base_data["secondary"].append(pileupread.alignment.is_secondary)
                base_data["supplementary"].append(pileupread.alignment.is_supplementary)
                base_data["unmapped"].append(pileupread.alignment.is_unmapped)
                base_data["qc_fail"].append(pileupread.alignment.is_qcfail)
                base_data["read1"].append(pileupread.alignment.is_read1)
                base_data["read2"].append(pileupread.alignment.is_read2)
        if len(alleles_list) == 0:
            continue
        allele = Alleles(alleles_list)
        balance = allele.allele_balance()
        count = allele.allele_count()
        depth = allele.depth()
        locus_data = pd.DataFrame({"ID": [id], "chromosome":[chromosome], "position": [position], "balance": [balance], "count": [count], "depth": [depth]})
        if len(locus_data) == 0 or len(base_data["position"]) == 0:
            full_data = None       
        else:
            full_data = locus_data.merge(pd.DataFrame(base_data), on="position")
        return full_data

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    chromosome_position_list = read_chromosome_position_file(args.chromosome_position_file)
    bamlist = read_bam_list(args.bamlist)
    outfile = args.outfile
    new_file = True
    include_header = True
    mode = 'w'
    for bam in bamlist:
        data = summarize_bam_sites(bam, chromosome_position_list)
        if new_file:
            data.to_csv(outfile, header = include_header, index = False, sep = "\t", mode = mode)   
            new_file = False
            include_header = False
            mode = 'a'
        else:
            data.to_csv(outfile, header = include_header, index = False, sep = "\t", mode = mode)   

#Unit Tests
class TestClass:
    def test_allele_class(self):
        aa = Alleles(["A", "a", "t", "g", "c"])
        assert aa.allele_balance() == 0.2
        assert aa.allele_count() == 4
        aa.depth() == 5