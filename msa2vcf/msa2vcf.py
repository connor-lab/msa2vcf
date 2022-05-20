#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys
from .readfq import readfq
from operator import itemgetter
from itertools import groupby


def iupac_to_base(base):
    
    lookup_iupac = { 'R' : ['A', 'G'],
                     'Y' : ['C', 'T'],
                     'S' : ['G', 'C'],
                     'W' : ['A', 'T'],
                     'K' : ['G', 'T'],
                     'M' : ['A', 'C'],
                     'B' : ['C', 'G', 'T'],
                     'D' : ['A', 'G', 'T'],
                     'H' : ['A', 'C', 'T'],
                     'V' : ['A', 'C', 'G'] }

    if lookup_iupac.get(base):
        return lookup_iupac.get(base)

    else:
        return base



def get_ref_seq(msa, refname):

    with open(msa) as f:
        for name, seq, qual in readfq(f):
            if name == refname:
                return seq.upper()
                break

def group_dels(dels):

    grouped_dels = []

    for k,g in groupby(enumerate(dels),lambda x:x[0]-x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        grouped_dels.append((group[0] - 1 ,len(group) + 1))

    return grouped_dels

def group_ins(ins):
    res =  [(el - 1, ins.count(el) + 1) for el in ins] 
    
    return set(res)

def fix_complex_vars(all_vars, rseq, qseq):

    positions = [i[2] for i in all_vars]

    duplicate_positions = set([x for x in positions if positions.count(x) > 1])

    for pos in duplicate_positions:
        duplicate_variants = [x for x in all_vars if x[2] == pos]
        for variant in duplicate_variants:
            if variant[0] == "snp":
                all_vars.remove(variant)

            elif variant[0] == "ins":
                variant[5] += 1
                variant[2] -= 1
                variant[4] -= 1
                variant[1] = rseq[variant[4]]
                variant[3] = qseq[variant[4]:variant[4]+variant[5]]

            elif variant[0] == "del":
                variant[5] += 1
                variant[2] -= 1
                variant[4] -= 1
                variant[3] = rseq[variant[4]]
                variant[1] = rseq[variant[4]:variant[4]+variant[5]]

    return all_vars

def update_snps(qseq, rseq):

    dels = []
    ins = []
    snps = []

    for qpos, qbase in enumerate(qseq):
        ins_aware_pos = qpos - len(ins)

        if qbase != rseq[qpos]:
            if qbase == "-":
                # deletion
                dels.append(ins_aware_pos)

            elif rseq[qpos] == "-":
                # insertion
                ins.append(ins_aware_pos)
                
            else:
                # snp
                snps.append((ins_aware_pos, qbase))

        elif qbase == "-" and rseq[qpos] == "-":
            # insertion but not in this sequence
            ins.append(ins_aware_pos)

    all_vars = []

    if dels:
        grouped_dels = group_dels(dels)
        for start,length in group_dels(dels):
            ins_correction = len([i for i in ins if i < start])
            var = ["del", rseq[start+ins_correction:start+ins_correction+length], start, rseq[start+ins_correction], start+ins_correction, length ]

            all_vars.append(var)

    if snps:
        for position,base in snps:
            ins_correction = len([i for i in ins if i < position])
            var = ["snp", rseq[position+ins_correction], position, base, position+ins_correction, 1]
            all_vars.append(var)

    if ins:
        for start,length in group_ins(ins):
            ins_correction = len([i for i in ins if i < start])
            if rseq[start+ins_correction:start+ins_correction+length] ==  qseq[start:start+length]:
                # insertion but not in this sequence
                continue
            var = ["ins", rseq[start+ins_correction], start, qseq[start:start+length], start+ins_correction, length]
            all_vars.append(var)

    if len(all_vars) >= 1:
         fixed_vars = fix_complex_vars(all_vars, rseq, qseq)

         for var in fixed_vars:
             deconvolute_IUPAC(var)
             #print(var)

         fixed_vars_s = sorted(fixed_vars, key=lambda x: x[2])

         return fixed_vars_s

    else:
        return None

def deconvolute_IUPAC(var):

    num_alts = 1

    for base in var[3]:
        no_iupac = iupac_to_base(base)
        if isinstance(no_iupac, list):
            num_alts = len(no_iupac)
            no_iupac.remove(var[1])
            var[3] = ','.join(no_iupac)

    var.append(round(1 / num_alts, 2))

    return var
    

def make_vcf(snps, qname, rname, keep_n):

    vcflines = []    

    # header
    vcflines.append("##fileformat=VCFv4.2")
    vcflines.append("##source="+os.path.basename(sys.argv[0]))
    vcflines.append("##contig=<ID=" + rname + ">")
    vcflines.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">")
    vcflines.append("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">")
    vcflines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+qname)

    # variants
    for line in snps:
         vcf_line = None 

         alt_no_n = re.sub(r'(N|-)+', '', line[3])

         alt_no_gap = re.sub(r'-+', '', line[3])

         if alt_no_n:
             if not line[1] == alt_no_n:
                 vcf_line = '\t'.join([rname, str(line[2] + 1 ), ".", line[1], alt_no_gap, ".", "PASS", "DP=1;AF="+str(line[6])])

         else:
             if keep_n:
                 vcf_line = '\t'.join([rname, str(line[2] + 1 ), ".", line[1], alt_no_gap, ".", "PASS", "DP=1;AF="+str(line[6])])

         if vcf_line:
             vcflines.append(vcf_line)

        
    return vcflines


def write_vcf(vcflines, qname):
    with open(qname +".vcf", 'w') as f:
        for line in vcflines:
            f.write(line+"\n")

def remove_terminal_gapns(seq):
    return re.sub(r'(N|-)*$', '', seq)

def go(args):

    refseq = get_ref_seq(args.msa, args.refname)

    if not refseq:
        print(args.refname + " not found\n") 
        sys.exit(1)

    with open(args.msa) as f:
        for name, seq, qual in readfq(f):
            if name == args.refname:
                continue
            else:
                print(name)
                snps = update_snps(remove_terminal_gapns(seq.upper()), refseq)

                if snps:
                    vcflines = make_vcf(snps, name, args.refname, args.keep_n)
                    write_vcf(vcflines, name)

