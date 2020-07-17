#!/usr/bin/env python
# coding: utf-8


#STEP 1 - MAKES INDIVIDUAL GENOTYPE FILES

###############
#   Imports  #
###############
import re
from collections import defaultdict
import argparse


################
#   Functions   #
################

def get_arguments():
    parser = argparse.ArgumentParser(description="This is a program that takes three arguments: a VCF file, a sample name, and a path for outfiles.  Returns a genotype file with 4 columns: sample name, marker, refernce genotype, and samples genotype")
    parser.add_argument("-vcf","--vcf_file", help="-vcf <path><file>, Requires a variant call format file, assumes sequences are in VCFv4 format. e.g. path/to/file/sample.vcf",required=True,type=str)
    parser.add_argument("-s","--sample_name", help="-s <path><file>, Requires a string, name of current sample. e.g. R20 ",required=True,type=str)
    parser.add_argument("-p","--out_path", help="-p <path><file>, Requires a path for output. e.g. path/to/file/",required=True,type=str)
    return parser.parse_args()


def create_out_file_name(vcf_file):
    file_parts = vcf_file.split("/")
    file_name = file_parts[-1]
    new_file_name_parts = file_name.split("vcf")
    new_file_name = new_file_name_parts[0]
    out_file_name = out_path + new_file_name + "fixed_genotype.txt"
    return out_file_name
        

####################
# Global Variables #
####################

args = get_arguments()
vcf_file = args.vcf_file #vcf file
sample_name = args.sample_name #sample name
out_path = args.out_path #path to vcf_file

ref_dict = {}
alt_dict = {}
called_geno_number_dict = {}
ref_geno_number_dict = {}
called_genos = defaultdict(list)
ref_genos = defaultdict(list)

###################
#  Opening Files  #
###################

vcf = open(vcf_file, "r")
out_file = create_out_file_name(vcf_file)
out_fh = open(out_file,"w")

############
#   Main   #
############

############################################################
#Getting lines in VCF file
for line in vcf:
    line = line.strip("\n")
    if line.startswith("#"):
        continue
    else:

############################################################
#Assigning vairables to line parts to create marker names and to get reference and alternative alleles. 

        parts = line.split("\t")
        plain_marker = parts[0]+"_"+parts[1]
        ref_allele = parts[3]
        geno_parts = parts[9]
        geno_parts = geno_parts.split(":")
        geno_parts = geno_parts[0]
        geno = geno_parts.split("/")
        called_geno = [int(i) for i in geno]
        ref_geno = [0 if i < 100 else i for i in called_geno]
        ref_geno = [int(i) for i in ref_geno]
        alt_alleles = parts[3]+","+parts[4]
        alt_alleles = alt_alleles.split(",")
        
############################################################        
# Creating 4 dictionaries that respectively contain: 
#marker and reference allele. e.g. {'Ro06_250016': ['T'], 'Ro06_250039': ['C']}
#marker and alternative alleles. e.g. {'Ro06_250016': ['C', 'T'], 'Ro06_250039': ['C', 'T']}
#marker and genotype from VCF. e.g. {'Ro06_250016': [0, 0, 1, 1, 1, 1], 'Ro06_250039': [0, 0, 0, 0, 0, 1]}
#marker and the reference genotypes from VCF. e.g. {'Ro06_250016': [0, 0, 0, 0, 0, 0], 'Ro06_250039': [0, 0, 0, 0, 0, 0]}

        ref_dict[plain_marker] = [ref_allele]
        alt_dict[plain_marker] = alt_alleles
        called_geno_number_dict[plain_marker] = called_geno
        ref_geno_number_dict[plain_marker] = ref_geno
        
############################################################
# Samples genotype - Substituting the called numbers in the genotype field with the actual reference and alternative allele calls.

for k1, v1 in called_geno_number_dict.items():
    for k2,v2 in alt_dict.items():
        if k1 == k2:
            for i in v1:
                called_genos[k1].append(v2[i])

############################################################
# Reference Genotype - Substituting the called numbers in the genotype field with the actual reference.

for k3, v3 in ref_geno_number_dict.items():
    for k4, v4 in ref_dict.items():
        if k3 == k4:
            for j in v3:
                ref_genos[k3].append(v4[j])

############################################################
# Creating the sample's genotype file. 


out_fh.write("Sample" + "\t" + "Marker" + "\t" + "Reference Genotype" + "\t" + "Sample Genotype" + "\n")

for k5, v5 in called_genos.items():
    for k6, v6 in ref_genos.items():
        if k5 == k6:
            sample_genotype = "/".join(v5)
            reference_genotype = "/".join(v6)
            out_fh.write(sample_name + "\t" + k5 + "\t" + reference_genotype + "\t" + sample_genotype + "\n")
            
out_fh.close()



