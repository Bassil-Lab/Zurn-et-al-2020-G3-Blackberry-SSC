#!/usr/bin/env python
# coding: utf-8

#STEP 2 - MAKES COMBINED GENOTYPE TABLE WITH DEPTHS ACCOUNTED FOR - NEW CODE

###############
#   Imports   #
###############
import re
import collections
from collections import defaultdict
import argparse

###############
#  Functions  #
###############

def get_arguments():
    parser = argparse.ArgumentParser(description="This is a program that takes four arguments: a file that contains the path and file names all genotype calls from step 1 for each sample, a file with the path and file names to all the deth files for each sample, a file with \
    all the name sample names (all names must be unique), a string with the path and file name for output. The input file must be in the same order. Returns a tab seperated genotype table with the marker and each sample's genotype as columns.")
    parser.add_argument("-g","--sample_genotypes", help="-g <path><file>, Requires a string for the text file that contains all paths and file names to the files created in step 1. All files must be in the same order. e.g. path/to/my/genotypes.txt",required=True,type=str)
    parser.add_argument("-d","--sample_depths", help="-d <path><file>, Requires a string for the text file that contains all paths and file names to the depth files created by Samtools. All files must be in the same order. e.g. path/to/my/depths.txt",required=True,type=str)
    parser.add_argument("-s","--sample_names", help="-s <path><file>, Requires a string for the text file that contains all the sample names. All files must be in the same order. e.g. path/to/my/sample/names.txt",required=True,type=str)
    parser.add_argument("-o","--out_file", help="-o <path><file>, Requires a string that contains the path and file name for output. e.g. path/for/out/file.tsv",required=True,type=str)
    return parser.parse_args()


# sample_genotypes="./masterfile.txt"
# sample_depths="./depth_masterfile.txt"
# sample_names="./depth_samples.txt"
# out_file="./testmastergenotable.tsv"

####################
# Global Variables #
####################
args = get_arguments()
sample_genotypes = args.sample_genotypes #sample genotypes file
sample_depths = args.sample_depths #sample depths file
sample_names = args.sample_names #sample names file
out_file = args.out_file #path and file name for out file


samp_list=[]
marker_list=[]
geno_list=[]
ref_list=[]
master_marker_list=[]
new_marker_list=[]
refsamplist=[]

ref_dict={}
depth_dict=defaultdict(list)
pos_dict={}
dict_entries={}
d = {}
missgenodict={}
sampledepthdict={}
refmarkerdict={}

###################
#  Opening Files  #
###################

filelist = open(sample_genotypes, 'r')
depthlist = open(sample_depths, 'r')
namefile=open(sample_names, 'r')
out_f=open(out_file,"w")

############
#   Main   #
############

######################################
#Getting lines in sample files.

for file in filelist:
    file=file.strip("\n")
    f=open(file, 'r')
    for file_line in f:
        if file_line.startswith("S"):
            continue
            
######################################
#Assigning Variables to parts of the samples line. 

        else:
            file_line=file_line.strip("\n")
            parts=file_line.split("\t")
            sample=parts[0]
            marker=parts[1]
            ref=parts[2]
            alt=parts[3]
            
######################################
#Creating lists for each vairable to use later.
                    
            samp_list.append(sample)
            marker_list.append(marker)
            geno_list.append(alt)
            
######################################
#Creating a dictionary that contains the marker and the reference allele.
            ref_dict[marker]=ref
            
######################################
#Creating a list that contains all the samples names.
uniq_samps=list(set(samp_list))

######################################
#Getting the lines in depth files.
for dfile, nline in zip(depthlist, namefile):
    nline=nline.strip()
    dfile=dfile.strip("\n")
    df=open(dfile, 'r')
    
######################################
#Creating variable for parts of the depth line. 

    for dline in df:
        nline=nline.strip()
        dline=dline.strip("\n")
        dline_parts=dline.split("\t")
        depth = dline_parts[2]
        depthmarker= dline_parts[0]+"_"+dline_parts[1]
        
######################################
#Removing all zero depths. Zero depths and non-present depths will be considered missing genotypes.

        if int(depth) != 0:

######################################
#Creating a dictionary with the marker(pos chrom) and all sample names from the depth files.
#e.g. {'Ro06_244699':['R20', '35'], 'Ro06_244700': ['R20', '35']}
        
            depth_dict[depthmarker].append(nline)
            
######################################
#Creating a dictionary for table format
#Dictionary contains all the unique makers from all files and a list of all sample names that will be used to subsitute in genotype inforamtion.
#e.g. {'Ro06_241281': ['35', 'R37', 'R20'], 'Ro06_241285': ['35', 'R37', 'R20']}

for marker in marker_list:
    if marker not in pos_dict:
        pos_dict[marker] = uniq_samps
        
######################################
#Creating a list of dictionaries that contain the sample name and the genotype.
#e.g. [{'35': 'A/A/A/A'}, {'35': 'AAA/AAA/AAA/AAA'}]
        
master_geno_list=([{k: v} for k,v in zip(samp_list, geno_list)])

######################################
#Creating a list of dictionaries that contain all uniq markers.
#e.g. [{'marker': 'Ro06_241281'}, {'marker': 'Ro06_241285'}]

for marker,geno in zip(marker_list,geno_list):

    dict_entries["marker"]=marker
    master_marker_list.append(dict(dict_entries))

######################################
#Creating a list of dictionaries that contains the marker and genotype for each sample.
#e.g. [{'marker': 'Ro06_241281', '35': 'A/A/A/A'}, {'marker': 'Ro06_241285', '35': 'AAA/AAA/AAA/AAA'}]

combine_list = [{**marker, **geno} for marker, geno in zip(master_marker_list, master_geno_list)]

######################################
#Creating a dictionary that contains dictionaries that have the all samples and genotype for each marker.
# e.g. {'Ro06_241281': {'marker': 'Ro06_241281', '35': 'A/A/A/A', 'R37': 'A/A/A/A'}, 'Ro06_241285': {'marker': 'Ro06_241285', '35': 'AAA/AAA/AAA/AAA', 'R37': 'AAA/AAA/AAA/AAA'}}

for s in combine_list:
    d.setdefault(s['marker'], {}).update(s)
    
######################################
#Updating dictionary table to contain genotype infromation for each marker. 
#e.g. {'Ro06_241281': ['A/A/A/A', 'A/A/A/A', 'R20'], 'Ro06_241285': ['AAA/AAA/AAA/AAA', 'AAA/AAA/AAA/AAA', 'R20'], 'Ro06_241302': ['CATATG/CATATG/CATATG/CATATG', 'CATATG/CATATG/CATATG/CATATG', 'R20']

d1 = {k: [d[k].get(i, i) for i in l] for k, l in pos_dict.items()}

######################################    
#Idenitifying missing genotypes from each marker.  
#e.g. R20 is missing for Ro06_241281, and R20 is missing for Ro06_241285

for k1,v1 in d1.items():
    missinggenos=[value for value in v1 if value in uniq_samps]

######################################    
#Creating a dictionary that contains all the missing genotypes for each marker.
#e.g. {'Ro06_241281': ['R20'], 'Ro06_241285': ['R20']}

    missgenodict[k1]=missinggenos

######################################
#Identifying which of the missing samples actully have depths.  
#e.g. {'Ro06_244736': ['R20'], 'Ro06_244742': ['R20']}

for k2,v2 in depth_dict.items():
    for k3,v3 in missgenodict.items():
        if k2 == k3:
            depthpresentsamps=[value for value in v2 if value in v3]
            sampledepthdict[k2]=depthpresentsamps

######################################
#Creating markers associated with the missing genotypes. 
#[{'marker': 'Ro06_244736'}, {'marker': 'Ro06_244742'}]
for k4,v4 in sampledepthdict.items():
    for k5,v5 in ref_dict.items():
        if k4 == k5:
            refsamplist.append("".join(v4))
            ref_list.append(v5)
            refmarkerdict["marker"]=k4
            new_marker_list.append(dict(refmarkerdict))
            
######################################
#Creating the reference genotypes for each missing sample.
#e.g. [{'R20': 'C/C/C/C'}, {'R20': 'A/A/A/A'}]
master_ref_list=([{k: v} for k,v in zip(refsamplist, ref_list)])

######################################
#Combining the marker and genotype information for each missing sample.
#e.g. [{'marker': 'Ro06_244736', 'R20': 'C/C/C/C'}, {'marker': 'Ro06_244742', 'R20': 'A/A/A/A'}]
ref_combine_list =[{**marker, **ref} for marker, ref in zip(new_marker_list, master_ref_list)]

######################################
#Getting multiple sample information per marker. Some markers might have two or more samples that need genotype inforamtion. 
#e.g. {'Ro06_244736': {'marker': 'Ro06_244736', 'R20': 'C/C/C/C'}, 'Ro06_244742': {'marker': 'Ro06_244742', 'R20': 'A/A/A/A'}}

d_updated = {}
for i in ref_combine_list:
    d_updated.setdefault(i['marker'], {}).update(i)

######################################
#Replacing and putting in reference genotype information if the sample had read depth present.
#before e.g.'Ro06_244736': ['A/A/A/A', 'R37', 'R20'], 'Ro06_244742': ['G/G/G/G', 'R37', 'R20'] 
#after e.g. {'Ro06_244736': ['A/A/A/A', 'R37', 'C/C/C/C'], 'Ro06_244742': ['G/G/G/G', 'R37', 'A/A/A/A']}

d1_updated = {k: [d_updated[k].get(i, i) for i in l] for k, l in d1.items()}

######################################
#Outputting genotype table with depth and reference genotype taken into account. 

out_f.write("marker"+"\t"+"\t".join(uniq_samps)+"\n")
for k6, v6 in d1_updated.items():
    out_f.write(k6+"\t"+"\t".join(v6)+"\n")
 