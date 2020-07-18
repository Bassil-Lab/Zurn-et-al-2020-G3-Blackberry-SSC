How to make genotype table from zurn et al. 2020. Genotype table is a tab seperated file with colums that conatin marker (chrom_position) and sample genotypes (A/A/A/A) for polyploids. The code takes into account missing samples and samples with reference genotypes.

Step 1. Create seperate VCF (variant call files) for each sample. File need to be in VCFv4 format. 

Step 2. Run make_sample_genotype_files.py on each VCF.

make_sample_genotype_files.py requires 3 argumrants:

-vcf --vcf_file This should be the path and name to your vcf file
-s --sample_name This should be a short name given to each sample. Each sample must have a unique sample    name.
-p --out_path This is a path for the output. The path needs to end in /. 

Example usage: /path/to/make_sample_genotype_files.py -vcf /path/to/vcf/file.vcf -s short_sample_name -p /path/to/output/

Step 3. Run make_combined_genotype_table.py

make_combined_genotype_table.py requires 4 arguments:

-g --sample_genotypes This the path and file name to a one column file that contains a text file that conatins the path and file name to all vcf files being used. All samples must be in the same order.
-d --sample_depths This the path and file name to a one column file that contains a text file that contains the path and file names to all the depth files being used. All samples must be in the same order.
-s --sample_names This the path and file name to a one column file that contains a text file that contains all the short sample names. All samples must be in the same order.
-o --out_file This is the path and file name to output.

Example useage: ./make_combined_genotype_table.py -g /path/to/sample_genotypes.txt -d /path/to/sample_depths.txt -s /path/to/sample_depths.txt -o /path/to/my/output_table.tsv

