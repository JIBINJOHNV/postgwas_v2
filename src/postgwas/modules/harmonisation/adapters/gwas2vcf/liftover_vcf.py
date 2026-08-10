import os,glob,argparse,json
import sys
import pandas as pd
import numpy as np


bcftools="/usr/local/bin/bcftools"

# Define command-line arguments
parser = argparse.ArgumentParser(description="Summarize Mosedepth TargetRegion.thresholds output")

parser.add_argument('-input_path', '--input_path', help="input vcf file path", required=True)
parser.add_argument('-output_path', '--output_path', help="output vcf file path", required=True)
parser.add_argument('-input_vcf', '--input_vcf', help="input vcf file name", required=True)
parser.add_argument('-output_vcf', '--output_vcf', help="output vcf file name", required=True)
parser.add_argument('-fastafile', '--fastafile', help="name of the target fasta file with full path")
parser.add_argument('-chain_file', '--chain_file', help="chain file", required=True)


'''
input_path="/mnt/disks/sdd/BrainImage_CLusterMetaSumstat/blood_pressure_pmid38689001/"
output_path="/mnt/disks/sdd/BrainImage_CLusterMetaSumstat/blood_pressure_pmid38689001/"
input_vcf=f"{f_name_prefix}_build_GRCh38.vcf.gz"
output_vcf=f"{f_name_prefix}_build_GRCh37.vcf"
target_fasta="/home/jjohn41/Softwares/Resourses/Fastafiles/human_g1k_v37.fasta"
chain_file="/home/jjohn41/Softwares/Resourses/Fastafiles/GRCh38_to_GRCh37.chain.gz"
'''



args=parser.parse_args()

input_path=args.input_path
output_path=args.output_path
input_vcf=args.input_vcf
output_vcf=args.output_vcf
fastafile=args.fastafile
chain_file=args.chain_file


os.system(f'''CrossMap vcf {chain_file} {input_path}{input_vcf} {fastafile} {output_path}{output_vcf}''')
os.system(f'''{bcftools}  sort --output-type v --output {output_path}temp.vcf {output_path}{output_vcf}''')
os.system(f'''mv {output_path}temp.vcf {output_path}{output_vcf}''')
os.system(f"bgzip -f {output_path}{output_vcf} ")
os.system(f"tabix -f -p vcf {output_path}{output_vcf}.gz ")
