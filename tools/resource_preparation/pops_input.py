import os
import subprocess



## pops

    # https://github.com/FinucaneLab/pops

## feature file downloaded from here
    # https://www.dropbox.com/scl/fo/ne7xhxkt4dwhvd52a59ub/AFKkJu7ACaun1uuE99kmTkc/data/PoPS.features.txt.gz?rlkey=ltdbcld1enyr1zefg1lfqm61i&e=1&dl=0

##  gene_annot_jun10.txt downloaded from https://github.com/FinucaneLab/pops/tree/master/example/data/utils

# Make sure output directory exists
os.makedirs("/Users/JJOHN41/Documents/software_resources/resourses/postgwas/pops/features_munged/", exist_ok=True)

# Command arguments as list (BEST PRACTICE)
cmd = [
    "python",
    "/Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/pops/munge_feature_directory.py",
    "--gene_annot_path", "/Users/JJOHN41/Documents/software_resources/resourses/postgwas/pops/GRCh37_gene_annot_jun10.txt",
    "--feature_dir", "/Users/JJOHN41/Documents/software_resources/resourses/postgwas/pops/features_raw/",
    "--save_prefix", "/Users/JJOHN41/Documents/software_resources/resourses/postgwas/pops/features_munged/pops_features",
    "--max_cols", "500",
]

subprocess.run(cmd, check=True)


df=pd.read_csv('GRCh37_gene_annot_jun10.txt',sep="\t")

df[['ENSGID','CHR', 'START', 'END','NAME']].to_csv('GRCh37_gene_annot_jun10.loc',sep="\t",index=None)