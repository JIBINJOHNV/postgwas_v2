


import pandas as pd 

path="~/Documents/software_resources/resourses/postgwas/imputation/pred-ld/ref/TOP_LD/EUR/SNV/"
for chromosome in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,22,'X']:
    df=pd.read_csv(f'EUR_chr{chromosome}_no_filter_0.2_1000000_info_annotation.csv.gz')
    df=df[['Position', 'rsID', 'MAF', 'REF', 'ALT', 'Uniq_ID']]
    df["rsID"] = ( str(chromosome)
        + "_" + df["Position"].astype(str)
        + "_" + df["REF"].astype(str)
        + "_" + df["ALT"].astype(str) )
    
    df.to_csv(
        f'new_folder/EUR_chr{chromosome}_no_filter_0.2_1000000_info_annotation.csv.gz',
        index=False,
        compression='gzip'
    )
