import pandas as pd
from concurrent.futures import ProcessPoolExecutor

def _process_chr(df_chr, func, **kwargs):
    return func(df_chr, **kwargs)

def parallel_split_and_run(df, func, threads=4, **kwargs):
    results = []
    with ProcessPoolExecutor(max_workers=threads) as ex:
        futures = {ex.submit(_process_chr, g, func, **kwargs): n for n, g in df.groupby("chr")}
        for f in futures:
            results.append(f.result())
    print("Parallel processing complete.")
    return pd.concat(results, ignore_index=True)
