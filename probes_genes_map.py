import pandas as pd
import wrenlab.normalize

# Quantile normalization + collapse from probes to genes

df = pd.read_csv("data/df_gsm.csv.gz", compression="gzip")
df.set_index('gsm', inplace=True)
df = df.T
df = wrenlab.normalize.quantile(df)
final_df = df.T

def fn():
    path = "data/GPL570.map.tsv"
    df = pd.read_table(path)
    o = {}
    for probe_id, symbols, gene_id in df.dropna().to_records(index=False):
        o[probe_id] = symbols.split(" /// ")[0]
    return o


M = fn()
final_df = final_df.groupby(M, axis=1).mean()
