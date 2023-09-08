import pandas as pd
import statsmodels.formula.api as smf
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("data/GSM_final_genes.csv.gz")

def replace_dot_with_underscore(df):
    df.rename(columns=lambda x: x.replace(".", "_"), inplace=True)
    df.rename(columns=lambda x: x.replace("-", "_"), inplace=True)
    df.rename(columns=lambda x: x.replace("/", "_"), inplace=True)
    df.rename(columns=lambda x: x.replace("@", " "), inplace=True)
    return df

# Specify covariates 
def fit_mixedlm(df, gene):
    formula = f"{gene} ~ AD + tissue + sex + age"
    endog = df[gene]
    groups = df["series"]
    md = smf.mixedlm(formula, data=df, groups=groups)
    return md

# Fit linear mixed model wiht lbfgs algorithm
def get_log2fc_pval(df, gene):
    md = fit_mixedlm(df, gene)
    mdf = md.fit(method='lbfgs')
    pval_final_filt = mdf.pvalues["AD"]
    log2fc_final_filt = mdf.params["AD"]
    return log2fc_final_filt, -np.log10(pval_final_filt)

df = replace_dot_with_underscore(df)
# First 6 columns were metadata
gene_columns = df.columns[6:].tolist()  
log2fc_pval_final_filt = [get_log2fc_pval(df, gene) for gene in gene_columns]
log2fc_final_filt, pval_final_filt = zip(*log2fc_pval_final_filt)

# Volcano plot
final_scatter = sns.scatterplot(x=log2fc_final_filt, y=pval_final_filt)
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 p value')
plt.show()

