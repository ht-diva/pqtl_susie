
from pathlib import Path
import pandas as pd


# Define input for the rules
# read loci list
lb = pd.read_csv(config["path_lb"])

# Create a new column by concatenating 
lb["locus"]  = lb["chr"].astype(str) + "_" + lb["start"].astype(str) + "_" + lb["end"].astype(str)
lb["locuseq"] = lb["seqid"].astype(str) + "_" + lb["locus"].astype(str)

data = (
    pd.DataFrame(lb, columns=["locuseq", "seqid", "chr", "locus", "SNPID", "gwas_path"])
    .set_index("locuseq", drop=False)
    .sort_index()
)


def ws_path(file_path):
    return str(Path(config.get("workspace_path"), file_path))

# return locus of locuseq
def get_locus(wildcards):
    return str(data.loc[wildcards, "locus"])

# return GWAS summary results 
def get_gwas(wildcards):
    seqid = data.loc[wildcards, "seqid"]
    #file_path = f"{seqid}/{seqid}.gwaslab.tsv.bgz"
    file_path = data.loc[wildcards, "gwas_path"]
    #return str(Path(config.get("path_gwas"), file_path))
    return str(Path(file_path))

# return genotype
def get_geno(wildcards):
    chrom = data.loc[wildcards, "chr"]
    path = config.get("genotype")
    filename = f"{path}{chrom}.pgen"
    return str(Path(filename))
