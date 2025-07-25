
from pathlib import Path
import pandas as pd


# Define input for the rules
# read loci list
lb = pd.read_csv(config["path_lb"])

# Create a new column by concatenating 
lb["locus"]  = lb["chr"].astype(str) + "_" + lb["start"].astype(str) + "_" + lb["end"].astype(str)
lb["locuseq"] = lb["phenotype_id"].astype(str) + "_" + lb["locus"].astype(str)

data = (
    pd.DataFrame(lb, columns=["locuseq", "phenotype_id", "chr", "locus", "SNPID"])
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
    seqid = data.loc[wildcards, "phenotype_id"]
    file_path = f"{seqid}/{seqid}.gwaslab.tsv.bgz"
    return str(Path(config.get("path_gwas"), file_path))

# return genotype
def get_geno(wildcards):
    chrom = data.loc[wildcards, "chr"]
    path = config.get("genotype")
    filename = f"{path}{chrom}.bed"
    return str(Path(filename))
