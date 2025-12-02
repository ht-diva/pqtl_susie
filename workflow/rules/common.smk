
from pathlib import Path
import pandas as pd
import os
import math

# Define input for the rules
# read loci list
lb = pd.read_csv(config["path_lb"])

# Create a new column by concatenating 
lb["locus"]  = lb["chr"].astype(str) + "_" + lb["start"].astype(str) + "_" + lb["end"].astype(str)
lb["locuseq"] = lb["seqid"].astype(str) + "_" + lb["locus"].astype(str)

data = (
    pd.DataFrame(lb, columns=["locuseq", "seqid", "chr", "locus", "SNPID"])
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
    file_path = f"{seqid}/{seqid}.gwaslab.tsv.gz"
    return str(Path(config.get("path_gwas"), file_path))

# return genotype
def get_geno(wildcards):
    chrom = data.loc[wildcards, "chr"]
    path = config.get("genotype")
    filename = f"{path}{chrom}.pgen"
    return str(Path(filename))

# Estimate memory needs for a SuSiE RSS job 
# based on the actual LD file size on disk.
def estimate_mem_mb(ld_file):
    """
    Empirical observation:
        - SuSiE uses ~2.5–3.5 × the file size in peak RAM
        - For large regions (>20k SNPs), use ×3.5 for safety
    """
    ld_size_bytes = os.path.getsize(ld_file)
    ld_size_gb = ld_size_bytes / 1e9

    # empirical multiplier:
    # R duplicates objects, SuSiE makes working copies, GC adds overhead
    multiplier = 7

    mem_gb = ld_size_gb * multiplier + 4   # add 4GB overhead for R
    mem_mb = int(mem_gb * 1024)

    # clamp to reasonable range
    return max(mem_mb, 4000)
