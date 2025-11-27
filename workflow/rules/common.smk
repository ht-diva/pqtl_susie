
from pathlib import Path
import pandas as pd


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

# Estimate memory needs for a SuSiE RSS job based on number of SNPs
def estimate_mem_mb(snp_list_file):
    """    
    Memory scales roughly with O(n^2) because SuSiE manipulates
    the LD matrix (double-precision float).
    
    This estimator computes memory for LD matrix: n^2 * 8 bytes
    """
    # Count SNPs
    with open(snp_list_file) as f:
        n_snps = sum(1 for _ in f)
    
    # LD matrix size (n*n doubles, 8 bytes each)
    ld_mb = (n_snps * n_snps * 8) / 1e6
    
    # SuSiE overhead + R session + buffers
    overhead_mb = 1024
    
    # Safety factor 25% overhead
    mem_estimate = ld_mb * 0.25 + overhead_mb
    
    # Clamp
    mem_estimate = max(2048, min(int(mem_estimate), 64000))
    
    return mem_estimate
