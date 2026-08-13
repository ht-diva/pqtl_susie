
from pathlib import Path
import pandas as pd
import os
import math

# Define input for the rules
isutdy = config["study"]
igeno = config.get("genotype")
gdir  = config.get("path_gwas")
ld_src = config["ld"]["source"]
ld_dir = config["ld"]["directory"]

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
def get_gwas(wc):
    seqid = data.loc[wc.locuseq, "seqid"]
    
    if istudy == "gnh":
        return f'{gdir}/{wc.locuseq}/{wc.locuseq}_genesandhealth_v010_quantitative_traits_median_values_f5ff31e8c6.csv.gz'
    #file_path = f"{seqid}/{seqid}.gwaslab.tsv.gz"
    #return str(Path(gdir, file_path))
    return f'{gdir}/{seqid}/{seqid}.gwaslab.tsv.gz'

# return genotype
def get_geno(wildcards):
    chrom = data.loc[wildcards, "chr"]
    filename = f"{igeno}{chrom}.pgen"
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
    multiplier = 5

    mem_gb = ld_size_gb * multiplier + 4   # add 4GB overhead for R
    mem_mb = int(mem_gb * 1024)

    # clamp to reasonable range
    return max(mem_mb, 4000)

# Accepted studies
# STUDY_LDFILE = {
#     "believe": "{chrom}_qced_new_id_alleles",
#     "interval": "",
#     "meta": "",
#     "Meta_Interval":,
#     "gnh":  
# }

def get_ld(wc):
    if ld_src == "external":
        locus = data.loc[wc.locuseq, "locus"]
        return f'{ld_dir}/{locus}_ld.matrix'
    return rules.compute_ld.output.ld

def get_header(wc):
    if ld_src == "external":
        locus = data.loc[wc.locuseq, "locus"]
        return f'{ld_dir}/{locus}_ld.header'
    return rules.compute_ld.output.headers

