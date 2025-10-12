import pandas as pd


def main(pheno_file, pheno, out_file):
    print(out_file)
    
    seqid = pheno.split("_")[0]
    print(f"Using seqid: {seqid}")

    pheno_df = pd.read_csv(pheno_file, sep="\t", header=0, usecols=["FID", "IID", seqid])
    pheno_df.rename(columns={seqid: "Y"}, inplace=True)
    pheno_df.to_csv(out_file, index=False, sep="\t")


if __name__ == "__main__":
    main(
        pheno_file = snakemake.input.phenotype,
        pheno      = snakemake.wildcards.locuseq,
        out_file   = snakemake.output[0]
        )
