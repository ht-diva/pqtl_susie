
rule subset_pgen:
    input:
        pgen = lambda wildcards: get_geno(wildcards.locuseq),
        snplist = ws_path("tmp/{locuseq}_snps.list"), #rules.subset_gwas.output.snplist
    output:
        dosage=ws_path("tmp/{locuseq}_dosage.pgen"),
    params:
        genotype=lambda wildcards, input: input.pgen.replace(".pgen", ""),
        ofile=lambda wildcards, output: output.dosage.replace(".pgen", ""),
    resources:
        runtime=lambda wc, attempt: 30 + attempt * 10,
    shell:
        """
    source /exchange/healthds/singularity_functions

    plink2 --pfile {params.genotype} \
    --keep-allele-order \
    --extract {input.snplist} \
    --make-pgen \
    --out {params.ofile} \
    --memory 6000
        """
