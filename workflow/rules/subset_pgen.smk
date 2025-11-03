
rule subset_pgen:
    input:
        pgen = lambda wildcards: get_geno(wildcards.locuseq),
        snplist = rules.subset_gwas.output.snplist
    output:
        pgen=temp(ws_path("tmp/{locuseq}_dosage.pgen")),
        pvar=temp(ws_path("tmp/{locuseq}_dosage.pvar")),
        psam=temp(ws_path("tmp/{locuseq}_dosage.psam")),
    params:
        genotype=lambda wildcards, input: input.pgen.replace(".pgen", ""),
        ofile= lambda wildcards, output: output.pgen.replace(".pgen", ""),
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
