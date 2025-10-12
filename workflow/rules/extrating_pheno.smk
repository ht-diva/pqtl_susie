
rule grep_pheno:
    message:
        "Extract protein residuals from the original GWAS"
    input:
        phenotype = config["phenotype"],
    output:
        ws_path("tmp/{locuseq}_pheno.csv"),
    resources:
        runtime=lambda wc, attempt: attempt * 10,
    conda:
        "../envs/pandas.yml"
    script:
        "../scripts/grep_pheno.py"
