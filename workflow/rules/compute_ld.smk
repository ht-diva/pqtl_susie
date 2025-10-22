
rule compute_ld:
    message:
        "Compute LD matrix using Plink2 with subset genotype in PGEN."
    input:
        # dosage  = rules.extract_dosage.output.dosage,
        dosage=ws_path("tmp/{locuseq}_dosage.pgen")
    output:
        ld=ws_path("tmp/{locuseq}_ld.matrix"),
        headers=ws_path("tmp/{locuseq}_ld.headers"),
    params:
        dosage=lambda wildcards, input, output: input.dosage.replace(".pgen", ""),
        prefix=lambda wildcards, input, output: output.ld.replace(".matrix", ""),
    # conda:
    #     "../envs/plink2.yml"
    resources:
        runtime=lambda wc, attempt: attempt * 30,
    shell:
        """
        source /exchange/healthds/singularity_functions

        plink2  \
          --pfile {params.dosage} \
          --r-unphased 'square' 'ref-based' \
          --out {params.prefix} \
          --threads 1  \
          --memory 2000
        
        mv {params.prefix}.unphased.vcor1  {output.ld}
        mv {params.prefix}.unphased.vcor1.vars  {output.headers}
        """