


rule subset_gwas:
    input:
        gwas = lambda wildcards: get_gwas(wildcards.locuseq),
    output:
        sumstat = ws_path("tmp/{locuseq}_sumstat.csv"),
        snplist = ws_path("tmp/{locuseq}_snps.list"),
    params:
        locus = lambda wildcards: get_locus(wildcards.locuseq),
        prefix = "{locuseq}",
        tail = config.get("susieR").get("extension"),
    #conda:
    #    "envs/environment.yml"
    resources:
        runtime=lambda wc, attempt: 120 + attempt * 60,
    shell:
        """
        source /exchange/healthds/singularity_functions

        echo "Genomic region: {params.locus}"
       
        # take region bounaries from locus string
        chr=$(echo {params.locus} | cut -d'_' -f1)
        beg=$(echo {params.locus} | cut -d'_' -f2)
        end=$(echo {params.locus} | cut -d'_' -f3)
        
        # extend the boundaries by +/- 100 kbp
        beg_ext=$((beg - {params.tail}))
        end_ext=$((end + {params.tail}))
        
        # reformat locus to be readable for locuzoom
        region=${{chr}}:${{beg_ext}}-${{end_ext}}
        
        echo "Extended region is: $region"

        tabix  {input.gwas} $region -h > {output.sumstat}
        tabix  {input.gwas} $region -h | cut -f3 > {output.snplist}
        """

