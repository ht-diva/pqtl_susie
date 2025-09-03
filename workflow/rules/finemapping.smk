

rule subset_gwas:
    input:
        gwas = lambda wildcards: get_gwas(wildcards.locuseq),
    output:
        sumstat = ws_path("fm/tmp/{locuseq}_sumstat.csv"),
        snplist = ws_path("fm/tmp/{locuseq}_snps.list"),
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


rule extract_dosage:
    input:
        bed = lambda wildcards: get_geno(wildcards.locuseq),
        snplist = ws_path("fm/tmp/{locuseq}_snps.list"), #rules.subset_gwas.output.snplist
    output:
        dosage=ws_path("fm/tmp/{locuseq}_dosage.raw"),
    params:
        genotype=lambda wildcards, input: input.bed.replace(".bed", ""),
        ofile=lambda wildcards, output: output.dosage.replace(".raw", ""),
    resources:
        runtime=lambda wc, attempt: 30 + attempt * 10,
    shell:
        """
    source /exchange/healthds/singularity_functions

    plink --bfile {params.genotype} \
    --keep-allele-order \
    --extract {input.snplist} \
    --recode 'A' \
    --out {params.ofile} \
    --memory 6000
        """


rule run_susieR:
    input:
        dosage  = rules.extract_dosage.output.dosage,
        sumstat = ws_path("fm/tmp/{locuseq}_sumstat.csv"),
    output:
        cs_summary = ws_path("fm/cs/{locuseq}.cssum"),
        cs_rds  = ws_path("fm/cs/{locuseq}_fit.rds"),
        cs_list = ws_path("fm/cs/{locuseq}.cslist"),
    params:
        iter=config["susieR"]["iter"],
        min_abs_corr=config["susieR"]["min_abs_corr"],
        chrcol = config.get("sumstat").get("chrcol"),
    resources:
        runtime=lambda wc, attempt: attempt * 10,
    conda:
        "../envs/susier.yml"
    script:
        "../scripts/susie_best_practice.R"
