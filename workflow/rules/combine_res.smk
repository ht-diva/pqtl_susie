
rule combine_results:
    input:
        expand(ws_path("susierss/cs_list/{locuseq}.cslist"), locuseq = data.locuseq),
    output:
        cslist = ws_path("susierss/combined_credibles.tsv"),
        report = ws_path("susierss/combined_reports.tsv"),
        sumstat = ws_path("susierss/combined_sumstats.tsv"),
    conda:
        "../envs/susier.yml"
    resources:
        runtime=lambda wc, attempt: 30 + attempt * 10,
    script:
        "../scripts/combine_cs_lists.R"
