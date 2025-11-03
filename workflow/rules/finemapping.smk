
rule run_susieR:
    input:
        pgen = rules.subset_pgen.output.pgen,
        pvar = rules.subset_pgen.output.pvar,
        psam = rules.subset_pgen.output.psam,
        sumstat = rules.subset_gwas.output.sumstat,
        ld = rules.compute_ld.output.ld,
        ld_snps = rules.compute_ld.output.headers,
    output:
        data_report = ws_path("susierss/cs_report/{locuseq}.report"),
        cs_summary = ws_path("susierss/cs_summary/{locuseq}.cssum"),
        cs_rds  = ws_path("susierss/cs_fitness/{locuseq}_fit.rds"),
        cs_annot= ws_path("susierss/cs_fitness/{locuseq}_annot.rds"),
        cs_list = ws_path("susierss/cs_list/{locuseq}.cslist"),
    log:
        ws_path("logs/susieR/{locuseq}.log"),
    params:
        iter=config["susieR"]["iter"],
        L=config["susieR"]["L"],
        min_abs_corr=config["susieR"]["min_abs_corr"],
        est_res_var =config["susieR"]["estimate_residual_variance"],
        chrcol = config.get("sumstat").get("chrcol"),
        ld_cor = config["run"]["ld_correlation"],
        study  = config["sumstat"]["study"],
        n_gwas = config["sumstat"]["n_samples"],
    resources:
        runtime=lambda wc, attempt: 6000 + attempt * 60,
    conda:
        "../envs/susier.yml"
    script:
        "../scripts/susie_best_practice.R"


rule collect_credible_sets:
    input:
        expand(ws_path("susierss/cs_list/{locuseq}.cslist"), locuseq = data.locuseq),
    output:
        ofile = ws_path("susierss/collected_credible_sets.tsv"),
    conda:
        "../envs/susier.yml"
    resources:
        runtime=lambda wc, attempt: 30 + attempt * 10,
    script:
        "../scripts/combine_cs_lists.R"
