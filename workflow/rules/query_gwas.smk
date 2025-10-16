
rule create_query:
    output:
        qregion = ws_path("tmp/{locuseq}_region.tsv"),
        qyaml = ws_path("tmp/{locuseq}_query.yaml"),
    params:
        locus = lambda wildcards: get_locus(wildcards.locuseq),
        locuseq = "{locuseq}",
        tail = config.get("susieR").get("extension"),
    #conda:
    #    "envs/environment.yml"
    resources:
        runtime=lambda wc, attempt: 5 + attempt * 5,
    shell:
        """
        source /exchange/healthds/singularity_functions

        echo "Genomic region: {params.locus}"
       
        # take region bounaries from locus string
        chr=$(echo {params.locus} | cut -d'_' -f1)
        beg=$(echo {params.locus} | cut -d'_' -f2)
        end=$(echo {params.locus} | cut -d'_' -f3)
        
        # extend the boundaries if needed
        beg_ext=$((beg - {params.tail}))
        end_ext=$((end + {params.tail}))
        
        region=${{chr}}:${{beg_ext}}-${{end_ext}}
        echo "Extended region is: $region"

        # create region file
        echo "$chr\t$beg\t$end" > {output.qregion}

        # create yaml file
        
        # Extract "8280.238" then replace dot with dash
        seqid=$(echo {params.locuseq} | awk -F'[._]' '{{print $2"-"$3}}')
        
        # Write YAML file
cat <<EOF > {output.qyaml}
project: pqtl
study: believe

trait:
  - seqid: $seqid

output:
  - build
  - trait.desc
EOF

        echo "YAML file created: {output.qyaml}"
        """


rule query_gwas:
    input:
        qregion = ws_path("tmp/{locuseq}_region.tsv"),
        qyaml = ws_path("tmp/{locuseq}_query.yaml"),
    output:
        sumstat = ws_path("tmp/{locuseq}/{locuseq}_sumstat.csv.gz")
    #conda:
    #    "envs/environment.yml"
    params:
        prefix=lambda wildcards, output: output.sumstat.replace("_sumstat.csv.gz", ""),
        locuseq = "{locuseq}",
    resources:
        runtime=lambda wc, attempt: 120 + attempt * 60,
    shell:
        """
        source /exchange/healthds/singularity_functions

        gwasstudio export  --search-file {input.qyaml}  --get-regions {input.qregion}  --output-prefix {params.prefix}

        seqid=$(echo {params.locuseq} | cut -d'_' -f1)
        odir=$(dirname {output.sumstat})

        mv {params.prefix}_$seqid.csv.gz {output.sumstat}
        # mv slurm*.out $odir
        """
