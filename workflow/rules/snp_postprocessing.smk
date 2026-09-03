# Module I post-processing


SNP_OUT = config["snp"]["output_dir"]
SAMPLES = config["snp"]["se_samples"] + config["snp"]["pe_samples"]
SE_SAMPLES = config["snp"]["se_samples"]
PE_SAMPLES = config["snp"]["pe_samples"]


rule snp_merge_se_pe:

    input:
        se=expand(rules.snp_gvcf_r3.output, sample=SE_SAMPLES),
        pe=expand(rules.snp_gvcf_r3.output, sample=PE_SAMPLES),
    output:
        f"{SNP_OUT}/merged_all_samples.vcf.gz"
    conda:
        "../envs/integration/bcftools.yaml"
    shell:
        """
        bcftools merge {input.se} {input.pe} -Oz -o {output}
        gatk IndexFeatureFile -I {output}
        """


rule snp_variant_filtration:
    input:
        rules.snp_merge_se_pe.output
    output:
        f"{SNP_OUT}/merged_filtered.vcf.gz"
    params:
        ref=config["snp"]["reference"]
    conda:
        "../envs/integration/gatk.yaml"
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input} -O {output} \
            --filter-expression "QD < 2.0" --filter-name "QD2" \
            --filter-expression "FS > 60.0" --filter-name "FS60" \
            --filter-expression "MQ < 40.0" --filter-name "MQ40" \
            --filter-expression "SOR > 3.0" --filter-name "SOR3" \
            --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
        """


rule snp_extract_pass:
    input:
        rules.snp_variant_filtration.output
    output:
        vcf=f"{SNP_OUT}/merged_PASS.vcf.gz"
    conda:
        "../envs/integration/bcftools.yaml"
    shell:
        "bcftools view -f PASS {input} -Oz -o {output.vcf}"


rule snp_variant_stats:
    input:
        rules.snp_extract_pass.output.vcf
    output:
        f"{SNP_OUT}/variant_stats.txt"
    conda:
        "../envs/integration/bcftools.yaml"
    shell:
        "bcftools stats {input} > {output}"


rule snp_only_extract:
    input:
        rules.snp_extract_pass.output.vcf
    output:
        vcf=f"{SNP_OUT}/SNP_only.vcf.gz"
    conda:
        "../envs/integration/bcftools.yaml"
    shell:
        """
        bcftools view -v snps {input} -Oz -o {output.vcf}
        bcftools index {output.vcf}
        """


rule snp_chromosome_distribution:
    input:
        rules.snp_only_extract.output.vcf
    output:
        f"{SNP_OUT}/chromosome_snp_final.txt"
    conda:
        "../envs/integration/bcftools.yaml"
    shell:
        """
        bcftools query -f '%CHROM\\n' {input} | sort | uniq -c \
          | grep "NC_" | awk '$1>1000' > {output}
        """


rule snp_tstv:
    input:
        rules.snp_only_extract.output.vcf
    output:
        f"{SNP_OUT}/tstv_stats.txt"
    conda:
        "../envs/integration/bcftools.yaml"
    shell:
        'bcftools stats {input} | grep "^TSTV" > {output}'


rule snp_snpeff_annotate:
    input:
        rules.snp_only_extract.output.vcf
    output:
        f"{SNP_OUT}/SNP_annotated.vcf"
    params:
        genome=config["snpeff"]["genome_name"],
        snpeff_config=config["snpeff"]["config_path"],
    conda:
        "../envs/integration/snpeff.yaml"
    shell:
        "snpEff -c {params.snpeff_config} -noStats -v {params.genome} {input} > {output}"


rule snp_annotation_summary:
    input:
        rules.snp_snpeff_annotate.output
    output:
        f"{SNP_OUT}/SNP_annotation_statistics.txt"
    shell:
        r"""
        {{
          for eff in missense_variant synonymous_variant stop_gained stop_lost \
                     upstream_gene_variant downstream_gene_variant intergenic_region \
                     intron_variant 5_prime_UTR_variant 3_prime_UTR_variant; do
            echo -n "$eff: "; grep -o "$eff" {input} | wc -l
          done
          for impact in HIGH MODERATE LOW MODIFIER; do
            echo -n "$impact impact: "; grep -o "|$impact|" {input} | wc -l
          done
        }} > {output}
        """


rule snp_extract_high_impact:
    input:
        rules.snp_snpeff_annotate.output
    output:
        vcf=f"{SNP_OUT}/HIGH_impact_variants.vcf",
        transcripts=f"{SNP_OUT}/HIGH_only_transcripts.txt",
    shell:
        """
        grep "|HIGH|" {input} > {output.vcf}
        awk -F'|' '$3=="HIGH"{{print $7}}' {input} | sort -u > {output.transcripts}
        """


rule snp_high_impact_proteins:
    input:
        transcripts=rules.snp_extract_high_impact.output.transcripts,
        gff=config["snp"]["annotation_gff"],
        proteome=config["snp"]["proteome_fasta"],
    output:
        ids=f"{SNP_OUT}/HIGH_impact_proteins.ids",
        faa=f"{SNP_OUT}/HIGH_impact_proteins.faa",
    conda:
        "../envs/integration/seqtk.yaml"
    shell:
        """
        awk -F'\\t' '$3=="mRNA"' {input.gff} \
          | grep -oE 'Genbank:[A-Z]{{2}}_[0-9]+\\.[0-9]+,Genbank:[A-Z]{{2}}_[0-9]+\\.[0-9]+' \
          > results/snp/transcript_protein_map.txt
        grep -Ff {input.transcripts} results/snp/transcript_protein_map.txt \
          | grep -oE '[A-Z]{{2}}_[0-9]+\\.[0-9]+' | sort -u > {output.ids}
        seqtk subseq {input.proteome} {output.ids} > {output.faa}
        """


rule snp_high_impact_interproscan:
    input:
        rules.snp_high_impact_proteins.output.faa
    output:
        f"{SNP_OUT}/HIGH_impact_proteins_interpro.tsv"
    conda:
        "../envs/integration/interproscan.yaml"
    shell:
        "interproscan.sh -i {input} -f tsv -goterms -pa -o {output}"
