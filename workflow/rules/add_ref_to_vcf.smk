
rule addRef:
    input:
        vcf="parabricks/deepvariant/{sample}.vcf",
        ref=config["reference"]["fasta"],
    output:
        temp("vcf_final/{sample}_ref.vcf"),
    log:
        "vcf_final/{sample}_add_ref.log",
    params:
        config["programdir"]["dir"],
    conda:
        "../envs/parabricks.yaml"
    shell:
        "( python {params}/scripts/ref_vcf.py {input.vcf} {input.ref} {output} ) &> {log}"


rule changeM2MT:
    input:
        "vcf_final/{sample}_ref.vcf",
    output:
        temp("vcf_final/{sample}.vcf"),
    log:
        "vcf_final/{sample}_chrMT.log",
    shell:
        """( awk '{{gsub(/chrM/,"chrMT"); print}}' {input} > {output} ) &> {log}"""


rule bgzipNtabix:
    input:
        "vcf_final/{sample}.vcf",
    output:
        "vcf_final/{sample}.vcf.gz",
        "vcf_final/{sample}.vcf.gz.tbi",
    log:
        "vcf_final/{sample}.bgzip-tabix.log",
    conda:
        "../envs/parabricks.yaml"
    shell:
        "( bgzip {input} && tabix {input}.gz ) &> {log}"
