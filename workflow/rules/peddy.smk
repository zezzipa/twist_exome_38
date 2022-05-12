
rule vcftools_vcfmerge:
    input:
        ["vcf_final/%s.vcf.gz" % sample for sample in get_samples(samples)],
    output:
        "peddy/all.vcf.gz",
    log:
        "peddy/vcftools_vcfmerge.log",
    shell:
        "vcf-merge {input} | bgzip -c > {output} &> {log}"


rule vcftools_plink:
    input:
        "peddy/all.vcf.gz",
    output:
        ped="peddy/all.ped",
        map="peddy/all.map",
    params:
        pre="peddy/all",
    log:
        "peddy/vcftools_plink.log",
    shell:
        "vcftools --gzvcf {input} --plink --out {params.pre} &> {log}"


rule peddy:
    input:
        vcf="peddy/all.vcf.gz",
        ped="peddy/all.ped",
    output:
        "peddy/all_peddy.peddy.ped"
    params:
        build=config.get("peddy", {}).get("p", ""),
        p=config.get("peddy", {}).get("build", "4"),
        pre="peddy/all_peddy",
    log:
        "peddy/peddy.log",
    shell:
        "peddy -p {params.p} {params.build} --plot --prefix {params.pre} {input.vcf} {input.ped} &> {log}"
