

rule exomedepth:
    input:
        bams=expand("/scratch/wp3/Exomedepth/{sex}/%s_N.bam" % sample for sample in get_samples(samples)),
        done="conifer/RPKM/Done.txt",
        ref=config["exomedepth"]["ref"]+{sex}RefCount.mat,
    output:
        count: "ExomeCount{sex}.txt",
        aed: "ExomeDepth_{sample}_{sex}.aed"),
        txt: "ExomeDepth_{sample}_{sex}.txt"),
        filtered: "ExomeDepth_{sample}_{sex}_filtered.aed"),
        sv: "ExomeDepth_{sample}_{sex}_SV.txt"),
    log:
        "exomedepth/{sample}.rpkm.log",
    benchmark:
        repeat(
            "exomedepth/{sample}.rpkm.benchmark.tsv",
            config.get("exomedepth", {}).get("benchmark_repeats", 1),
        )
    params:
        dir: config["programdir"]["dir"],
    threads: config.get("exomedepth", {}).get("threads", config["default_resources"]["threads"]),
    resources:
        mem_mb=config.get("exomedepth", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("exomedepth", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("exomedepth", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("exomedepth", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("exomedepth", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("exomedepth", {}).get("container", config["default_container"]),
    message:
        "{rule}: Generate exomedepth files for {wildcards.sample}",
    shell:
        "(Rscript {params.dir}/scripts/ExomeDepth{sex}.R {input.ref} {input.bams}) &> {log}"
