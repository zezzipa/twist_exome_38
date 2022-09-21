rule conifer_mapq20:
    input:
        bam="parabricks/pbrun_fq2bam/{sample}_N.bam",
        bai="parabricks/pbrun_fq2bam/{sample}_N.bam.bai",
    output:
        temp("conifer/MAPQ20/{sample}.MAPQ20.bam"),
    resources:
        mem_mb=config.get("conifer", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("conifer", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("conifer", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("conifer", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get(config["default_container"]),
    shell:
        "samtools view -bq 20 {input.bam} > {output}"


rule conifer_mapq20_index:
    input:
        "conifer/MAPQ20/{sample}.MAPQ20.bam",
    output:
        temp("conifer/MAPQ20/{sample}.MAPQ20.bam.bai"),
    resources:
        mem_mb=config.get("conifer", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("conifer", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("conifer", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("conifer", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get(config["default_container"]),
    shell:
        "samtools index {input}"


rule conifer_rpkm:
    input:
        bam="conifer/MAPQ20/{sample}.MAPQ20.bam",
        bai="conifer/MAPQ20/{sample}.MAPQ20.bam.bai",
        ref=config["reference"]["conifer"],
    output:
        RPKM=temp("conifer/RPKM/{sample}.rpkm"),
    log:
        "conifer/RPKM/{sample}.rpkm.log",
    benchmark:
        repeat(
            "conifer/RPKM/{sample}.rpkm.benchmark.tsv",
            config.get("conifer", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
    params:
        dir="/beegfs-storage/projects/wp3/nobackup/Workspace/CoNIFER",
    resources:
        mem_mb=config.get("conifer", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("conifer", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("conifer", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("conifer", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("conifer", {}).get("container", config["default_container"]),
    message:
        "{rule}: Generate conifer rpkm files for {wildcards.sample}",
    shell:
        "python {params.dir}/conifer.py rpkm "
        "--probes {input.ref} "
        "--input {input.bam} "
        "--output {output.RPKM} &> {log}"


rule conifer_sex:
    input:
        rpkm=expand("conifer/RPKM/%s.rpkm" % sample for sample in get_samples(samples)),
        idx="qc/multiqc/multiqc_DNA.html",
    output:
        done="conifer/RPKM/Done.txt",
    log:
        "conifer/RPKM/divide_sexes.log",
    params:
        config["programdir"]["dir"],
    resources:
        mem_mb=config.get("conifer", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("conifer", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("conifer", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("conifer", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get(config["default_container"]),
    shell:
        "python {params}/scripts/divide_sexes.py &> {log}"


rule exomedepth_sex:
    input:
        bam="parabricks/pbrun_fq2bam/{sample}_N.bam",
        bai="parabricks/pbrun_fq2bam/{sample}_N.bam.bai",
        done="conifer/RPKM/Done.txt",
    output:
        bam="cnv_sv/exomedepth_{sex}/{sample}_N.bam",
    log:
        "cnv_sv/exomedepth_{sex}/{sample}_moved.log",
    params:
        config["programdir"]["dir"],
    resources:
        mem_mb=config.get("conifer", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("conifer", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("conifer", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("conifer", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get(config["default_container"]),
    script:
        "{params}/scripts/move_ed.py"


rule conifer_analyze:
    input:
        done="conifer/RPKM/Done.txt",
        dir=config["conifer"]["dir"]+"RPKM_{sex}/",
        ref=config["reference"]["conifer"],
    output:
        hdf5=temp("conifer/SVD-ZRPKM/analyse_{sex}_svd6.hdf5"),
        svtxt=temp("conifer/SVD-ZRPKM/singular_{sex}_values.txt"),
        sdtxt=temp("conifer/SVD-ZRPKM/sd_{sex}_values.txt"),
    log:
        "conifer/SVD-ZRPKM/conifer_{sex}_analyze.log",
    benchmark:
        repeat(
            "conifer/SVD-ZRPKM/conifer_{sex}_analyze.benchmark.tsv",
            config.get("conifer", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
    params:
        dir="/beegfs-storage/projects/wp3/nobackup/Workspace/CoNIFER",
    resources:
        mem_mb=config.get("conifer", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("conifer", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("conifer", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("conifer", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("conifer", {}).get("container", config["default_container"]),
    message:
        "{rule}: Analysing conifer rpkm files",
    shell:
        "python {params.dir}/conifer.py analyze --rpkm_dir {input.dir} --svd 6 "
        "--probes {input.ref} "
        "--write_svals {output.svtxt} "
        "--write_sd {output.sdtxt} "
        "--output {output.hdf5} &> {log}"

# Can add
# png=temp("conifer/SVD-ZRPKM/screenplot.png"),  in output and
# "--plot_scree {output.png} "  in shell
# but matplotlib and pylab is needed

# To make screenplot of region to check
# python conifer.py plot --input SVD-ZRPKM/analyse.hdf5 --region chr1:1565000-1670000
# --output image1.b.png --sample D20-05543-ready_8


rule conifer_call:
    input:
        "conifer/SVD-ZRPKM/analyse_{sex}_svd6.hdf5",
    output:
        "conifer/calls/calls_{sex}_svd6.txt",
    log:
        "conifer/calls/conifer_{sex}_call.log",
    benchmark:
        repeat(
            "conifer/calls/conifer_{sex}_call.benchmark.tsv",
            config.get("conifer", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
    params:
        dir="/beegfs-storage/projects/wp3/nobackup/Workspace/CoNIFER",
    resources:
        mem_mb=config.get("conifer", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("conifer", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("conifer", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("conifer", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("conifer", {}).get("container", config["default_container"]),
    message:
        "{rule}: Calling CNVs from conifer data",
    shell:
        "python {params.dir}/conifer.py call --threshold 1.75 "
        "--input {input} "
        "--output {output} &> {log}"


rule conifer_individual:
    input:
        F="conifer/calls/calls_F_svd6.txt",
        M="conifer/calls/calls_M_svd6.txt",
    output:
        "conifer/samples/{sample}.txt",
    log:
        "conifer/samples/{sample}.log",
    benchmark:
        repeat(
            "conifer/samples/{sample}.benchmark.tsv",
            config.get("conifer", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
    params:
        dir="/beegfs-storage/projects/wp3/nobackup/Workspace/CoNIFER",
    resources:
        mem_mb=config.get("conifer", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("conifer", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("conifer", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("conifer", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("conifer", {}).get("container", config["default_container"]),
    message:
        "{rule}: Get calls for each sample",
    shell:
        """awk "BEGIN{{ print "sample\\tchr\\tstart\\tstop\\tcnv"}}" > {output} | grep {wildcards.sample} {input} >> {output}"""


rule conifer_aed:
    input:
        txtfile="conifer/samples/{sample}.txt",
    output:
        aedfile="conifer/samples/{sample}.aed",
    log:
        log="conifer/samples/{sample}_aed.log",
    benchmark:
        repeat(
            "conifer/samples/{sample}_aed.benchmark.tsv",
            config.get("conifer", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
    params:
        config["programdir"]["dir"],
    resources:
        mem_mb=config.get("conifer", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("conifer", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("conifer", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("conifer", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("conifer", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("exomedepth", {}).get("container", config["default_container"])
    message:
        "{rule}: Get calls for each sample",
    script:
        "/projects/wp3/nobackup/Workspace/Jessika_playground/fastq_med_cnv/twist_exome_hg37/workflow/scripts/conifer_aed.R"
