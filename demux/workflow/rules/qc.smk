
rule qc_undetermined_nreads:
    input:
        "<out_dir>/{run_id}/Reports/Top_Unknown_Barcodes.csv",
    output:
        temp("temp/{run_id}/qc/undetermined/nreads.tsv"),
    log:
        "logs/{run_id}/qc/undetermined/nreads.log",
    params:
        extra="""filter -e '$index != "GGGGGGGGGG" && $index2 != "GGGGGGGGGG" && $index != "NNNNNNNNNN" && $index2 != "NNNNNNNNNN" && ${# Reads} > 5e6'""",
    wrapper:
        "v9.15.0/utils/miller"


use rule qc_undetermined_nreads as qc_fastq_nreads with:
    input:
        "<out_dir>/{run_id}/Reports/Demultiplex_Stats.csv",
    output:
        temp("temp/{run_id}/qc/fastq/nreads.tsv"),
    log:
        "logs/{run_id}/qc/fastq/nreads.log",
    params:
        extra="""filter -e '${# Reads} < 5e6'""",


rule qc_cross_contamination:
    input:
        idx_counts="<out_dir>/{run_id}/Reports/Index_Hopping_Counts.csv",
        idx_list=config["idx_list"],
    output:
        multiext(
            "<out_dir>/{run_id}/Reports/Index_Hopping_Counts/{pool}",
            ".counts.tsv",
            ".counts.html",
            ".cross_contam.tsv",
            ".cross_contam.html",
        ),
    log:
        stdout="logs/{run_id}/qc/cross_contamination/{pool}.log",
        stderr="logs/{run_id}/qc/cross_contamination/{pool}.err",
    conda:
        Path(workflow.basedir) / "scripts/cross_contamination/environment.yaml"
    params:
        lanes=lambda w: ",".join(ss_pools[w.pool]),
        outdir=lambda w, output: Path(output[0]).with_suffix("").with_suffix(""),
        extra="--rpm-warn 100",
    shell:
        "{workflow.basedir}/scripts/cross_contamination/cross_contamination.py --index-counts {input.idx_counts} --index-known {input.idx_list} --lanes {params.lanes} {params.extra} --out-prefix {params.outdir} > {log.stdout} 2> {log.stderr}"


rule qc_fastq_md5sum:
    input:
        fastq_sample,
    output:
        temp(
            "temp/{run_id}/qc/fastq/md5/{Sample_Project}/{Sample_ID}_S{sample_n}_L00{Lane}_R{read}.md5"
        ),
    shell:
        "md5sum {input} > {output}"


rule qc_fastq_md5sum_cat:
    input:
        # expand(rules.qc_fastq_md5sum.output, )
        expand(
            "temp/{run_id}/qc/fastq/md5/{item.Sample_Project}/{item.Sample_ID}_S{item.sample_n}_L00{item.Lane}_R{item.read}.md5",
            item=lookup(
                query="Sample_Project == '{Sample_Project}'", within=ss_samples
            ),
            allow_missing=True,
        ),
    output:
        md5sum="<out_dir>/{run_id}/{Sample_Project}/{run_id}.md5",
    log:
        "logs/{run_id}/qc/fastq/md5sum/{Sample_Project}.log",
    message:
        "Concatenating MD5SUM files."
    shell:
        r"cat {input} | sed 's:\S*/::g' | tee {output.md5sum} | cut --delimiter ' ' --fields 1 | sort | uniq -d | sed '1i # Duplicated MD5 checksums' > {log}"
