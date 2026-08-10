
rule qc_undetermined_adapters:
    input:
        out_dir / "{run_id}/Reports/Top_Unknown_Barcodes.csv",
    output:
        temp("temp/{run_id}/undetermined/adapters.tsv"),
    params:
        extra="""filter -e '$index != "GGGGGGGGGG" && $index2 != "GGGGGGGGGG" && $index != "NNNNNNNNNN" && $index2 != "NNNNNNNNNN" && ${# Reads} > 5e6'""",
    wrapper:
        "v9.15.0/utils/miller"


rule qc_cross_contamination:
    input:
        idx_counts=out_dir / "{run_id}/Reports/Index_Hopping_Counts.csv",
        idx_list=Path(workflow.basedir)
        / ".."
        / "resources"
        / "eDNA_index_list_UDP097-UDP288_UDI001-UDI096_250807.txt",
    output:
        multiext(
            str(out_dir / "{run_id}/Reports/Index_Hopping_Counts/{pool}"),
            ".counts.tsv",
            ".counts.html",
            ".cross_contam.tsv",
            ".cross_contam.html",
        ),
    log:
        "logs/{run_id}/qc_cross_contamination/{pool}.log",
    conda:
        Path(workflow.basedir) / ".." / "envs" / "cross_contamination.yaml"
    params:
        lanes=lambda w: ",".join(ss_pools[w.pool]),
        outdir=lambda w, output: Path(output[0]).with_suffix("").with_suffix(""),
        extra="--rpm-warn 100",
    shell:
        "{workflow.basedir}/scripts/cross_contamination.py --index-counts {input.idx_counts} --index-known {input.idx_list} --lanes {params.lanes} {params.extra} --out-prefix {params.outdir} 2> {log}"


rule qc_fastq_md5:
    input:
        rules.fastq_sample.output,
    output:
        temp(
            "temp/{run_id}/qc/fastq/md5/{Sample_Project}/{Sample_ID}_S{sample_n}_L00{Lane}_R{read}.md5"
        ),
    shell:
        "md5sum {input} > {output}"


rule qc_fastq_md5sum_cat:
    input:
        # expand(rules.qc_fastq_md5.output, )
        expand(
            "temp/{run_id}/qc/fastq/md5/{item.Sample_Project}/{item.Sample_ID}_S{item.sample_n}_L00{item.Lane}_R{item.read}.md5",
            item=lookup(
                query="Sample_Project == '{Sample_Project}'", within=ss_samples
            ),
            allow_missing=True,
        ),
    output:
        out_dir / "{run_id}/{Sample_Project}/{run_id}.md5",
    log:
        "logs/{run_id}/qc/fastq/md5sum/{Sample_Project}.log",
    shell:
        r"cat {input} | sed 's:\S*/::g' | tee {output} | cut --delimiter ' ' --fields 1 | sort | uniq -d > {log}"


rule qc_fastq_size:
    input:
        # expand(rules.fastq_sample.output, )
        expand(
            out_dir
            / "{run_id}/{item.Sample_Project}/{item.Sample_ID}_S{item.sample_n}_L00{item.Lane}_R{item.read}_001.fastq.gz",
            item=lookup(
                query="Sample_Project == '{Sample_Project}'", within=ss_samples
            ),
            allow_missing=True,
        ),
    output:
        temp("temp/{run_id}/qc/fastq/size/{Sample_Project}.tsv"),
    log:
        "logs/{run_id}/qc/fastq/size/{Sample_Project}.log",
    params:
        expr=lambda w: "$1 > 1e6" if w.Sample_Project == "Undetermined" else "$1 < 1e6",
    shell:
        "stat -c '%s %n' {input} | awk {params.expr:q} > {output} 2> {log}"
