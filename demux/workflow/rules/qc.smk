
rule qc_undetermined_top:
    input:
        out_dir / "Reports" / "Top_Unknown_Barcodes.csv",
    output:
        "stats/undetermined_top.tsv",
    params:
        extra="""filter '$index != "GGGGGGGGGG" && $index2 != "GGGGGGGGGG" && $index != "NNNNNNNNNN" && $index2 != "NNNNNNNNNN" && ${# Reads} > 5e6""",
    wrapper:
        "v9.15.0/utils/miller"


rule qc_cross_contamination:
    input:
        idx_counts=out_dir / "Reports" / "Index_Hopping_Counts.csv",
        idx_list=Path(workflow.basedir)
        / ".."
        / "resources"
        / "eDNA_index_list_UDP097-UDP288_UDI001-UDI096_250807.txt",
    output:
        directory(out_dir / "Reports" / "Index_Hopping_Counts" / "{pool}"),
    log:
        "logs/qc_cross_contamination/{pool}.log",
    conda:
        "../envs/cross_contamination.yaml"
    params:
        lanes=lambda w: ",".join(ss_pools[w.pool]),
        extra="--rpm-warn 100",
    shell:
        "./cross_contamination.py --index-counts {input.idx_counts} --index-known {input.idx_list} --lanes {params.lanes} {params.extra} --out-prefix {output} > {log}"


rule qc_md5sum:
    input:
        rules.fastq_sample.output,
    output:
        temp("temp/qc/md5sum/{Sample_Project}/{Sample_ID}_S{sample_n}_L{Lane}_R{read}.md5"),
    threads: 1
    shell:
        "md5sum {input} > {output}"


rule qc_md5sum_cat:
    input:
        # expand(rules.qc_md5sum.output, )
        expand("temp/qc/md5sum/{item.Sample_Project}/{item.Sample_ID}_S{item.sample_n}_L{item.Lane}_R{item.read}.md5",
            item=lookup(query="Sample_Project == '{Project_ID}'", within=ss_samples),
        ),
    output:
        out_dir / "{Project_ID}" / f"{run_id}.md5",
    threads: 5
    shell:
        "cat {input} | tee {output} | cut --delimiter --fields 1 | sort | uniq -d"
