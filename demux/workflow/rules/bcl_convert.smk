rule bcl_convert:
    input:
        bcl=in_dir,
        sample_sheet=sample_sheet,
    output:
        dir=directory(out_dir),
        reports=expand(
            out_dir / "Reports" / "{report}",
            report=[
                "RunCompletionStatus.xml",
                "RunInfo.xml",
                "RunParameters.xml",
                "SampleSheet.csv",
                "Demultiplex_Stats.csv",
                "Demultiplex_Tile_Stats.csv",
                "Index_Hopping_Counts.csv",
                "IndexMetricsOut.bin",
                "Top_Unknown_Barcodes.csv",
                "Adapter_Cycle_Metrics.csv",
                "Adapter_Metrics.csv",
                "Quality_Metrics.csv",
                "Quality_Tile_Metrics.csv",
            ],
        ),
    log:
        "logs/bcl_convert.log",
    threads: 4
    resources:
        mem="12 GiB",
    params:
        extra="--bcl-num-parallel-tiles 1 --bcl-sampleproject-subdirectories true --force",
    shell:
        "bcl-convert --bcl-num-conversion-threads {threads} --bcl-num-compression-threads 2 --bcl-num-decompression-threads 1 --bcl-input-directory {input.bcl} --sample-sheet {input.sample_sheet} {params.extra} --output-directory {output} > {log} 2>&1"


rule fastq_undetermined:
    input:
        rules.bcl_convert.output.dir,
    output:
        out_dir / "Undetermined_S0_L{Lane}_R{read}_001.fastq.gz",
    localrule: True


rule fastq_sample:
    input:
        rules.bcl_convert.output.dir,
    output:
        out_dir
        / "{Sample_Project}"
        / "{Sample_ID}_S{sample_n}_L{Lane}_R{read}_001.fastq.gz",
    localrule: True
