
ruleorder: bcl_convert > fastq_sample


rule bcl_convert:
    input:
        bcl=in_dir,
        sample_sheet=sample_sheet,
    output:
        # Not created by current version of BCL-convert
        touch(out_dir / "{run_id}/Reports/Demultiplex_Detailed_Stats.csv"),
        touch(out_dir / "{run_id}/Reports/RunParameters.xml"),
        touch(out_dir / "{run_id}/Reports/RunCompletionStatus.xml"),
        # Created by BCL-convert
        fq_samples=[
            expand(
                out_dir
                / "{run_id}/{Sample_Project}/{Sample_ID}_S{sample_n}_L00{Lane}_R{read}_001.fastq.gz",
                **sample,
                allow_missing=True,
            )
            for sample in ss_samples[
                ss_samples["Sample_Project"].ne("Undetermined")
            ].to_dict("records")
        ],
        fq_undetermined=expand(
            out_dir / "{run_id}/Undetermined_S0_L00{Lane}_R{read}_001.fastq.gz",
            Lane=ss_lanes,
            read=ss_reads,
            allow_missing=True,
        ),
        reports=expand(
            out_dir / "{run_id}/Reports/{report}",
            report=[
                "RunInfo.xml",
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
                "fastq_list.csv",
            ],
            allow_missing=True,
        ),
        logs=expand(
            out_dir / "{run_id}/Logs/{log}",
            log=["Errors.log", "FastqComplete.txt", "Info.log", "Warnings.log"],
            allow_missing=True,
        ),
    log:
        "logs/{run_id}/bcl_convert.log",
    envmodules:
        "bcl-convert/4.4.6",
    threads: 4
    resources:
        mem="12 GiB",
    params:
        outdir=lambda w, output: Path(output.fq_undetermined[0]).parent,
        extra="--bcl-num-parallel-tiles 1 --bcl-sampleproject-subdirectories true --force",
    shell:
        "bcl-convert --bcl-num-conversion-threads {threads} --bcl-num-compression-threads 2 --bcl-num-decompression-threads 1 --bcl-input-directory {input.bcl} --sample-sheet {input.sample_sheet} {params.extra} --output-directory {params.outdir} > {log} 2>&1"


rule fastq_sample:
    input:
        out_dir / "{run_id}/{Sample_ID}_S{sample_n}_L00{Lane}_R{read}_001.fastq.gz",
    output:
        out_dir
        / "{run_id}/{Sample_Project}/{Sample_ID}_S{sample_n}_L00{Lane}_R{read}_001.fastq.gz",
    wildcard_constraints:
        Sample_Project="Undetermined",
    localrule: True
    shell:
        "ln {input} {output}"
