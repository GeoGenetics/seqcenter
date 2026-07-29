### REPORTS
rule reports_run:
    input:
        in_dir / "{report}",
    output:
        out_dir / "Reports" / "{report}",
    threads: 1
    resources:
        mem="100 MiB",
    shell:
        "rsync -a {input} {output} > {log} 2>&1"


use rule reports_run as reports_demux with:
    input:
        in_analysis / "Demux" / "{report}",


use rule reports_run as reports_bclconvert with:
    input:
        in_analysis / "BCLConvert" / "fastq" / "Reports" / "{report}",


use rule reports_run as reports_fastqc with:
    input:
        in_analysis
        / "{Sample_Project}"
        / "BCLConvert"
        / "{Sample_ID}"
        / "fastqc"
        / "{Sample_ID}.{type}.csv",
    output:
        out_dir / "{Sample_Project}" / "fastqc" / "{Sample_ID}.{type}.csv",


### FASTQ
rule fastq_undetermined:
    input:
        in_analysis
        / "BCLConvert"
        / "fastq"
        / "Undetermined_S0_L{Lane}_R{read}_001.fastq.gz",
    output:
        out_dir / "Undetermined_S0_L{Lane}_R{read}_001.fastq.gz",
    threads: 1
    resources:
        mem="100 MiB",
    shell:
        "rsync -a --stats {input} {output} > {log} 2>&1"


use rule fastq_undetermined as fastq_sample with:
    input:
        in_analysis
        / "{Sample_Project}"
        / "BCLConvert"
        / "fastq"
        / "{Sample_Project}"
        / "{Sample_ID}_S{sample_n}_L{Lane}_R{read}_001.fastq.gz",
    output:
        out_dir
        / "{Sample_Project}"
        / "{Sample_ID}_S{sample_n}_L{Lane}_R{read}_001.fastq.gz",
