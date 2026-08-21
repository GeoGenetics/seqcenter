### LOGS


rule logs_errors:
    output:
        "<out_dir>/{run_id}/Logs/{log}.log",
    wildcard_constraints:
        log="Errors",
    priority: 10
    threads: 1
    params:
        path=f"{in_analysis}/logs/BCLConvert",
        files=lambda w: "{MSLB2F,P2FSWT2}-stderr_*.txt",
    shell:
        "cat {params.path}/{params.files} | grep -iv warning | sort -u > {output}"


rule logs_warnings:
    output:
        "<out_dir>/{run_id}/Logs/{log}.log",
    wildcard_constraints:
        log="Warnings",
    priority: 10
    threads: 1
    params:
        path=f"{in_analysis}/logs/BCLConvert/",
        files=lambda w: "{MSLB2F,P2FSWT2}-stderr_*.txt",
    shell:
        "cat {params.path}/{params.files} | grep -i warning | sort -u > {output}"


rule logs_info:
    output:
        "<out_dir>/{run_id}/Logs/{log}.log",
    wildcard_constraints:
        log="Info",
    priority: 10
    threads: 1
    params:
        path=f"{in_analysis}/logs/BCLConvert/",
        files=lambda w: "{MSLB2F,P2FSWT2}-stdout_*.txt",
    shell:
        "cat {params.path}/{params.files} > {output}"


### REPORTS
ruleorder: reports_ss > reports_run > reports_demux > reports_bclconvert


rule reports_run:
    input:
        "<in_dir>/{report}",
    output:
        "<out_dir>/{run_id}/Reports/{report}",
    priority: 10
    threads: 1
    resources:
        mem="100 MiB",
    shell:
        "rsync --owner --group --times --links --chmod=D755,F644 {input} {output}"


use rule reports_run as reports_ss with:
    input:
        sample_sheet,
    wildcard_constraints:
        report="SampleSheet.csv",


use rule reports_run as reports_demux with:
    input:
        "<in_analysis>/Demux/{report}",


use rule reports_run as reports_bclconvert with:
    input:
        "<in_analysis>/BCLConvert/fastq/Reports/{report}",


use rule reports_run as reports_fastqc with:
    input:
        "<in_analysis>/{Sample_Project}/BCLConvert/{Sample_ID}/fastqc/{Sample_ID}.{type}.csv",
    output:
        "<out_dir>/{run_id}/{Sample_Project}/fastqc/{Sample_ID}.{type}.csv",


### Logs
use rule reports_run as logs with:
    input:
        "<in_analysis>/Logs/{log}.log",
    output:
        "<out_dir>/{run_id}/Logs/{log}.log",


### FASTQ
rule fastq_sample:
    input:
        "<in_analysis>/{Sample_Project}/BCLConvert/fastq/{Sample_Project}/{Sample_ID}_S{sample_n}_L00{Lane}_R{read}_001.fastq.gz",
    output:
        "<out_dir>/{run_id}/{Sample_Project}/{Sample_ID}_S{sample_n}_L00{Lane}_R{read}_001.fastq.gz",
    log:
        "logs/{run_id}/fastq/{Sample_Project}/{Sample_ID}_S{sample_n}_L00{Lane}_R{read}_001.fastq.gz",
    threads: 1
    resources:
        mem="100 MiB",
    shell:
        "rsync --owner --group --times --links --chmod=D755,F644 --stats {input} {output} > {log} 2>&1"


use rule fastq_sample as fastq_undetermined with:
    input:
        "<in_analysis>/BCLConvert/fastq/{Sample_ID}_S{sample_n}_L00{Lane}_R{read}_001.fastq.gz",
    output:
        "<out_dir>/{run_id}/{Sample_ID}_S{sample_n}_L00{Lane}_R{read}_001.fastq.gz",
    log:
        "logs/{run_id}/fastq/{Sample_ID}_S{sample_n}_L00{Lane}_R{read}_001.fastq.gz",
