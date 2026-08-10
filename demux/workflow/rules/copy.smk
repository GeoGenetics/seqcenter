def search_ss(wildcards, w_keys=None):
    if w_keys is None:
        w_keys = wildcards.keys()

    mask = pd.Series(True, index=ss_samples.index)
    for key in w_keys:
        mask &= ss_samples[key].eq(wildcards.get(key))
    return ss_samples.loc[mask]


def ss2sss(wildcards, w_keys=None, out_keys=None):
    x = search_ss(wildcards, w_keys)

    if not x[out_keys.keys()].isna().any().any():
        x = x[out_keys.keys()].rename(columns=out_keys)

    x = x.drop_duplicates()
    assert x.shape[0] == 1, f"invalid SSS search: {x}"
    x = x.iloc[0]

    return x.to_dict()


### REPORTS
ruleorder: reports_ss > reports_run > reports_demux > reports_bclconvert


rule reports_run:
    input:
        in_dir / "{report}",
    output:
        out_dir / "{run_id}/Reports/{report}",
    threads: 1
    resources:
        mem="100 MiB",
    shell:
        "rsync --owner --group --times --links --chmod=D555,F444 {input} {output}"


use rule reports_run as reports_ss with:
    input:
        sample_sheet,
    wildcard_constraints:
        report="SampleSheet.csv",


use rule reports_run as reports_demux with:
    input:
        in_analysis / "Demux/{report}",


use rule reports_run as reports_bclconvert with:
    input:
        in_analysis / "BCLConvert/fastq/Reports/{report}",


use rule reports_run as reports_fastqc with:
    input:
        lambda w: str(
            in_analysis
            / "BCLConvert"
            / "{Sample_ID}"
            / "fastqc"
            / "{Sample_ID}.{{type}}.csv"
        ).format(**ss2sss(w, ["Sample_ID"], out_keys={"Sample_ID_demux": "Sample_ID"})),
    output:
        out_dir / "{run_id}/{Sample_Project}/fastqc/{Sample_ID}.{type}.csv",


### Logs
use rule reports_run as logs with:
    input:
        in_analysis / "Logs/{log}.log",
    output:
        out_dir / "{run_id}/Logs/{log}.log",


### FASTQ
rule fastq_sample:
    input:
        lambda w: in_analysis
        / "BCLConvert"
        / "fastq"
        / "{Sample_ID}_S{sample_n}_L00{{Lane}}_R{{read}}_001.fastq.gz".format(
            **ss2sss(
                w,
                ["Sample_Project", "Sample_ID", "sample_n", "Lane", "read"],
                out_keys={
                    "Sample_ID_demux": "Sample_ID",
                    "sample_n_demux": "sample_n",
                },
            )
        ),
    output:
        out_dir
        / "{run_id}/{Sample_Project}/{Sample_ID}_S{sample_n}_L00{Lane}_R{read}_001.fastq.gz",
    log:
        "logs/{run_id}/fastq/{Sample_Project}/{Sample_ID}_S{sample_n}_L00{Lane}_R{read}_001.fastq.gz",
    threads: 1
    resources:
        mem="100 MiB",
    shell:
        "rsync --owner --group --times --links --chmod=D555,F444 --stats {input} {output} > {log} 2>&1"
