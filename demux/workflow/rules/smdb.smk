
rule smdb_upload:
    input:
        run_info=out_dir / "Reports" / "RunInfo.xml",
        demux_stats=out_dir / "Reports" / "Demultiplex_Stats.csv",
        sample_sheet=sample_sheet,
    log:
        "logs/smdb_upload.log",
    conda:
        "scripts/smdb-upload/environment.yaml"
    params:
        db_password=os.environ.get("DB_PASSWORD"),
        mail="julie.bitz-thorsen@sund.ku.dk",
        extra="--db_host dandypdb01fl --db_port 5432 --db_name smdb --schema_name uploaded_data --db_user upload_user --table_name flowcell",
    shell:
        "./scripts/smdb-upload/smdb_upload.py --path_to_demultiplex_stats {input.demux_stats} --path_to_run_info {input.run_info} --path_to_sample_sheet {input.sample_sheet} {params.extra} --db_password {params.db_password:q} --send_upload_receipts_to {params.mail:q} > {log} 2>&1"
