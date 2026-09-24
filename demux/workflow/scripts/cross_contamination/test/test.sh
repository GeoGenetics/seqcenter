

TYPE=20251013_A00706_0986_AHHTKKDSXF
# With index list
python3 ../cross_contamination.py --loglevel ERROR --index-counts $TYPE.Index_Hopping_Counts.csv --index-known ../../../../resources/eDNA_index_list_20250807.tsv --lanes 1 --rpm-warn 100 --out-prefix /tmp/cc1_idx.$TYPE
python3 ../cross_contamination.py --loglevel ERROR --index-counts $TYPE.Index_Hopping_Counts.csv --index-known ../../../../resources/eDNA_index_list_20250807.tsv --lanes 2 --rpm-warn 100 --out-prefix /tmp/cc2_idx.$TYPE
python3 ../cross_contamination.py --loglevel ERROR --index-counts $TYPE.Index_Hopping_Counts.csv --index-known ../../../../resources/eDNA_index_list_20250807.tsv --lanes 3,4 --rpm-warn 100 --out-prefix /tmp/cc34_idx.$TYPE

# Without index list
python3 ../cross_contamination.py --loglevel ERROR --index-counts $TYPE.Index_Hopping_Counts.csv --lanes 1 --rpm-warn 100 --out-prefix /tmp/cc1_noidx.$TYPE
python3 ../cross_contamination.py --loglevel ERROR --index-counts $TYPE.Index_Hopping_Counts.csv --lanes 2 --rpm-warn 100 --out-prefix /tmp/cc2_noidx.$TYPE
python3 ../cross_contamination.py --loglevel ERROR --index-counts $TYPE.Index_Hopping_Counts.csv --lanes 3,4 --rpm-warn 100 --out-prefix /tmp/cc34_noidx.$TYPE

TYPE=custom_indexes
python3 ../cross_contamination.py --loglevel ERROR --index-counts $TYPE.Index_Hopping_Counts.csv --rpm-warn 100 --out-prefix /tmp/cc_idx.$TYPE

# Check output
md5sum -c test.md5
rm /tmp/cc_*
