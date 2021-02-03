wget -O gencode.gtf.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.basic.annotation.gtf.gz

python get_tss_from_gencode.py gencode.gtf.gz ../test_input/regions/promoter.bed 1000
