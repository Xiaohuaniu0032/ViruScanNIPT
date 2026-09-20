#perl extract_accession_from_gbff.pl \
#    -i /ehpcdata/fulongfei/project/AnkangCentralHospital/viral.1.genomic.gbff.gz \
#    -o viral_accession.list

perl extract_accession_taxid.pl \
    -a viral_accession.list \
    -t /ehpcdata/fulongfei/database/taxonomy/nucl_gb.accession2taxid.gz \
    -o viral_accession_taxid.tsv

