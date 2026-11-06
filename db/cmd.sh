#perl extract_gbff_metadata.pl /ehpcdata/fulongfei/project/AnkangCentralHospital/viral.1.genomic.gbff.gz >viral.genomic.metadata.tsv

#cd /ehpcdata/fulongfei/project/AnkangNIPT_20260828_reanalysis/parse_viral_gbff

#perl add_genus_to_metadata.pl \
#    viral.genomic.metadata.tsv \
#    /ehpcdata/fulongfei/database/taxonomy/nodes.dmp \
#    /ehpcdata/fulongfei/database/taxonomy/names.dmp \
#    viral.genomic.metadata.with_genus.tsv

#perl filter_viral_genomic_dna_host.pl \
#viral.genomic.metadata.with_genus.tsv \
#viral.genomic.metadata.with_genus.genomicDNA.tsv

#(base) [fulongfei@manager parse_viral_gbff]$ perl filter_phage_organism.pl

#perl annotate_virushostdb_host.pl \
#    viral.genomic.metadata.with_genus.genomicDNA.tsv \
#    /ehpcdata/fulongfei/database/VirusHostDB/virushostdb.virus_host.tsv \
#    viral.genomic.metadata.with_genus.genomicDNA.host_annotated.tsv

perl filter_virushostdb_human_host.pl \
    viral.genomic.metadata.with_genus.genomicDNA.host_annotated.tsv \
    viral.genomic.metadata.with_genus.genomicDNA.host_human.tsv
