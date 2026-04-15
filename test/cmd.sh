ls $PWD/data/*.fq.gz >fq.list

perl ../ViruScanNIPT.pl --config ../conf/ViruScanNIPT.conf --fq-list fq.list --outdir ViruScanNIPT_result

