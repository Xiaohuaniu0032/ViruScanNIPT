ls $PWD/simulated_FASTQ/*.fq >fq.list
perl ../ViruScanNIPT.pl --config ViruScanNIPT.conf.example --fq-list fq.list --outdir ViruScanNIPT_result
perl ../bin/summarizeGroupnameResults.pl --input ViruScanNIPT_result --group-table ../database/annotation/viralGroupname.tsv --outdir ViruScanNIPT_result/Summary --positive-hit-cutoff 2

