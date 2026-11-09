ls $PWD/simulated_FASTQ/*.fq >fq.list
perl ../ViruScanNIPT.pl --config ViruScanNIPT.conf.example --fq-list fq.list --outdir ViruScanNIPT_result
