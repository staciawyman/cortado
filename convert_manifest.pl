#!/usr/bin/perl

$threads = 24;
$count = 0;
$window_size = 6;   
#$window_size = 20; 
#print "if [ -z \"\$1\" ]\nthen\necho \"FASTQ file directory is required as command line argument.\"\nexit 1\nfi\n";
while (<>) {
    if (/^#/) { next; }

    chomp;
    ($samp,$ref,$refseq,$donorseq,$main_site,$guideseq) = split(/\t/);
    if ($samp =~ /(JVNGS\d\d\d)_/) { $run = $1; }
    $name = $samp."_".$ref;
    $samp .= "_";
    if ($donorseq eq "-") {
        print "/usr/bin/python /data/applications/cortado/cortado.py -r1 /data/amplicon_fastqs/$run/$samp\*R1_001.fastq.gz -r2 /data/amplicon_fastqs/$run/$samp\*R2_001.fastq.gz -o output -n $name -a $refseq -g $guideseq --trim_sequences --trimmomatic_options_string ILLUMINACLIP:/data/applications/Trimmomatic-0.36/adapters/adapters.fa:2:30:10 --keep_intermediate  --min_identity_score 58 --window_around_sgrna $window_size  --min_frequency_alleles_around_cut_to_plot 0.1 --max_rows_alleles_around_cut_to_plot 250";
	if ($count < $threads) {
		print " & \n";
		$count++
	} else {
		print "  \n";
		$count = 0;
	}
    } else {
        print "/usr/bin/python /data/applications/cortado/cortado.py -r1 /data/amplicon_fastqs/$run/$samp\*R1_001.fastq.gz -r2 /data/amplicon_fastqs/$run/$samp\*R2_001.fastq.gz -o output -n $name -a $refseq -e $donorseq -g $guideseq --trim_sequences --trimmomatic_options_string ILLUMINACLIP:/data/applications/Trimmomatic-0.36/adapters/adapters.fa:2:30:10 --keep_intermediate --all_edits --main_site $main_site --min_identity_score 58 --window_around_sgrna $window_size --min_frequency_alleles_around_cut_to_plot 0.1 --max_rows_alleles_around_cut_to_plot 250 ";
	if ($count < $threads) {
		print " & \n";
		$count++
	} else {
		print "  \n";
		$count = 0;
	}
    }
}
