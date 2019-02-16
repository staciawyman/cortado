#!/usr/bin/perl

use Getopt::Long;

my $threads = 24;
my $window_size = 6;   
my $cortado_path = `which cortado.py`;
chomp($cortado_path);
my $fastq_path = "/data/amplicon_fastqs";

GetOptions ("t=i" => \$threads,    # numeric
	    "window=i" => \$window_size,
            "c=s"   => \$cortado_path,
            "f=s"    => \$fastq_path)    
  or die("USAGE: perl convert_manifest.pl [options]  <manifest text file>\n");


$count = 0;
#print "if [ -z \"\$1\" ]\nthen\necho \"FASTQ file directory is required as command line argument.\"\nexit 1\nfi\n";
while (<>) {
    if (/^#/) { next; }

    chomp;
    ($samp,$ref,$refseq,$donorseq,$main_site,$guideseq) = split(/\t/);
    $name = $samp."_".$ref;
    $samp .= "_";
    if ($donorseq eq "-") {
        print "/usr/bin/python $cortado_path -r1 $fastq_path/$samp\*R1_001.fastq.gz -r2 $fastq_path/$samp\*R2_001.fastq.gz -o output -n $name -a $refseq -g $guideseq --trim_sequences --trimmomatic_options_string ILLUMINACLIP:/data/applications/Trimmomatic-0.36/adapters/adapters.fa:2:30:10 --keep_intermediate  --min_identity_score 58 --window_around_sgrna $window_size  --min_frequency_alleles_around_cut_to_plot 0.1 --max_rows_alleles_around_cut_to_plot 250";
	if ($count < $threads) {
		print " & \n";
		$count++
	} else {
		print "  \n";
		$count = 0;
	}
    } else {
        print "/usr/bin/python $cortado_path -r1 $fastq_path/$samp\*R1_001.fastq.gz -r2 $fastq_path/$samp\*R2_001.fastq.gz -o output -n $name -a $refseq -e $donorseq -g $guideseq --trim_sequences --trimmomatic_options_string ILLUMINACLIP:/data/applications/Trimmomatic-0.36/adapters/adapters.fa:2:30:10 --keep_intermediate --all_edits --main_site $main_site --min_identity_score 58 --window_around_sgrna $window_size --min_frequency_alleles_around_cut_to_plot 0.1 --max_rows_alleles_around_cut_to_plot 250 ";
	if ($count < $threads) {
		print " & \n";
		$count++
	} else {
		print "  \n";
		$count = 0;
	}
    }
}
