# **Cortado**
## ABOUT
Cortado is a reimplementation of CRISPResso https://github.com/lucapinello/CRISPResso (version 1.0.8) with the following primary changes:

-Takes as input a sample manifest (format below) from which a run file of shell commands is created with a call to cortado for each sample so that unlimited numbers of samples can be run.

-A summary of all samples in a run is created in one file with one sample per line.

-HDR is reported for each edited location in the donor, with the total based on an indicated primary edit site

-A read is counted as an HDR read if the exact nucleotide change indicated in the donor sequence is acheived.

Cortado has much of the same functionality as CRISPResso, but some functionality was lost in streamlining to implement our scroring scheme.  Command line options of CRISPResso can be modified in the perl script and others can be added, though they are not gauranteed to work with cortado.


## INSTALLATION REQUIREMENTS: 
Flash and Needle (part of EMBOSS) must be in your path to run cortado.

Flash: http://ccb.jhu.edu/software/FLASH/

Needle from the EMBOSS suite: ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.5.7.tar.gz

64-bit Linux executables are included in the bin directory and may work for you. Be sure to include both 
the cortado directory and bin subdirectory to your path to use these.

Cortado runs with Python version 3.8.2 and above and has not beed tested with others. 

The following Python libraries must also be installed: 

	pandas>=0.15
	numpy>=1.9
	matplotlib>=1.3.1
	biopython>=1.6.5
	argparse>=1.3

## USAGE: 
Cortado takes as input a sample manifest that has sample, reference (and donor sequence if HDR), and guide information for each sample. The sequencing fastq files should all be in one directory (and given as an argument to the convert script). These should be 300bp paired-end sequencing reads (though cortado also works on 150bp PE reads if the target is short enough for the reads to overlap).

A sample manifest is created by the experimenter (that edited the samples) with the following format:

<CENTER>

![manifest](https://github.com/staciawyman/cortado/blob/master/cortado_manifest_dirs.png)

</CENTER>

The manifest can be created in Excel by the experimenter and then copied into a text file in the directory where the editing outcomes will be analyzed. 

The manifest text file is given as input to the convert_manifest_reorient.pl script which creates a shell script of cortado commands. If you keep all your fastq files in the same directory, you can set that in the script, but otherwise you can add it to the command line.

	`perl convert_manifest_reorient.pl -f /path/to/fastq/files [options] manifest.txt > manifest.sh`

###	OPTIONS:

	-t int	Number of threads to use (runs in batches), default 24
	-w int	Window size around cut site within which to check for indels, default 3 (on each side of cut site)
	-f /path/to/fastq/dir
	-c /path/to/cortado.py
	-h 	Print help message


The number of threads is how many simultaneous versions of cortado will run at a time in batches. 
The convert script will try and infer the location of cortado.py from your path, but it can also be given as a command-line argument. 
Then the manifest.sh file is executed to run the cortado commands. I recommend running under "tmux" or "screen"  if you have a lot of samples.

	`sh manifest.sh`



## OUTPUT: 
A summary file called "output_summary.txt" is produced with output for one sample per line in the directory that cortado is called from. If you rerun manifest.sh, you should first delete this file or the new results will be concatenated onto the old. I usually copy the output_summary.txt file contents into the Excel manifest workbook to return to the experimenter. 

The output summary (example outputsummary.txt below) file has the following columns:

	-TotalReads - number of reads in the input fastq files
        -MergedReads - cortado reorients reads to all be in the same direction to simplify alignment. This is based on a subsequence of nucleotides at the 5' of the reference sequence. 
        -PercentMerged - If there happen to be a lot of errors at that at the 5' end of the read, or the reference sequence doesn't start at the primer (extra sequence in the reference), then the percent merged might be low and you should review your reference sequence to make sure it is correct.
	-AlignedReads - this is used as the denominator for %NHEJ and %HDR.
	-PercentAligned	- this should be very high (>90) and if it is not, there may be a problem with your sequencing data or reference sequence.
	-Unmodified _ number of unmodified reads
	-%Unmodified - number of unmodified reads divided by the number of aligned reads
	-CutsiteSubs - SNP within the cut site window (6bp by default) are not counted as indels and thus are just recorded here but do not affect editing outcome percents. If this number is very large, then it may be that your sample has a SNP at the cut site or there were PCR errors at the cut site.
	-%CutsiteSubs	
	-NHEJ - number of reads with insertion or deletion occurring withing 3bp to either side of cutsite (with default windowsize of 6bp, windowsize can be changed by user)	
	-%NHEJ - number of NHEJ reads/aligned reads
	-MainSite - number of reads edited at the primary site (given by number in the analysis manifest)	
	-Main% - number of MainSite reads/aligned reads	
	-EditSiteN - number of reads with site at N edited	
	-Edit% - number of EditSiteN reads/aligned reads
	-All_HDR_Pos - number of reads where all possible sites are edited	
	-All_Edit% - number of All_HDR_Pos reads/aligned reads

Running manifest.sh creates a directory called "output" where the cortado output files for each run are placed for each sample in its own directory. Here you can find many helpful files if a sample failed to run. You can also find a pdf image of the aligned alleles with indels and SNPs marked. This is called "9.Alleles_around_cut_site_for_<sample_name>.pdf." 

![sample_output](https://github.com/staciawyman/cortado/blob/master/sample_output.png)

## BASE EDITING FUNCTIONALITY:
Cortado can be used to assess base editing, though it's a bit of a hack. If you put all the expected sites of base editing into the donor sequence with the expected nucleotide change (as though you are doing HDR), then the number of successful edits at all of those locations will be output (counted independently). You should also put a "1" in the MainSite column.

## ERROR CHECKING:
Errors are output interleaved with the rest of the output.  Some common errors include not getting enough sequence reads to analyze the sample (<1K reads). Errors in the manifest will be output when trying to convert it. Those usually are that the guide is not found in the reference sequence or the donor sequence is not the same length as the reference sequence. The first place to look to start debugging is the log file in the output/<sample_name>/cortado_running_log.txt directory for the failed sample.

