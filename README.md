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

Cortado runs with Python version 2.7.5 and has not beed tested with others. To temporarily set your version to 2.7, create a virtual environment to run cortado:

	$ conda create -n py27 python=2.7 anaconda

This will prompt you to install many packages for this version of python. Then to start the virtual environment, making python2.7 your default version:

	$ source activate py27

To deactivate:

	$ source deactivate

(logging out will also deactivate your virtual environment)

The following Python libraries must also be installed: pandas, numpy, matplotlib, pylab, Bio, seaborn

## USAGE: 
Cortado takes as input a sample manifest that has sample, reference and donor sequence, and guide information for each sample. The sequencing fastq files should all be in one directory (and given as an argument to the convert script). These should be 300bp paired-end sequencing reads.

A sample manifest is created by the experimenter (that edited the samples) with the following format:

<CENTER>

![manifest](https://github.com/staciawyman/cortado/blob/master/cortado_manifest_dirs.png)

</CENTER>

The manifest can be created in Excel by the experimenter and then copied into a text file in the directory where the editing outcomes will be analyzed. 

The manifest text file is given as input to the convert_manifest.pl script which creates a shell script of cortado commands. If you keep all your fastq files in the same directory, you can set that in the script, but otherwise you can add it to the command line.

	`perl convert_manifest.pl [options] manifest.txt > manifest.sh`

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

The output file has the following columns:

	-Total reads
	-Aligned reads
	-PercentAligned	
	-Unmodified reads	
	-%Unmodified	
	-CutsiteSubs	
	-%CutsiteSubs	
	-NHEJ-number of reads with insertion or deletion occurring withing 3bp to either side of cutsite (windowsize can be changed by user)	
	-%NHEJ (NHEJ reads/aligned reads)
	-MainSite-reads edited at the primary site	
	-Main% (MainSite reads/aligned reads)	
	-EditSiteN-number of reads with this site edited	
	-Edit% (EditSiteN reads/aligned reads)	
	-All_HDR_Pos-number of reads where all possible sites are edited	
	-All_Edit% (All_HDR_Pos/aligned reads)

Running manifest.sh also creates a directory called "output" where the cortado output files for each run are placed for each sample in its own directory. Here you can find many helpful file if a sample failed to run. You can also find a pdf image of the aligned alleles with indels and SNPs marked. This is called "9.Alleles_around_cut_site_for_<sample_name>.pdf." Note that there is a bug in this representation that shows the allele from the forward read alignment and reverse read alignment (of the same allele) separately. The abundance should be summed for these two.

![sample_output](https://github.com/staciawyman/cortado/blob/master/sample_output.png)

## ERROR CHECKING:
Sometimes there are errors in the input manifest and it causes of one of the instances of cortado to fail. I've found the best way to check for these is to grep for "ERROR" in the log file:

`grep -i err output/*/cor*`

If you get no results from that command, then you're good. Some common errors include the guide not being found in the reference sequence or the donor sequence not being the same length as the reference sequence. The first place to look to start debugging is the log file in the output/<sample_name>/cortado_running_log.txt directory for the failed sample.

