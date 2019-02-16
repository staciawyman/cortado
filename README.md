# **Cortado**
## ABOUT
Cortado is a reimplementation of CRISPResso https://github.com/lucapinello/CRISPResso (version 1.0.8) with the following primary changes:

-Takes as input a sample manifest (format below) from which a run file of shell commands is created with a call to cortado for each sample so that unlimited numbers of samples can be run.

-A summary of all samples in a run is created in one file with one sample per line.

-HDR is reported for each edited location in the donor, with the total based on an indicated primary edit site


## REQUIREMENTS: 
Flash and Needle (part of EMBOSS) must be in your path to run cortado.

## USAGE: 
To use cortado, a sample manifest is first created with the following format:

<CENTER>

![manifest](https://github.com/staciawyman/cortado/blob/master/cortado_manifest_dirs.png)

</CENTER>

The manifest can be created in Excel by the experimenter and then copied into a text file in the directory where the editing outcomes will be analyzed. 

The manifest text file is given as input to the convert_manifest.pl script which creates a shell script of cortado commands:

`perl convert_manifest.pl manifest.txt > manifest.sh`

Then the manifest.sh file can be executed to run the cortado commands. I recommend "tmux" or "screen"  if you have a lot of samples.

`sh manifest.sh`

## DEFAULT ARGUMENTS:

-window size: 3


## OUTPUT: 
A summary file is produced with output for one sample per line. I usually copy the output_summary.txt file contents into the Excel manifest workbook to return to the experimenter. 

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

## ERROR CHECKING:
Sometimes there are errors in the input manifest. I've found the best way to check for these is to grep for "ERROR" in the log file:

`grep -i err output/*/CR*`

If you get no results from that command, then you're good. Some common errors include the guide not being found in the reference sequence or the donor sequence not being the same length as the reference sequence.

