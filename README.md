Cortado is a reimplementation of CRISPResso https://github.com/lucapinello/CRISPResso (version 1.0.8) with the following primary changes:

-Takes as input a sample manifest (format below) from which a run of file of shell commands are created with a cortado run for each sample so that unlimited numbers of samples can be run.

-A summary of all samples in a run is created in one file with the following columns: 

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

-HDR is reported for each edited location in the donor, with the total based on an indicated primary edit site


## REQUIREMENTS: Flash and Needle (part of EMBOSS) must be in your path to run cortado.

## USAGE: To use cortado, a sample manifest is first created with the following format:

<CENTER>

![manifest](https://github.com/staciawyman/cortado/blob/master/manifest.png)

</CENTER>

## OUTPUT: A summary file is produced with output for one sample per line.

