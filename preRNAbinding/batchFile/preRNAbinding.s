#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time 2:00:00              # Runtime in D-HH:MM:SS
#SBATCH --mem=40GB               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name=preRNAbinding
#SBATCH --mail-user=[youremail] # This needs to be changed to the email to which notifications will be sent
#SBATCH --mail-type=END,FAIL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH -e /home/[yourusername]/preRNAbinding/batchFile/preRNAbinding.err # This needs to be changed to the path to the subdirectory "batchFile" of the directory that the TAR is untarred into.
#SBATCH -o /home/[yourusername]/preRNAbinding/batchFile/preRNAbinding.out # This needs to be changed to the path to the subdirectory "batchFile" of the directory that the TAR is untarred into.


##########################
# SECTION 1: DOCUMENTATION
##########################

# NAME: preRNAbinding.s - Plot distance histograms of protein binding sites to all types of pre-RNA junctions (intron-exon and exon-intron, etc.), bedtools closest requires post-processing to isolate pertinent information.  Instructions for this one-command batch script working example are provided below.  "preRNAbinding.s" also uses ngs.plot.r to create exon and genebody distance graphs.

# Background: The plots of these calculations are included and explained in my masters thesis "Lipscomb_Thomas_Thesis.pdf" included in directory "batchFile".  To understand the structure of pre-RNA see "Gene_structure_eukaryote_2_annotated.png" in "batchFile".

# The following code is an example, which generates graphs of QKI binding frequency near pre-RNA junctions, used in the attached thesis "Lipscomb_Thomas_Thesis.pdf".

# DESCRIPTION:
# Stdin: Files in batchFile directory.
# Stdout: The same as the files in the afterRunningBatchFile example output directory.
# Named input file: None.
# Named output file: None.

# USED BY: command line: $ sbatch preRNAbinding.s
# USES:
# • lastexons.awk
# • onlylastexons.py
# • removelastexons.py

# SYNOPSIS:
# Download from this directory: https://drive.google.com/open?id=1PnudvWR_xbyfqaJJve5R44PHeCyQBBRr file preRNAbinding.tar by right-clicking the TAR file and selecting "Download".  Untarring that file creates a subdirectory "preRNAbinding", which contains the contents of the TAR file.
# To untar:
# $ tar -xvf preRNAbinding.tar
# $ cd preRNAbinding
# $ ls
# Expect to see a working directory "batchFile" and an example of a working directory after a completed run "afterRunningBatchFile".
# $ cd batchFile
# $ ls
# $ vi README.txt
# This readme contains the contents of this forum post plus more background information.

# "batchFile" contains the batch script "preRNAbinding.s" and the input files needed to generate output data files, after which the contents of directory "batchFile" should be the same as the contents of "afterRunningBatchFile".

# PREREQUISITES ON PERSONAL COMPUTER:
## (1) Install a VPN client to log into your institution's VPN.  Example instructions: https://www.nyu.edu/life/information-technology/getting-started/network-and-connectivity/vpn.html
## (2) Install a remote terminal and ssh client to run commands and download the results.  Instructions: https://mobaxterm.mobatek.net/
# OPTIONAL (EITHER INSTALL ON PC OR ENTERPRISE SERVER):
## (3) Install R 3.0.1+ from http://cran.rstudio.com/ as a corequisite of RStudio Desktop (immediately below) on either the enterprise server or your personal computer.
## (4) Install RStudio Desktop from https://www.rstudio.com/products/rstudio/download/#download on either the enterprise server or your personal computer.
## (5) If you installed RStudio Desktop on the Enterprise Server and will run RStudio in interactive mode on the server terminal from a Windows machine, you will need Xming and Putty or their equivalents installed on your Windows machine. Instructions: https://wikis.nyu.edu/display/NYUHPC/Running+RStudio+on+the+HPC+Prince+Cluster

# PREREQUISITES ON ENTERPRISE SERVER:
## (1) Install the "Environment Modules" management tool in order to do "module avail" and "module load" (how to install: http://modules.sourceforge.net/) (commands: https://wikis.nyu.edu/display/NYUHPC/Slurm+Tutorial#SlurmTutorial-SoftwareandEnvironmentModules).
## (2) Use "module avail" to list all available modules and check if the immediately below modules are available.  Also check if the version numbers are compatible with each other (e.g. Module ngsplot/20170524 was setup based on module r/intel/3.3.2 not the more recent version of r/intel).
# $ module avail
# LIST OF REQUIRED MODULES INCLUDING THE VERSION NUMBERS I USED:
# rstudio/1.1.383
# r/intel/3.3.2
# Module ngsplot/20170524 was setup based on module r/intel/3.3.2, which is why the following was done: "module swap r/intel  r/intel/3.3.2".
# ngsplot/20170524
# bedtools/intel/2.26.0
# samtools/intel/1.6
# python/intel/2.7.12
## (3) At the top of "preRNAbinding.s" change mail-user=<youremail> to your email address.
## (4) At the top of "preRNAbinding.s" change "SBATCH -e ..." to the path to the .out and .err files that should be generated in the subdirectory "batchFile" of the directory that the TAR was untarred into.
## (5) At the top of "preRNAbinding.s" change "SBATCH -o ..." to the path to the .out and .err files that should be generated in the subdirectory "batchFile" of the directory that the TAR was untarred into.
## (6) Change the explicit path to singularity in the two commands in the following BEGIN NGSPLOT section to the explicit path on your machine.  Test to see if singularity is installed:
# $ ls /share/apps/singularity/3.0.1/bin/singularity
# One should see the singularity directory which ngs.plot.r uses.
## (7) Change the explicit path to the ngsplot .img file in the two commands in the following BEGIN NGSPLOT section to the explicit path on your machine.  Test to see if ngsplot .img is installed:
# $ ls /beegfs/work/public/singularity/ngsplot*.img
# One should see the .img which ngs.plot.r uses (the .img file used in this code was ngsplot-2.61.img).
# Immediately below is how to arrange the two explicit paths in the ngs.plot.r.  To see what the [arguments] should be, see when this is used in the BEGIN NGSPLOT section:
# /share/apps/singularity/3.0.1/bin/singularity exec /beegfs/work/public/singularity/ngsplot-2.61.img ngs.plot.r [arguments]
# If the singularity and .img explicit paths are not changed and the version of the .img is not checked, then ngs.plot.r may not work and in that case one will not get the exon and genebody graphs from ngs.plot.r, but one will still get the post-processing of bedtools closest that goes into the RStudio distance histograms of RNA binding protein binding sites to junctions of all types that a pre-RNA can have.

# Before running the below sbatch command, check to see that the number of starting files is 26, by executing:
# $ ls | wc
# The first number is how many files are in the directory.
# If the number of files is not 26, then the directory may not have all of the required starting files in batchFile needed to run the below sbatch successfully.

# ONE-COMMAND EXECUTION:
# $ sbatch preRNAbinding.s
# sbatch should take less than an hour on an enterprise server.

# The sbatch command above automatically runs these scripts:
# lastexons.awk
# removelastexons.py
# onlylastexons.py

# To monitor sbatch progress, type squeue -u [username], e.g.: squeue -u thl312.

# Wait for sbatch to finish.  When sbatch is finished, check to see that the number of ending files is 290, by executing:
# $ ls | wc
# The first number is how many files are in the directory.
# If the number of files is not 290, then the above sbatch may have only partially run and ran into an error.  Check "preRNAbinding.err" and "preRNAbinding.out" for verbose information about what might have gone wrong.  Also run the below to find the last file produced by the sbatch:
# $ ls -ltr
# The most recent file should be "hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_regionstart_lastexon.bed".  If that is not the most recent file, the sbatch probably terminated early.

# The final step is to run RStudio Knit to produce the PNG histograms of how far away QKI binds from RNA junctions.  There are three possibilities to use RStudio Knit to generate "Hg38_QKI_closest_dist_freqV9.html":
# • On your personal computer the "Hg38_QKI_closest_dist_freqV9.Rmd" file may be opened in the RStudio Desktop and click "Knit" to generate the HTML file that contains the PNG histograms that are generated from the data files.
# • Run RStudio as a module on the Enterprise Server in interactive mode:
# $ srun --x11 --pty /bin/bash
# $ module load rstudio/1.1.383
# $ rstudio
# Click "Knit" to generate the HTML file that contains the PNG histograms that are generated from the data files.
# • Run RStudio as a module on the Enterprise Server in command line mode (however this produces smaller PNG histograms so their titles are cut off):
# $ module load rstudio/1.1.383
# $ Rscript -e "knitr::stitch_rmd('Hg38_QKI_closest_dist_freqV9.Rmd')"

# Open "Hg38_QKI_closest_dist_freqV9.html" to view the PNG histograms produced above.
# Open "PARCLIP_QKI_Hafner2010c_hg19_noWeirdChr_sorted.exon.avgprof.pdf" to view the exon ngs.plot.r.
# Open "PARCLIP_QKI_Hafner2010c_hg19_noWeirdChr_sorted.genebody.avgprof.pdf" to view the genebody ngs.plot.r.

# Compare the RStudio histograms in "Hg38_QKI_closest_dist_freqV9.html" to the histograms in my masters thesis and compare the ngs.plot.r exon and genebody graphs to the exon and genebody graphs in my masters thesis to verify that the run is correct and complete:

# The junction histograms in my Masters thesis are:
# • Figure 14 (left): intron-exon (intron and exon are measured by Solution to Limitation 1 and Solution to Limitation 2 (described below)).
# • Figure 14 (right): exon-intron (intron and exon are measured by Solution to Limitation 1 and Solution to Limitation 2).
# • Figure 16 (left): intron-start of cassette exon junction (intron is measured by bedtools closest -io -id and cassette exon is measured by Solution to Limitation 1 and Solution to Limitation 2).
# • Figure 16 (right): end of cassette exon-intron junction (intron is measured by bedtools closest -io -iu and cassette exon is measured by Solution to Limitation 1 and Solution to Limitation 2).
# • Figure 17 (left): end of 5' UTR-CDS junction (CDS is measured by bedtools closest -io -iu and 5' UTR is measured by Solution to Limitation 1 and Solution to Limitation 2).
# • Figure 17 (right): CDS-start of 3' UTR junction (CDS is measured by bedtools closest -io -id and 5' UTR is measured by Solution to Limitation 1 and Solution to Limitation 2).
# The "Should be the same as Figure 17" histogram in the RMD file double-checks the methodology of "Solution to Limitation 2" (below), hoping the results will be the same as "Figure 17"'s use of bedtools closest.  The two histograms do not agree, possibly because "Should be the same as Figure 17" may confuse whether exons are part of the UTR (see "Gene_structure_eukaryote_2_annotated.png" the exon is part of the UTR).  For more information see the following DOUBLE-CHECK section.
# • Figure 18 (left): start of 5' UTR junction (Transcription Start Site) (Upstream of the start of the 5' UTR was not measured because that is beyond the pre-RNA and 5' UTR is measured by Solution to Limitation 1 and Solution to Limitation 2).
# • Figure 18 (right): end of 3' UTR junction (Transcription End Site) (Downstream of the end of the 3' UTR was not measured because that is beyond the pre-RNA and 3' UTR is measured by Solution to Limitation 1 and Solution to Limitation 2).
# • Figure 20 (left): intron-start of circRNA junction (intron is measured by bedtools closest -io -id and circRNA is measured by Solution to Limitation 1 and Solution to Limitation 2).
# • Figure 20 (right): end of a circRNA-intron junction (intron is measured by bedtools closest -io -iu and circRNA is measured by Solution to Limitation 1 and Solution to Limitation 2).

# The ngs.plot.r plots in my Masters thesis are:
# • Figure 15: Idealized Exon and regions 500nt upstream and 500nt downstream.
# • Figure 19: Idealized Gene and regions 2000nt upstream and 2000nt downstream.
# The ngsplot heatmap does not work because there are too many exons or genes to visualize on one page (the red lines are too thin).

# BEDTOOLS LIMITATIONS AND HOW SBATCH FILE "preRNAbinding.s" OVERCOMES THEM:

# • BEDTOOLS CLOSEST LIMITATIONS:

# Limitation 1: Bedtools closest returns a mixture of needed and unneeded information when trying to find the distance from an RNA binding protein site to an a junction.  Below is one example of Limitation 1 regarding creating data for the intron-exon junction histogram and the exon-intron junction histogram.  If the nearest RNA binding protein binding site is beyond the flanking region then that protein binding site is erroneously mixed in with the flanking region results (see Figure 7A of the thesis PDF "Lipscomb_Thomas_Thesis.pdf" in the batchFile folder).  Also, bedtools closest returns a zero for distance if the RNA binding protein site and preRNA region overlap, instead of the needed distance to either the start or end of region, whichever is closer:

# Limitation 1A: "bedtools closest -io -id -a <RNA binding protein binding sites BED file> -b <exon locations BED file>" finds the distance to the nearest upstream RNA binding protein binding site from the intron-exon junciton, but does not limit the results to only RNA binding protein binding sites within the flanking upstream intron, which is what is needed.  If the nearest RNA binding protein binding site is beyond the flanking upstream intron then that is mixed in with the flanking upstream intron results (see Figure 7A of the thesis PDF "Lipscomb_Thomas_Thesis.pdf" in in the preRNAbinding folder).

# Limitation 1B: "bedtools closest -a <RNA binding protein binding sites BED file> -b <exon locations BED file>" returns unwanted upstream and downstream values, which should be removed leaving the zeroes, which signify overlap (the RNA binding protein site is inside the junction).  Also, the zeroes are not the distance from the RNA binding protein binding site to the start or end of the exon (whichever is closer), which is what is needed.

# Limitation 1C: "bedtools closest -iu -io -a <RNA binding protein binding sites BED file> -b <exon locations BED file>" finds the distance to the nearest downstream RNA binding protein binding site from the exon-intron junciton, but does not limit the results to only RNA binding protein binding sites within the flanking downstream intron, which is what is needed.  If the nearest RNA binding protein binding site is beyond the flanking downstream intron then that is mixed in with the flanking downstream intron results (see Figure 7A of the thesis PDF in in the preRNAbinding folder).

# Limitation 2: Erroneous double-counting of QKI binding sites inside the region can happen if the region is short (e.g. exons can be short), causing the two junction histograms to overlap (see Figure 8A of the thesis PDF in in the preRNAbinding folder).

# • HOW BEDTOOLS CLOSEST LIMITATIONS ARE OVERCOME by "preRNAbinding.s":

# About Solutions: Post-processing, using AWK and python, separates the needed from the unneeded information returned by bedtools closest.  RStudio generates distance histograms of RNA binding protein binding sites to junctions of all types by importing files from post-processing.  ngs.plot.r generates exon and genebody RNA binding protein binding site frequency graphs by importing the BAM file from post-processing.


# Solution to Limitation 1: Where bedtools closest returns not zero, throw those results out.  One cannot run "bedtools closest -iu -id" because it returns an error, so run "bedtools closest" without "-iu" and "-id".  When bedtools closest returns zero, signifying that the RNA binding protein bound within the region (either an exon or intron), AWK code provided measures the distance to the start or the end of the region, whichever is closer (see Figure 7B of the thesis PDF in in the preRNAbinding folder).  Looking within adjacent regions prevents erroneous inclusion of more distant protein binding sites.

# Solution to Limitation 2: Split the center region in half and split the flanking regions in half.  This prevents erroneous double-counting, because no matter what the nucleotide distance on the x-axis of the junction histogram is, the junction histogram is given only the distances that were closer to either the start or the end of the center region, whichever the histogram includes (see Figure 8B of the thesis PDF in the batchFile folder).  This solution fully applied to intron-exon and exon-intron histograms.  This solution only partially applied to the other histograms, which split only the center region in half because the size of the flanking regions is unknown.

# The scope of the solution to Limitation 1 and Limitation 2 is complete for intron-exon and exon-intron junction histograms.  Both introns and exons can be split in half, because the size of the flanking introns are known because they are in the introns BED file and the size of the exon in the middle is known because it is in the exons BED file.

# However, the scope of the solution to Limitation 1 and Limitation 2 is more limited for junction histograms other than intron-exon and exon-intron (e.g. upstream-cassette exon and cassette exon-downstream), because the upstream and downstream cannot be split in half because the size of the flanking upstream region and flanking downstream region are unknown (In this case those regions are probably the flanking introns, which are a subset of all introns in the intron BED file).  Instead of attempting Future Research 2 (in "README.txt" file in directory "batchFile") to try to limit the results to only flanking introns, one good choice is -300nt upstream and +300nt downstream as the x-axis of my junction histograms, because that is too small to go beyond the flanking introns (Limitation 2) except possibly in rare cases.


# POTENTIAL PROBLEM WITH ALTERNATE PATH TO CALCULATE FIGURE 17 WITHOUT BEDTOOLS CLOSEST UPSTREAM AND DOWNSTREAM:

# "Figure 17" should be correct anyway because it used bedtools closest to find the upstream and downstream regions and that should work in any case, but "Should be the same as Figure 17" (see the RMD file) is not correct because it does not produce the same histograms as "Figure 17".
# • The UTR regions do agree (as expected because both paths did the same thing: bedtools closest to find QKI binding sites overlapping with the UTR, then split the UTR in half).
# • However, the first and last exon regions do not agree (bedtools closest upstream and downstream was used for "Figure 17" and trimmed BED files containing only the first exon or only the last exon were used for "Should be the same as Figure 17").
# • A possible explanation why is that the first exon in the exon BED file may contain both the 5' UTR and the first exon after the UTR and the last exon in the exon BED file may contain both the last exon before the UTR and the 3' UTR.  See the pictured immature RNA in "Gene_structure_eukaryote_2_annotated.png", the first exon includes the 5' UTR and the last exon includes the 3' UTR.  If that is the case, to produce a correct "Should be the same as Figure 17" histogram, the exon BED file would need to trimmed so that the 5' UTR is removed from the first exon and 3' UTR is removed from the last exon.

# "Should be the same as Figure 17" is generated by post-processing explained in "Solution to Limitation 1" and "Solution to Limitation 2" instead of bedtools closest upstream (-io -id) and downstream (-iu -io).
# • To generate the file that contains only first exons: "awk '$10 ~ /exon_0/' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_firstexon.bed".  This was combined with the 5' UTR file to generate "Should be the same as Figure 17 (left)".
# • To generate the file with the last exons, lastexons.awk, removelastexons.py, and onlylastexons.py were used to generate: "hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_lastexon_sorted_dist.bed".  This was combined with the 3' UTR file to generate "Should be the same as Figure 17 (right)".


# INSTRUCTIONS ON HOW TO RUN:

# preRNAbinding.s is compute intensive so it needs to be run in sbatch or srun.  Run in sbatch not srun, because srun executes one line at a time on the command line (one must paste the lines into the command line one at a time), while sbatch refers to a file that can run any number of lines, so sbatch is much more convenient.
# ssh -XC thl312@prince.hpc.nyu.edu
# srun --cpus-per-task=20 --mem=62GB --time=72:00:00 --pty /bin/bash

#################################
# END OF SECTION 1: DOCUMENTATION
#################################


############################################
# SECTION 2: MANUALLY DOWNLOADING DATA FILES
############################################

# THE WAY I DOWNLOADED THE INPUT FILES FOR THE BATCHFILE DIRECTORY:

# Download genes.
# Go to https://genome.ucsc.edu/cgi-bin/hgTables
# Select assembly: hg38
# Select table:knownGene
# Select output format: all fields from selected table
# Select output file: hg38_genes.bed
# Press get output
# Create one BED record per: Whole Gene
# Press get BED
# NOTE: hg38_genes.bed starts with the Transcription Start Site (col 2) and ends with the Transcription End Site (col 3).
# Download hg38_genes.bed.

# Download exons (includes both constitutive exons and skipped exons).
# Go to https://genome.ucsc.edu/cgi-bin/hgTables
# Select assembly: hg38
# Select track: GENCODE V24
# Select table: knownGene
# Select output format: BED - browser extensible data
# Select output file: hg38_exons.bed
# Press get output
# Create one BED record per: Exons
# Press get BED
# Download hg38_exons.bed.

# Download introns.
# Go to https://genome.ucsc.edu/cgi-bin/hgTables
# Select assembly: hg38
# Select track: GENCODE V24
# Select table: knownGene
# Select output format: BED - browser extensible data
# Select output file: hg38_introns.bed
# Press get output
# Create one BED record per: Introns
# Press get BED
# Download hg38_introns.bed.

# Download constitutive exons.
# Go to https://genome.ucsc.edu/cgi-bin/hgTables
# Select assembly: hg38
# Select track: GENCODE V24
# Select table: knownGene
# Select output format: BED - browser extensible data
# Select output file: hg38_CDS.bed
# Press get output
# Create one BED record per: Coding Exons
# Press get BED
# Download hg38_CDS.bed.

# Download 5' UTRs.
# Go to https://genome.ucsc.edu/cgi-bin/hgTables
# Select assembly: hg38
# Select track: GENCODE V24
# Select table: knownGene
# Select output format: BED - browser extensible data
# Select output file: hg38_5p.bed
# Press get output
# Create one BED record per: 5' UTR Exons
# Press get BED
# Download hg38_5p.bed.

# Download 3' UTRs.
# Go to https://genome.ucsc.edu/cgi-bin/hgTables
# Select assembly: hg38
# Select track: GENCODE V24
# Select table: knownGene
# Select output format: BED - browser extensible data
# Select output file: hg38_3p.bed
# Press get output
# Create one BED record per: 3' UTR Exons
# Press get BED
# Download hg38_3p.bed.

# Download cassette exons.
# Go to http://hexevent.mmg.uci.edu/cgi-bin/HEXEvent/HEXEventWEB.cgi
# It is hg38
# Select Whole chromosome: all
# Select Type of exon: cassette
# Select Would you like a strict definition of exon types?: No
# Select Result display: file download
# Click Submit
# Download result.txt.
# Rename result.txt to cassette_exons_notstrict.txt.
mv result.txt cassette_exons_notstrict.txt
# TO DO: do this with a rename instead of mv.
# https://www.maketecheasier.com/rename-files-in-linux/

# Rearrange cassette_exons_notstrict.txt to a bed file by:
# Remove header from cassette_exons_notstrict.txt and call it cassette_exons_notstrict_nohead.txt.
tail -n +2 cassette_exons_notstrict.txt >> cassette_exons_notstrict_nohead.txt
# sed '1d' $file >> headerless.txt is the other way to do it besides tail
# https://askubuntu.com/questions/862804/which-is-faster-to-delete-first-line-in-file-sed-or-tail
# Convert cassette_exons_notstrict_nohead.txt to BED format.
awk -F $'\t' '{print $1, $3, $4, "exon"NR, 0, $2}' cassette_exons_notstrict_nohead.txt > cassette_exons_notstrict_nohead.bed
# Replace all space with tab, because that is proper BED format.
sed -e "s/\ /\t/g" cassette_exons_notstrict_nohead.bed > temp$$
cat temp$$ > cassette_exons_notstrict_nohead.bed
rm -f temp$$
# Remove 			exon74082	0	 from the final line of the bed file.
# https://unix.stackexchange.com/questions/52779/delete-last-line-from-the-file
sed '$d' cassette_exons_notstrict_nohead.bed > temp$$
cat temp$$ > cassette_exons_notstrict_nohead.bed
rm -f temp$$

# Download circular RNA.
# Go to http://www.circbase.org/cgi-bin/downloads.cgi
# Download http://www.circbase.org/download/hsa_hg19_circRNA.bed
# hsa_hg19_circRNA.bed

# Converts hsa_hg19_circRNA.bed to hsa_hg38_circRNA.bed.
# http://genome.ucsc.edu/cgi-bin/hgLiftOver
# Original Genome: Human
# Original Assembly: hg19
# New Genome: Human
# New Assembly: hg38
# Drag BED file over Choose File box near the bottom of the webpage.
# Click Submit File
# Wait for Upload to finish.
# Click View Conversions to download converted BED file.
# Rename hglft_genome_b55_e32980.bed to hsa_hg38_circRNA.bed.
mv hglft_genome_*.bed hsa_hg38_circRNA.bed
# TO DO: do this with a rename instead of mv.
# https://www.maketecheasier.com/rename-files-in-linux/

# Go to http://dorina.mdc-berlin.de/regulators
# Select hg19
# Download http://dorina.mdc-berlin.de/api/v1.0/download/regulator/hg19/PARCLIP_QKI_Hafner2010c_hg19
# This is the QKI binding site BED file that came from this paper: "Transcriptome-wide identification of RNA-binding protein and microRNA target sites by PAR-CLIP."

# Converts PARCLIP_QKI_Hafner2010c_hg19.bed to PARCLIP_QKI_Hafner2010c_hg38.bed
# http://genome.ucsc.edu/cgi-bin/hgLiftOver
# Original Genome: Human
# Original Assembly: hg19
# New Genome: Human
# New Assembly: hg38
# Drag BED file over Choose File box near the bottom of the webpage.
# Click Submit File
# Wait for Upload to finish.
# Click View Conversions to download converted BED file.
# Rename hglft_genome_19a6_e340a0.bed to PARCLIP_QKI_Hafner2010c_hg38.bed.
mv hglft_genome_*.bed PARCLIP_QKI_Hafner2010c_hg38.bed
# TO DO: do this with a rename instead of mv.
# https://www.maketecheasier.com/rename-files-in-linux/

# Download human.hg19.genome for use in bedToBam to convert PARCLIP_QKI_Hafner2010c_hg19.bed to bam file to use in ngs.plot.r.
# https://github.com/arq5x/bedtools/tree/master/genomes
# Go to https://raw.githubusercontent.com/arq5x/bedtools/master/genomes/human.hg19.genome then do Control+S or Right Click Save As to download the file.
mv human.hg19.genome.txt human.hg19.genome
# TO DO: do this with a rename instead of mv.
# https://www.maketecheasier.com/rename-files-in-linux/

###################################################
# END OF SECTION 2: MANUALLY DOWNLOADING DATA FILES
###################################################


#######################################################
# SECTION 3: ONE-COMMAND SBATCH SCRIPT EXECUTABLE LINES
#######################################################

# Running command sbatch preRNAbinding.s above causes the following to be executed:

# In order to do "module avail" and "module load" one must have the "Environment Modules" management tool (example: https://wikis.nyu.edu/display/NYUHPC/Slurm+Tutorial#SlurmTutorial-SoftwareandEnvironmentModules).
module load rstudio/1.1.383
# Module ngsplot/20170524 was setup based on module r/intel/3.3.2, so please try to use R 3.3.2 on the server.
# NOTE: module swap may not work if the server does not have r/intel loaded by default.  In that case do "module load r/intel/3.3.2".  If that does not work, install r/intel.
module swap r/intel  r/intel/3.3.2
module load ngsplot/20170524
module load bedtools/intel/2.26.0
module load samtools/intel/1.6
module load python/intel/2.7.12
# Create bam file that will be used in ngs.plot.r NOTE: this bam file has to be hg19 not hg38 like everything else because ngs.plot.r has -G hg19 but not -G hg38.
# The first bedToBam not used which is why it is commented out, use the second one.
# bedToBam -i PARCLIP_QKI_Hafner2010c_hg19.bed -g human.hg19.genome > PARCLIP_QKI_Hafner2010c_hg19.bam
sortBed -i PARCLIP_QKI_Hafner2010c_hg19.bed > PARCLIP_QKI_Hafner2010c_hg19_sorted.bed
bedToBam -i PARCLIP_QKI_Hafner2010c_hg19_sorted.bed -g human.hg19.genome > PARCLIP_QKI_Hafner2010c_hg19_sorted.bam
samtools index PARCLIP_QKI_Hafner2010c_hg19_sorted.bam

# BEGIN NGSPLOT

# CREATE FIGURE 19 in thesis.
# -G is Genome name.
# -R is Genomic regions to plot.
#  PARCLIP_QKI_Hafner2010c_hg19_sorted.bam is the input file.
#  PARCLIP_QKI_Hafner2010c_hg19_sorted.genebody is the output file.
# -P is #CPUs to be used.
# -T is title of output graph.
# NOTE: the paths of singularity and ngsplot*.img may need to be changed to the ones on the particular server.
/share/apps/singularity/3.0.1/bin/singularity exec /beegfs/work/public/singularity/ngsplot-2.61.img ngs.plot.r -G hg19 -R genebody -C PARCLIP_QKI_Hafner2010c_hg19_sorted.bam -O PARCLIP_QKI_Hafner2010c_hg19_sorted.genebody -P 28 -T QKI

# CREATE FIGURE 15 in thesis.
# -G is Genome name.
# -R is Genomic regions to plot.
#  PARCLIP_QKI_Hafner2010c_hg19_sorted.bam is the input file.
#  PARCLIP_QKI_Hafner2010c_hg19_sorted.exon is the output file.
# -P is #CPUs to be used.
# -T is title of output graph.
# NOTE: the paths of singularity and ngsplot*.img may need to be changed to the ones on the particular server.
/share/apps/singularity/3.0.1/bin/singularity exec /beegfs/work/public/singularity/ngsplot-2.61.img ngs.plot.r -G hg19 -R exon -C PARCLIP_QKI_Hafner2010c_hg19_sorted.bam -O PARCLIP_QKI_Hafner2010c_hg19_sorted.exon -P 28 -T QKI

# https://github.com/shenlab-sinai/ngsplot
# https://github.com/shenlab-sinai/ngsplot/wiki/ProgramArguments101

# END NGSPLOT

# BEGIN BEDTOOLS CLOSEST

# Find middle of PAR-CLIP binding site (is usually 40nt wide) and make start and end of QKI binding region that number, so that bedtools closest does not return 0 for an overlap when the midpoint is not in the region.  The middle is supposed to be where the T to C conversion is that signifies QKI binding to the RNA, although I think the paper says it is not always exactly in the middle.
awk '{two=$2; $2=int((two+(($3-two)/2))); $3=int((two+(($3-two)/2))); print $0}' PARCLIP_QKI_Hafner2010c_hg38.bed > PARCLIP_QKI_Hafner2010c_hg38_centered.bed
# Replace all space with tab, because that is proper BED format.
sed -e "s/\ /\t/g" PARCLIP_QKI_Hafner2010c_hg38_centered.bed > temp$$
cat temp$$ > PARCLIP_QKI_Hafner2010c_hg38_centered.bed
rm -f temp$$

awk '$1 !~ /\_/' hg38_exons.bed > hg38_exons_noWeirdChr.bed
awk '$1 !~ /\_/' hg38_genes.bed > hg38_genes_noWeirdChr.bed
awk '$1 !~ /\_/' hg38_5p.bed > hg38_5p_noWeirdChr.bed
awk '$1 !~ /\_/' hg38_3p.bed > hg38_3p_noWeirdChr.bed
awk '$1 !~ /\_/' hsa_hg38_circRNA.bed > hsa_hg38_circRNA_noWeirdChr.bed
awk '$1 !~ /\_/' PARCLIP_QKI_Hafner2010c_hg38_centered.bed > PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr.bed
awk '$1 !~ /\_/' hg38_introns.bed > hg38_introns_noWeirdChr.bed
awk '$1 !~ /\_/' hg38_CDS.bed > hg38_CDS_noWeirdChr.bed

# This was part of a quick test to see if QKI binds more often near large or small circRNAs.
# Split the circRNA BED file into 100-3000nt and >3000nt.
awk '{if ($3-$2 >= 100 && $3-$2 <= 3000) {print $0}}' hsa_hg38_circRNA_noWeirdChr.bed > hsa_hg38_circRNA_noWeirdChr_100to3000.bed
awk '{if ($3-$2 > 3000) {print $0}}' hsa_hg38_circRNA_noWeirdChr.bed > hsa_hg38_circRNA_noWeirdChr_greaterthan3000.bed

# Sort the BED files by chromosome and other criteria.
sortBed -i hg38_exons_noWeirdChr.bed > hg38_exons_noWeirdChr_sorted.bed
sortBed -i hg38_genes_noWeirdChr.bed > hg38_genes_noWeirdChr_sorted.bed
sortBed -i hg38_5p_noWeirdChr.bed > hg38_5p_noWeirdChr_sorted.bed
sortBed -i hg38_3p_noWeirdChr.bed > hg38_3p_noWeirdChr_sorted.bed
sortBed -i hsa_hg38_circRNA_noWeirdChr.bed > hsa_hg38_circRNA_noWeirdChr_sorted.bed
sortBed -i hg38_introns_noWeirdChr.bed > hg38_introns_noWeirdChr_sorted.bed
sortBed -i hg38_CDS_noWeirdChr.bed > hg38_CDS_noWeirdChr_sorted.bed
sortBed -i cassette_exons_notstrict_nohead.bed > cassette_exons_notstrict_nohead_sorted.bed
sortBed -i hsa_hg38_circRNA_noWeirdChr_100to3000.bed > hsa_hg38_circRNA_noWeirdChr_100to3000_sorted.bed
sortBed -i hsa_hg38_circRNA_noWeirdChr_greaterthan3000.bed > hsa_hg38_circRNA_noWeirdChr_greaterthan3000_sorted.bed
sortBed -i PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr.bed > PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed

# BEGIN OVERLAP MEASUREMENT (WITHIN BEDTOOLS CLOSEST) FROM MIDDLE OF QKI BINDING SITE TO START OR END OF REGION

# This first step finds whether the QKI binding site (-a) is upstream (-), inside (0), or downstream (+) of the region of interest (-b).
bedtools closest -D b -s -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_genes_noWeirdChr_sorted.bed > hg38_QKI_closest_gene_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_5p_noWeirdChr_sorted.bed > hg38_QKI_closest_5p_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_3p_noWeirdChr_sorted.bed > hg38_QKI_closest_3p_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hsa_hg38_circRNA_noWeirdChr_sorted.bed > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_exons_noWeirdChr_sorted.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_introns_noWeirdChr_sorted.bed > hg38_QKI_closest_introns_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_CDS_noWeirdChr_sorted.bed > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b cassette_exons_notstrict_nohead_sorted.bed > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hsa_hg38_circRNA_noWeirdChr_100to3000_sorted.bed > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP.bed
bedtools closest -D b -s -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hsa_hg38_circRNA_noWeirdChr_greaterthan3000_sorted.bed > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP.bed

# Remove the part of the QKI PAR-CLIP that overlaps with exons, to be used with circRNA to see if that makes peaks.
awk '{if ($NF!=0) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}}' hg38_QKI_closest_exons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_notexons_noWeirdChr_PARCLIP.bed
# Remove all duplicates because the QKI PAR-CLIP will be fed back into another bedtools closest, unlike all other cases.
awk '{if (two!=$2 || three!=$3) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}
two=$2
three=$3}' hg38_QKI_closest_notexons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_notexons_noWeirdChr_PARCLIP.bed
rm -f temp$$
# don't have to do this, can skip to -io and -id
bedtools closest -D b -s -a hg38_QKI_closest_notexons_noWeirdChr_PARCLIP.bed -b hsa_hg38_circRNA_noWeirdChr_sorted.bed > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP.bed

# Remove the part of the QKI PAR-CLIP that overlaps with introns, to be used with 5p UTR and 3p UTR to see if that makes peaks.
awk '{if ($NF!=0) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}}' hg38_QKI_closest_introns_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_notintrons_noWeirdChr_PARCLIP.bed
# Remove all duplicates because the QKI PAR-CLIP will be fed back into another bedtools closest, unlike all other cases.
awk '{if (two!=$2 || three!=$3) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}
two=$2
three=$3}' hg38_QKI_closest_notintrons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_notintrons_noWeirdChr_PARCLIP.bed
rm -f temp$$

# Remove -1 error messages with .
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_gene_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_gene_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_5p_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_5p_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_3p_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_3p_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_circgene_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_exons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_exons_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_introns_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_introns_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_CDS_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP.bed
rm -f temp$$

# Actually I shouldn't have negated the upstream or downstream in the first step immediately below, because that gives the wrong sign.  I didn't bother changing all of the names of the files to remove neg from the filename.
# To fix this I reversed the sign in the second step immediately below "undo the previous upneg and downneg".
# See Figure 7B of the thesis this is the code that solves the problem in Figure 7A by doing Figure 7B.
# The first step immediately below is the second overall step (first step was bedtools closest) and it both removes upstream and downstream (not 0) values and measures the distance of QKI binding site to either the start or the end of the region, whichever is closer.
awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream < absdownstream) {print $0 "\t" (upstream)}}}' hg38_QKI_closest_exons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg.bed
awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream > absdownstream) {print $0 "\t" (downstream)}}}' hg38_QKI_closest_exons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg.bed

awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream < absdownstream) {print $0 "\t" (upstream)}}}' hg38_QKI_closest_5p_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_upneg.bed
awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream > absdownstream) {print $0 "\t" (downstream)}}}' hg38_QKI_closest_5p_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg.bed

awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream < absdownstream) {print $0 "\t" (upstream)}}}' hg38_QKI_closest_3p_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg.bed
awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream > absdownstream) {print $0 "\t" (downstream)}}}' hg38_QKI_closest_3p_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_downneg.bed

awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream < absdownstream) {print $0 "\t" (upstream)}}}' hg38_QKI_closest_gene_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_upneg.bed
awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream > absdownstream) {print $0 "\t" (downstream)}}}' hg38_QKI_closest_gene_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_downneg.bed

awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream < absdownstream) {print $0 "\t" (upstream)}}}' hg38_QKI_closest_circgene_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_upneg.bed
awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream > absdownstream) {print $0 "\t" (downstream)}}}' hg38_QKI_closest_circgene_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_downneg.bed

awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream < absdownstream) {print $0 "\t" (upstream)}}}' hg38_QKI_closest_introns_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg.bed
awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream > absdownstream) {print $0 "\t" (downstream)}}}' hg38_QKI_closest_introns_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_downneg.bed

awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream < absdownstream) {print $0 "\t" (upstream)}}}' hg38_QKI_closest_CDS_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_upneg.bed
awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream > absdownstream) {print $0 "\t" (downstream)}}}' hg38_QKI_closest_CDS_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_downneg.bed

awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream < absdownstream) {print $0 "\t" (upstream)}}}' hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_upneg.bed
awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream > absdownstream) {print $0 "\t" (downstream)}}}' hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_downneg.bed

awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream < absdownstream) {print $0 "\t" (upstream)}}}' hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP.bed > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_upneg.bed
awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream > absdownstream) {print $0 "\t" (downstream)}}}' hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP.bed > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_downneg.bed

awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream < absdownstream) {print $0 "\t" (upstream)}}}' hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP.bed > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_upneg.bed
awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream > absdownstream) {print $0 "\t" (downstream)}}}' hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP.bed > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_downneg.bed

awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream < absdownstream) {print $0 "\t" (upstream)}}}' hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_upneg.bed
awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream > absdownstream) {print $0 "\t" (downstream)}}}' hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_downneg.bed

# Look at uc001atl.2_cds_0_0_chr1_11982576_r in Integrative Genomics Viewer for possible circRNA formation.

# Undo the negation in the previous upneg and downneg (reverse sign of last column of BED file) to get the proper mathematical sign (- means the QKI is upstream of the start or end of the region and + means the QKI is downstream of the start or end of the region).
awk '$NF *= -1' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk '$NF *= -1' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk '$NF *= -1' hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk '$NF *= -1' hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk '$NF *= -1' hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk '$NF *= -1' hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk '$NF *= -1' hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk '$NF *= -1' hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk '$NF *= -1' hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk '$NF *= -1' hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk '$NF *= -1' hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk '$NF *= -1' hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk '$NF *= -1' hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk '$NF *= -1' hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk '$NF *= -1' hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk '$NF *= -1' hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk '$NF *= -1' hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk '$NF *= -1' hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk '$NF *= -1' hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk '$NF *= -1' hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk '$NF *= -1' hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk '$NF *= -1' hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

# Immediately below is another way to make all spaces tabs instead of using sed -e "s/\ /\t/g".
awk -v OFS="\t" '$1=$1' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk -v OFS="\t" '$1=$1' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk -v OFS="\t" '$1=$1' hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk -v OFS="\t" '$1=$1' hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk -v OFS="\t" '$1=$1' hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk -v OFS="\t" '$1=$1' hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk -v OFS="\t" '$1=$1' hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk -v OFS="\t" '$1=$1' hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk -v OFS="\t" '$1=$1' hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk -v OFS="\t" '$1=$1' hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk -v OFS="\t" '$1=$1' hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk -v OFS="\t" '$1=$1' hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk -v OFS="\t" '$1=$1' hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk -v OFS="\t" '$1=$1' hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk -v OFS="\t" '$1=$1' hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk -v OFS="\t" '$1=$1' hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk -v OFS="\t" '$1=$1' hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk -v OFS="\t" '$1=$1' hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk -v OFS="\t" '$1=$1' hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk -v OFS="\t" '$1=$1' hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_downneg.bed
rm -f temp$$

awk -v OFS="\t" '$1=$1' hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_upneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_upneg.bed
rm -f temp$$
awk -v OFS="\t" '$1=$1' hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
cat temp$$ > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_downneg.bed
rm -f temp$$

# Remove exon_0 and intron_0.
# Used later for histograms: hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_exonstart.bed and hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_exonend.bed
awk '$10 !~ /exon_0/' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_nofirstexon.bed
# awk '$10 !~ /intron_0/' hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg.bed > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg_nointron0.bed

# Select only exon_0/ so those can be cat with the 5' UTR file to make the 5' UTR - first exon junction histogram.
awk '$10 ~ /exon_0/' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_firstexon.bed

# Reverse file to help in finding the last exon.
tac hg38_exons.bed > hg38_exons_rev.bed

# Make lastexons.awk executable.
chmod +x lastexons.awk
# Find last exons.
cat hg38_exons_rev.bed | ./lastexons.awk > hg38_exons_rev_lastexons.bed
# I am not copying and pasting the body of lastexons.awk in the terminal because it loses some of my text.
# I am not making the awk code a single line here to solve that because then it is too hard to understand and comment.

# Print out geneexon id.
awk '{print $4}' hg38_exons_rev_lastexons.bed > hg38_exons_rev_lastexons_id.bed

# Remove last exons from hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg.bed.
# Make removelastexons.py executable.
chmod +x removelastexons.py
# Make onlylastexons.py executable.
chmod +x onlylastexons.py
# Remove last exons from the file containing the QKI that is closer to the end of the exon.
# The immediately below python script creates hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_nolastexon.bed.
python removelastexons.py
# Select only last exons so those can be cat with the 3' UTR file to make the last exon - 3' UTR junction histogram.
# The immediately below python script creates hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_onlylastexon.bed.
python onlylastexons.py

# Sort file.
sort hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed
sort hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed
sort hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_nofirstexon.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_nofirstexon_sorted.bed
sort hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_nolastexon.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_nolastexon_sorted.bed
sort hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_firstexon.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_firstexon_sorted.bed
sort hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_onlylastexon.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_lastexon_sorted.bed
sort hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_upneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed
sort hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed
sort hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed
sort hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_downneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed
sort hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_upneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed
sort hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_downneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed
sort hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_upneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed
sort hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_downneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed
sort hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed
sort hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_downneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed
sort hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_upneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed
sort hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_downneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed
sort hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_upneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed
sort hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_downneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed
sort hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_upneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_upneg_sorted.bed
sort hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_downneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_downneg_sorted.bed
sort hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_upneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_upneg_sorted.bed
sort hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_downneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_downneg_sorted.bed
sort hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_upneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed
sort hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_downneg.bed -k1,1 -k10,10 -k2,2n > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed

# Prints the last column of the bedtools closest file, which is the distance to either the start (upneg) or end (downneg) of the region.
awk '{print $NF}' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_nofirstexon_sorted.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_nofirstexon_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_nolastexon_sorted.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_nolastexon_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_firstexon_sorted.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_firstexon_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_lastexon_sorted.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_lastexon_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_upneg_sorted.bed > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_upneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_downneg_sorted.bed > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_downneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_upneg_sorted.bed > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_upneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_downneg_sorted.bed > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_downneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_upneg_sorted.bed > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed
awk '{print $NF}' hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_downneg_sorted.bed > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed

# Find the frequency of distances (not used in the project), e.g. if a distance of 10 happened 3 times it would print two columns: 3	10
sort hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist_freq.bed
sort hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist_freq.bed
sort hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_nofirstexon_sorted_dist.bed | uniq -c > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_nofirstexon_sorted_dist_freq.bed
sort hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_nolastexon_sorted_dist.bed | uniq -c > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_nolastexon_sorted_dist_freq.bed
sort hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_firstexon_sorted_dist.bed | uniq -c > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_firstexon_sorted_dist_freq.bed
sort hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_lastexon_sorted_dist.bed | uniq -c > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_lastexon_sorted_dist_freq.bed
sort hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist_freq.bed
sort hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist_freq.bed
sort hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist_freq.bed
sort hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist_freq.bed
sort hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist_freq.bed
sort hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist_freq.bed
sort hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist_freq.bed
sort hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist_freq.bed
sort hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist_freq.bed
sort hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist_freq.bed
sort hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist_freq.bed
sort hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist_freq.bed
sort hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist_freq.beeg_sorted_dist_freq.bed
sort hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist_freq.bed
sort hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_upneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_upneg_sorted_dist_freq.bed
sort hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_downneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_downneg_sorted_dist_freq.bed
sort hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_upneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_upneg_sorted_dist_freq.bed
sort hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_downneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_downneg_sorted_dist_freq.bed
sort hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist_freq.bed
sort hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed | uniq -c > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist_freq.bed

# END OVERLAP MEASUREMENT (WITHIN BEDTOOLS CLOSEST) FROM MIDDLE OF QKI BINDING SITE TO START OR END OF REGION

# BEGIN -io -iu -id (WITHIN BEDTOOLS CLOSEST), which will be combined with upneg and downneg to create the histograms

bedtools closest -D b -s -io -id -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_genes_noWeirdChr_sorted.bed > hg38_QKI_closest_io_id_gene_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -iu -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_genes_noWeirdChr_sorted.bed > hg38_QKI_closest_io_iu_gene_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -id -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hsa_hg38_circRNA_noWeirdChr_sorted.bed > hg38_QKI_closest_io_id_circgene_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -iu -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hsa_hg38_circRNA_noWeirdChr_sorted.bed > hg38_QKI_closest_io_iu_circgene_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -id -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_exons_noWeirdChr_sorted.bed > hg38_QKI_closest_io_id_exons_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -iu -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_exons_noWeirdChr_sorted.bed > hg38_QKI_closest_io_iu_exons_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -id -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_5p_noWeirdChr_sorted.bed > hg38_QKI_closest_io_id_5p_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -iu -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_5p_noWeirdChr_sorted.bed > hg38_QKI_closest_io_iu_5p_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -id -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_3p_noWeirdChr_sorted.bed > hg38_QKI_closest_io_id_3p_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -iu -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_3p_noWeirdChr_sorted.bed > hg38_QKI_closest_io_iu_3p_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -id -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_introns_noWeirdChr_sorted.bed > hg38_QKI_closest_io_id_introns_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -iu -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_introns_noWeirdChr_sorted.bed > hg38_QKI_closest_io_iu_introns_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -id -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_CDS_noWeirdChr_sorted.bed > hg38_QKI_closest_io_id_CDS_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -iu -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_CDS_noWeirdChr_sorted.bed > hg38_QKI_closest_io_iu_CDS_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -id -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b cassette_exons_notstrict_nohead_sorted.bed > hg38_QKI_closest_io_id_cassette_exons_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -iu -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b cassette_exons_notstrict_nohead_sorted.bed > hg38_QKI_closest_io_iu_cassette_exons_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -id -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hsa_hg38_circRNA_noWeirdChr_100to3000_sorted.bed > hg38_QKI_closest_io_id_circgene_noWeirdChr_100to3000_PARCLIP.bed
bedtools closest -D b -s -io -iu -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hsa_hg38_circRNA_noWeirdChr_100to3000_sorted.bed > hg38_QKI_closest_io_iu_circgene_noWeirdChr_100to3000_PARCLIP.bed
bedtools closest -D b -s -io -id -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hsa_hg38_circRNA_noWeirdChr_greaterthan3000_sorted.bed > hg38_QKI_closest_io_id_circgene_noWeirdChr_greaterthan3000_PARCLIP.bed
bedtools closest -D b -s -io -iu -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hsa_hg38_circRNA_noWeirdChr_greaterthan3000_sorted.bed > hg38_QKI_closest_io_iu_circgene_noWeirdChr_greaterthan3000_PARCLIP.bed
bedtools closest -D b -s -io -id -a hg38_QKI_closest_notexons_noWeirdChr_PARCLIP.bed -b hsa_hg38_circRNA_noWeirdChr_sorted.bed > hg38_QKI_closest_io_id_circgene_notexons_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -iu -a hg38_QKI_closest_notexons_noWeirdChr_PARCLIP.bed -b hsa_hg38_circRNA_noWeirdChr_sorted.bed > hg38_QKI_closest_io_iu_circgene_notexons_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -id -a hg38_QKI_closest_notexons_noWeirdChr_PARCLIP.bed -b cassette_exons_notstrict_nohead_sorted.bed > hg38_QKI_closest_io_id_cassette_exons_notexons_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -iu -a hg38_QKI_closest_notexons_noWeirdChr_PARCLIP.bed -b cassette_exons_notstrict_nohead_sorted.bed > hg38_QKI_closest_io_iu_cassette_exons_notexons_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -iu -a hg38_QKI_closest_notintrons_noWeirdChr_PARCLIP.bed -b hg38_5p_noWeirdChr_sorted.bed > hg38_QKI_closest_io_iu_5p_notintrons_noWeirdChr_PARCLIP.bed
bedtools closest -D b -s -io -id -a hg38_QKI_closest_notintrons_noWeirdChr_PARCLIP.bed -b hg38_3p_noWeirdChr_sorted.bed > hg38_QKI_closest_io_id_3p_notintrons_noWeirdChr_PARCLIP.bed

# Remove -1 error messages with .
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_id_gene_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_id_gene_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_iu_gene_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_iu_gene_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_id_circgene_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_id_circgene_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_iu_circgene_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_iu_circgene_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_id_exons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_id_exons_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_iu_exons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_iu_exons_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_id_5p_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_id_5p_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_iu_5p_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_iu_5p_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_id_3p_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_id_3p_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_iu_3p_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_iu_3p_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_id_introns_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_id_introns_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_iu_introns_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_iu_introns_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_id_CDS_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_id_CDS_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_iu_CDS_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_iu_CDS_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_id_cassette_exons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_id_cassette_exons_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_iu_cassette_exons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_iu_cassette_exons_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_id_circgene_noWeirdChr_100to3000_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_id_circgene_noWeirdChr_100to3000_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_iu_circgene_noWeirdChr_100to3000_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_iu_circgene_noWeirdChr_100to3000_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_id_circgene_noWeirdChr_greaterthan3000_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_id_circgene_noWeirdChr_greaterthan3000_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_iu_circgene_noWeirdChr_greaterthan3000_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_iu_circgene_noWeirdChr_greaterthan3000_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_id_circgene_notexons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_id_circgene_notexons_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_iu_circgene_notexons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_iu_circgene_notexons_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_id_cassette_exons_notexons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_id_cassette_exons_notexons_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_iu_cassette_exons_notexons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_iu_cassette_exons_notexons_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_iu_5p_notintrons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_iu_5p_notintrons_noWeirdChr_PARCLIP.bed
rm -f temp$$
awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_io_id_3p_notintrons_noWeirdChr_PARCLIP.bed > temp$$
cat temp$$ > hg38_QKI_closest_io_id_3p_notintrons_noWeirdChr_PARCLIP.bed
rm -f temp$$

# Print distances upstream (io_id) or downstream (io_iu) of the region.
awk '{print $NF}' hg38_QKI_closest_io_id_gene_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_id_gene_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_iu_gene_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_iu_gene_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_id_circgene_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_id_circgene_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_iu_circgene_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_iu_circgene_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_id_exons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_id_exons_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_iu_exons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_iu_exons_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_id_5p_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_id_5p_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_iu_5p_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_iu_5p_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_id_3p_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_id_3p_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_iu_3p_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_iu_3p_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_id_introns_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_id_introns_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_iu_introns_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_iu_introns_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_id_CDS_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_id_CDS_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_iu_CDS_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_iu_CDS_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_id_cassette_exons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_id_cassette_exons_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_iu_cassette_exons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_iu_cassette_exons_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_id_circgene_noWeirdChr_100to3000_PARCLIP.bed > hg38_QKI_closest_io_id_circgene_noWeirdChr_100to3000_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_iu_circgene_noWeirdChr_100to3000_PARCLIP.bed > hg38_QKI_closest_io_iu_circgene_noWeirdChr_100to3000_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_id_circgene_noWeirdChr_greaterthan3000_PARCLIP.bed > hg38_QKI_closest_io_id_circgene_noWeirdChr_greaterthan3000_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_iu_circgene_noWeirdChr_greaterthan3000_PARCLIP.bed > hg38_QKI_closest_io_iu_circgene_noWeirdChr_greaterthan3000_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_id_circgene_notexons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_id_circgene_notexons_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_iu_circgene_notexons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_iu_circgene_notexons_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_id_cassette_exons_notexons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_id_cassette_exons_notexons_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_iu_cassette_exons_notexons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_iu_cassette_exons_notexons_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_iu_5p_notintrons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_iu_5p_notintrons_noWeirdChr_PARCLIP_dist.bed
awk '{print $NF}' hg38_QKI_closest_io_id_3p_notintrons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_io_id_3p_notintrons_noWeirdChr_PARCLIP_dist.bed

# Find the frequency of distances (not used in the project), e.g. if a distance of 10 happened 3 times it would print two columns: 3	10
sort hg38_QKI_closest_io_id_gene_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_id_gene_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_iu_gene_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_iu_gene_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_id_circgene_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_id_circgene_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_iu_circgene_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_iu_circgene_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_id_exons_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_id_exons_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_iu_exons_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_iu_exons_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_id_5p_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_id_5p_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_iu_5p_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_iu_5p_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_id_3p_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_id_3p_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_iu_3p_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_iu_3p_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_id_introns_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_id_introns_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_iu_introns_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_iu_introns_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_id_CDS_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_id_CDS_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_iu_CDS_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_iu_CDS_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_id_cassette_exons_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_id_cassette_exons_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_iu_cassette_exons_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_iu_cassette_exons_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_id_circgene_noWeirdChr_100to3000_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_id_circgene_noWeirdChr_100to3000_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_iu_circgene_noWeirdChr_100to3000_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_iu_circgene_noWeirdChr_100to3000_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_id_circgene_noWeirdChr_greaterthan3000_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_id_circgene_noWeirdChr_greaterthan3000_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_iu_circgene_noWeirdChr_greaterthan3000_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_iu_circgene_noWeirdChr_greaterthan3000_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_id_circgene_notexons_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_id_circgene_notexons_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_iu_circgene_notexons_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_iu_circgene_notexons_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_id_cassette_exons_notexons_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_id_cassette_exons_notexons_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_iu_cassette_exons_notexons_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_iu_cassette_exons_notexons_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_iu_5p_notintrons_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_iu_5p_notintrons_noWeirdChr_PARCLIP_dist_freq.bed
sort hg38_QKI_closest_io_id_3p_notintrons_noWeirdChr_PARCLIP_dist.bed | uniq -c > hg38_QKI_closest_io_id_3p_notintrons_noWeirdChr_PARCLIP_dist_freq.bed

# END -io -iu -id (WITHIN BEDTOOLS CLOSEST), which will be combined with upneg and downneg to create the histograms

# END BEDTOOLS CLOSEST

# BEGIN biologically relevant histograms

# BEGIN OLD there is a better version immediately below

cat hg38_QKI_closest_io_id_exons_noWeirdChr_PARCLIP_dist.bed hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_regionstart.bed
cat hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_io_iu_exons_noWeirdChr_PARCLIP_dist.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_regionend.bed

cat hg38_QKI_closest_io_id_introns_noWeirdChr_PARCLIP_dist.bed hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_regionstart.bed
cat hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_io_iu_introns_noWeirdChr_PARCLIP_dist.bed > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_regionend.bed

# END OLD there is a better version immediately below

# exonstart and exonend are better than the four above because I combine inside intron and inside exon without needing io_iu and io_id and I get rid of the 5' UTR - exon junction and exon - 3' UTR junction which the above 4 do not do.
# I also don't need an intron start and end because the exon start and exon end also include the end of the intron and the start of the intron.
# CREATE FIGURE 14 (left) in thesis
cat hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_nofirstexon_sorted_dist.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_exonstart.bed
# CREATE FIGURE 14 (right) in thesis
cat hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_nolastexon_sorted_dist.bed hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_exonend.bed

# BEGIN cat io_id with upneg and downneg with io_iu (this is Figure 8B of the thesis)

# CREATE FIGURE 18 (left) in thesis.
cat hg38_QKI_closest_io_id_5p_noWeirdChr_PARCLIP_dist.bed hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_regionstart.bed
cat hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_io_iu_5p_noWeirdChr_PARCLIP_dist.bed > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_regionend.bed
cat hg38_QKI_closest_io_id_3p_noWeirdChr_PARCLIP_dist.bed hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_regionstart.bed
# CREATE FIGURE 18 (right) in thesis.
cat hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_io_iu_3p_noWeirdChr_PARCLIP_dist.bed > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_regionend.bed
cat hg38_QKI_closest_io_id_gene_noWeirdChr_PARCLIP_dist.bed hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_regionstart.bed
cat hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_io_iu_gene_noWeirdChr_PARCLIP_dist.bed > hg38_QKI_closest_gene_noWeirdChr_PARCLIP_overlap_regionend.bed
cat hg38_QKI_closest_io_id_circgene_noWeirdChr_PARCLIP_dist.bed hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_regionstart.bed
cat hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_io_iu_circgene_noWeirdChr_PARCLIP_dist.bed > hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_regionend.bed
cat hg38_QKI_closest_io_id_CDS_noWeirdChr_PARCLIP_dist.bed hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_regionstart.bed
cat hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_io_iu_CDS_noWeirdChr_PARCLIP_dist.bed > hg38_QKI_closest_CDS_noWeirdChr_PARCLIP_overlap_regionend.bed
cat hg38_QKI_closest_io_id_cassette_exons_noWeirdChr_PARCLIP_dist.bed hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_regionstart.bed
cat hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_io_iu_cassette_exons_noWeirdChr_PARCLIP_dist.bed > hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_regionend.bed
cat hg38_QKI_closest_io_id_circgene_noWeirdChr_100to3000_PARCLIP_dist.bed hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_regionstart.bed
cat hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_io_iu_circgene_noWeirdChr_100to3000_PARCLIP_dist.bed > hg38_QKI_closest_circgene_noWeirdChr_100to3000_PARCLIP_overlap_regionend.bed
cat hg38_QKI_closest_io_id_circgene_noWeirdChr_greaterthan3000_PARCLIP_dist.bed hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_regionstart.bed
cat hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_io_iu_circgene_noWeirdChr_greaterthan3000_PARCLIP_dist.bed > hg38_QKI_closest_circgene_noWeirdChr_greaterthan3000_PARCLIP_overlap_regionend.bed
# NOTE: DO NOT cat hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed or hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed, because that would remove exons from the circRNA.  Use the two cat commands immediately below.
# CREATE FIGURE 20 (left) in thesis.
cat hg38_QKI_closest_io_id_circgene_notexons_noWeirdChr_PARCLIP_dist.bed hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_regionstart.bed
# CREATE FIGURE 20 (right) in thesis.
cat hg38_QKI_closest_circgene_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_io_iu_circgene_notexons_noWeirdChr_PARCLIP_dist.bed > hg38_QKI_closest_circgene_notexons_noWeirdChr_PARCLIP_overlap_regionend.bed
# CREATE FIGURE 16 (left) in thesis.
cat hg38_QKI_closest_io_id_cassette_exons_notexons_noWeirdChr_PARCLIP_dist.bed hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_cassette_exons_notexons_noWeirdChr_PARCLIP_overlap_regionstart.bed
# CREATE FIGURE 16 (right) in thesis.
cat hg38_QKI_closest_cassette_exons_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_io_iu_cassette_exons_notexons_noWeirdChr_PARCLIP_dist.bed > hg38_QKI_closest_cassette_exons_notexons_noWeirdChr_PARCLIP_overlap_regionend.bed
# Two of the BED files immediately below had their introns removed (called notintrons) so that short exons downstream of 5' UTR or upstream of 3' UTR would not spill over into intronic regions if the histogram distance of 300nt was larger than the exon distance.  If that was not the case, the 5' UTR - exon 1 junction could spill over into introns past exon 1 and the last exon - 3' UTR junction could have an intron to the left of the last exon if the exons are short.  See Figure 8B of my thesis for why short regions are a problem when making the histograms.
# CREATE FIGURE 17 (left) in thesis.
cat hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_io_iu_5p_notintrons_noWeirdChr_PARCLIP_dist.bed > hg38_QKI_closest_5p_notintrons_noWeirdChr_PARCLIP_overlap_regionend.bed
# CREATE FIGURE 17 (right) in thesis.
cat hg38_QKI_closest_io_id_3p_notintrons_noWeirdChr_PARCLIP_dist.bed hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_3p_notintrons_noWeirdChr_PARCLIP_overlap_regionstart.bed

# END cat io_id with upneg and downneg with io_iu (this is Figure 8B of the thesis)

# BEGIN Double-check of Figure 17

# exonstart and exonend removed the 5' UTR - exon junction and exon - 3' UTR junction properly.  The UTR regions of the alternate path histogram is identical to Figure 17 (as expected because both paths did the same thing: bedtools closest to find QKI binding sites overlapping with the UTR, then split the UTR in half), but the first and last exon regions of this histogram is not the same (bedtools closest upstream and downstream was used for Figure 17 and trimmed BED files containing only the first exon or only the last exon were used for the two figures immediately below called 'Should be the same as FIGURE 17'").  A possible explanation why is that the 5' UTR is part of the first exon and the 3' UTR is part of the last exon in the exons, 5' UTR, or 3' UTR BED files I got from https://genome.ucsc.edu/cgi-bin/hgTables.  See the pictured immature RNA in Gene_structure_eukaryote_2_annotated.png, the first exon includes the 5' UTR and the last exon includes the 3' UTR.
# https://en.wikiversity.org/wiki/WikiJournal_of_Medicine/Eukaryotic_and_prokaryotic_gene_structure
# The immediately below should be the same as FIGURE 17 (left) in thesis (only in "Hg38_QKI_closest_dist_freqV9.html" not published in my thesis)
cat hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_downneg_sorted_dist.bed hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_firstexon_sorted_dist.bed > hg38_QKI_closest_5p_noWeirdChr_PARCLIP_overlap_regionend_firstexon.bed
# The immediately below should be the same as FIGURE 17 (right) in thesis (only in "Hg38_QKI_closest_dist_freqV9.html" not published in my thesis)
cat hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_upneg_lastexon_sorted_dist.bed hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_3p_noWeirdChr_PARCLIP_overlap_regionstart_lastexon.bed

# END Double-check of Figure 17

# END biologically relevant histograms

# BEGIN count how many QKI sites are within 300nt of junction

# never done
# awk '($11>0 && $11<0.05){ ++count } END{ print count }' file

# END count how many QKI sites are within 300nt of junction

##############################################################
# END OF SECTION 3: ONE-COMMAND SBATCH SCRIPT EXECUTABLE LINES
##############################################################


# AUTHOR:
# Thomas H. Lipscomb
# In partial fulfillment of New York University 2018 MS thesis.

# REPORTING BUGS:
# https://groups.google.com/forum/#!forum/bedtools-discuss
# post titled "How to make distance histograms of RNA-binding protein binding frequency near all types of pre-RNA junctions (intron-exon and exon-intron, etc.)"
# I may not be able to assist in debugging because I will no longer have access to the server I used to generate the example output directory "afterRunningBatchFile".

# COPYRIGHT: freeware copyleft