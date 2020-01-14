To plot distance histograms of protein binding sites to all types of pre-RNA junctions (intron-exon and exon-intron, etc.), closestBed (also called bedtools closest) requires post-processing to isolate pertinent information, limit results to flanking regions, and not double-count with the histogram if the regions are short. Also ngs.plot.r is used to create exon and genebody distance graphs.

For instructions and a one-command batch script working example, open the preRNAbinding folder then open the batchFile folder and open preRNAbinding.s

"batchFile" contains the batch script "preRNAbinding.s" and the input files needed to generate output data files, after which the contents of directory "batchFile" should be the same as the contents of "afterRunningBatchFile".

Step by step instructions to run "preRNAbinding.s" are self-documented in comments in the SYNOPSIS and PREREQUISITES sections.  Batch file "preRNAbinding.s" post-processes the files in directory "batchFile" to create the data files that it passes to R Studio to plot the distance histograms.  "preRNAbinding.s" also runs ngs.plot.r to create exon and genebody RNA binding protein binding frequency graphs.

To create the histograms: After running "preRNAbinding.s", R Studio must be manually run to use the "preRNAbinding.s" output files to produce the HTML file containing the PNG histograms.  Step by step instructions are documented in "preRNAbinding.s".


Background: The plots of these calculations are included and explained in my masters thesis "Lipscomb_Thomas_Thesis.pdf" included in directory "batchFile".  To understand the structure of pre-RNA see "Gene_structure_eukaryote_2_annotated.png" in "batchFile".


BEDTOOLS LIMITATIONS AND HOW SBATCH FILE "preRNAbinding.s" OVERCOMES THEM:

• BEDTOOLS CLOSEST LIMITATIONS:

Limitation 1: Bedtools closest returns negative numbers for RNA binding protein binding sites upstream of a region and positive numbers for RNA binding protein binding sites downstream of a region (e.g. an exon and the regions upstream and downstream), but does not limit the results to only RNA binding protein binding sites within the regions flanking the center region (e.g. intron-exon-intron, the intronic regions flank the center exon region), which is what is needed.  If the nearest RNA binding protein binding site is beyond the flanking region then that protein binding site is erroneously mixed in with the flanking region results (see Figure 7A of the thesis PDF "Lipscomb_Thomas_Thesis.pdf" in the batchFile folder).  Also, bedtools closest returns a zero for distance if the RNA binding protein site and preRNA region overlap, instead of the needed distance to either the start or end of region, whichever is closer.

Limitation 2: Erroneous double-counting of QKI binding sites inside the region can happen if the region is short (e.g. exons can be short), causing the two junction histograms to overlap (see Figure 8A of the thesis PDF in the batchFile folder).


• HOW BEDTOOLS CLOSEST LIMITATIONS ARE OVERCOME by "preRNAbinding.s":

Solution to Limitation 1: When bedtools closest returns zero, signifying that the RNA binding protein bound within the region (e.g. the exon or the flanking introns), AWK code provided in "preRNAbinding.s" measures the distance to the start or the end of the region, whichever is closer (see Figure 7B of the thesis PDF in the batchFile folder).  Looking within adjacent regions prevents erroneous inclusion of more distant protein binding sites.

Solution to Limitation 2: Split the center region in half and split the flanking regions in half.  This prevents erroneous double-counting, because no matter what the nucleotide distance on the x-axis of the junction histogram is, the junction histogram is given only the distances that were closer to either the start or the end of the center region, whichever the histogram includes (see Figure 8B of the thesis PDF in the batchFile folder).  This solution fully applied to intron-exon and exon-intron histograms.  This solution only partially applied to the other histograms, which split only the center region in half because the size of the flanking regions is unknown.


FUTURE RESEARCH:

Future Research 1: Fix edge case for Figure 16 and Figure 20 where the cassette exon or circRNA is the first or last exon, causing the flanking regions to not be introns even though they are labeled as such in Figure 16 and Figure 20.  If it is the first exon, then upstream actually is the 5' UTR not the upstream intron and if it is the last exon then downstream actually is the 3' UTR not the downstream intron.  "bedtools closest" was used to find the upstream and downstream regions for Figure 16 and Figure 20, but to fix this edge case more AWK and python code may be required to remove the cassette exons and circRNAs that are the first or last exon and make separate figures for those edge cases.

Future Research 2: The complete solutions to Limitation 1 and Limitation 2 could be extended from only intron-exon and exon-intron junction histograms to also be the complete solutions for the cassette exons and circRNA instead of only being partial solutions for those two.  Instead of using "bedtools closest" to find upstream and downstreasm RNA binding protein binding sites without being able to limit them to the flanking introns (Limitation 1), code might be written to try to include only introns that have a start or end position that is the same as the start or end of the cassette exon or circRNA, thereby identifying them as flanking introns, but first one should check if there are edge cases where the start or end of the intron is not exactly the same location as the start or end of the cassette exon (e.g. they are separated by a few nucleotides).

Future Research 3: Although Figure 17 is correct it does not agree with the histogram in the R Studio "Hg38_QKI_closest_dist_freqV9.html" file called: "Should be the same as Figure 17".  The 5' UTR and 3' UTR BED files should be checked to see whether or not they include the first and last exon (respectively) as "Gene_structure_eukaryote_2_annotated.png" suggests they might.


AUTHOR:

Thomas H. Lipscomb
In partial fulfillment of New York University 2018 MS thesis.

REPORTING BUGS:

https://groups.google.com/forum/#!forum/bedtools-discuss
post titled "Answer: Making distance histograms of protein binding frequency near all types of pre-RNA junctions".

NOTES: 
• INCORRECT POST TITLE IN COMMENTS: The executable files preRNAbinding.s, lastexons.awk, onlylastexons.py, and removelastexons.py have an incorrect post title "How to make distance histograms of RNA-binding protein binding frequency near all types of pre-RNA junctions (intron-exon and exon-intron, etc.)" because I did not find out that title is too long to post on bedtools-discuss before I lost access to the server I ran that code on.
• DEPERSONALIZATION: After testing everything with my personal email and personal directory, in preRNAbinding.s I changed my email to the metadata [youremail] and my directory to contain the metadata [yourusername], as follows:
First change:
#SBATCH --mail-user=[youremail]
Second change:
#SBATCH -e /home/[yourusername]/preRNAbinding/batchFile/preRNAbinding.err
#SBATCH -o /home/[yourusername]/preRNAbinding/batchFile/preRNAbinding.out

I may not be able to assist in debugging because I will no longer have access to the server I used to generate the example output directory "afterRunningBatchFile".

COPYRIGHT: freeware copyleft
