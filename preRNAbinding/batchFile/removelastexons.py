# NAME: removelastexons.py - take bedtools closest BED file that has QKI binding sites closer to the start of the exon and remove the lines that contain the last exon of a gene

# SYNOPSIS: python removelastexons.py
# Named input file 1: [reversed file containing ID of the last exons of all genes], from column 4 of hg38_exons_rev_lastexons.bed.
# Named input file 2: [measure of the negative distance from the end of the exon to the binding site of interest inside of the exon], a result of bedtools closest to find only binding sites inside of the exon (e.g. QKI binding sites inside the exon) and an awk script that measures the distance from only the sites that are closer to the end of the exon.
# Named output file: [Named input file 2 trimmed so binding sites in the last exon of genes are not included], because the last exon is right before the 3' UTR, not an intron, so it should not be included in the exon-intron junction histogram.

# DESCRIPTION:
# Stdin: none
# Stdout: none
# Named input file 1: hg38_exons_rev_lastexons_id.bed
# Named input file 2: hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg.bed
# Named output file: hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_nolastexon.bed
# Output file used in "cat hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_nolastexon_sorted_dist.bed hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_upneg_sorted_dist.bed > hg38_QKI_closest_introns_noWeirdChr_PARCLIP_overlap_exonend.bed", which is used to create the histogram of the end of the exon.

# USED BY: preRNAbinding.s
# USES: none

# EXAMPLE:

# # Go to http://dorina.mdc-berlin.de/regulators
# # Select assembly: hg19
# # Download http://dorina.mdc-berlin.de/api/v1.0/download/regulator/hg19/PARCLIP_QKI_Hafner2010c_hg19

# # Go to https://genome.ucsc.edu/cgi-bin/hgTables
# # Select assembly: hg38
# # Select track: GENCODE V24
# # Select table: knownGene
# # Select output format: BED - browser extensible data
# # Press get output
# # Select option Exons
# # Press get BED
# # Download hg38_exons.bed

# $ srun --cpus-per-task=20 --mem=62GB --time=04:00:00 --pty /bin/bash
# # Make removelastexons.py executable in preRNAbinding.s.
# $ chmod +x removelastexons.py
# # Remove chromosomes that contain an underscore because those are not standard.  Don't need to remove mitochondrial chromosome (chrM) even though QKI is nuclear, so that is commented out.
# $ awk '$1 !~ /\_/' PARCLIP_QKI_Hafner2010c_hg38.bed > PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr.bed
# $ sortBed -i PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr.bed > PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed
# # START OVERLAP MEASUREMENT FROM MIDDLE OF QKI BINDING SITE TO START OR END OF REGION
# $ bedtools closest -D b -s -a PARCLIP_QKI_Hafner2010c_hg38_noWeirdChr_sorted.bed -b hg38_exons_noWeirdChr_sorted.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP.bed

# # Remove -1 error messages with .
# $ awk '!/\.\t-1$/ {print $0}' hg38_QKI_closest_exons_noWeirdChr_PARCLIP.bed > temp$$
# $ cat temp$$ > hg38_QKI_closest_exons_noWeirdChr_PARCLIP.bed
# $ rm -f temp$$

# # Measure the distance from the start of the exon to the QKI inside of exons if that is closer than the distance of the QKI to the end of the exon.
# $ awk '{if ($NF == 0) {upstream = $8-$2; downstream = $9-$3; absupstream = upstream < 0 ? -upstream : upstream; absdownstream = downstream < 0 ? -downstream : downstream; if (absupstream > absdownstream) {print $0 "\t" (downstream)}}}' hg38_QKI_closest_exons_noWeirdChr_PARCLIP.bed > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg.bed

# Undo the previous upneg and downneg.
# $ awk '$NF *= -1' hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg.bed > temp$$
# $ cat temp$$ > hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg.bed
# $ rm -f temp$$

# # Reverse the file to help in finding the last exon is done in preRNAbinding.s.
# $ tac hg38_exons.bed > hg38_exons_rev.bed
# # Find last exons is continued from preRNAbinding.s.
# $ cat hg38_exons_rev.bed | ./lastexons.awk > hg38_exons_rev_lastexons.bed
# # Print out geneexon id is continued from preRNAbinding.s.
# $ awk '{print $4}' hg38_exons_rev_lastexons.bed > hg38_exons_rev_lastexons_id.bed
# $ module load python/intel/2.7.12
# $ python removelastexons.py

exons_to_remove = open("hg38_exons_rev_lastexons_id.bed","r").readlines()
ends_of_exons = open("hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg.bed","r").readlines()

# exons_to_remove = open("test_lastexons.bed","r").readlines()
# ends_of_exons = open("test_exon_downneg.bed","r").readlines()

found = False

# https://www.reddit.com/r/learnpython/comments/2emugh/how_do_i_delete_lines_in_a_text_file/?st=jf71ooay&sh=2251d309
# infile = open('input.txt','r').readlines()
# https://stackoverflow.com/questions/25044938/in-python-how-can-i-print-lines-that-do-not-contain-a-certain-string-rather-th?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
with open('hg38_QKI_closest_exons_noWeirdChr_PARCLIP_overlap_downneg_nolastexon.bed','w') as outfile:
    for lineQKI in ends_of_exons:
        found = False
        for lineRem in exons_to_remove:
            # print found
            # print lineQKI
            # print lineRem
            lineRemR = lineRem.rstrip()
            # print lineRemR
            if lineRemR in lineQKI:
            # if -1 != lineQKI.find(lineRem):
                found = True
                # print "XXXXXXXXXXXXXXXXXXXXX ", found
		# if this is not the last exon, print the line
        if found == False:
            outfile.write(lineQKI)

# AUTHOR:
# Thomas H. Lipscomb
# New York University 2018 MS thesis

# REPORTING BUGS:
# https://groups.google.com/forum/#!forum/bedtools-discuss
# post titled "How to make distance histograms of RNA-binding protein binding frequency near all types of pre-RNA junctions (intron-exon and exon-intron, etc.)"
# I may not be able to assist in debugging because I will no longer have access to the server I used to generate the example output directory "afterRunningBatchFile".

# COPYRIGHT: freeware copyleft