# NAME: lastexons.awk - awk program to print out BED file containing only last exons of genes (genes can have more than one exon)

# SYNOPSIS: cat [input reversed exon BED file] | ./lastexons.awk > [output BED file of only last exons of genes]

# DESCRIPTION:
# Stdin: hg38_exons_rev.bed
# Stdout: hg38_exons_rev_lastexons.bed
# Named input file: none
# Named output file: none
# Stdout is teed to removelastexons.py and onlylastexons.py

# USED BY: preRNAbinding.s
# USES: none

# EXAMPLE:

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
# # make lastexons.awk executable in preRNAbinding.s
# $ chmod +x lastexons.awk
# # take starting file hg38_exons.bed and reverse it for use in this awk program
# $ tac hg38_exons.bed > hg38_exons_rev.bed
# # find last exons in preRNAbinding.s
# $ cat hg38_exons_rev.bed | ./lastexons.awk > hg38_exons_rev_lastexons.bed

awk '
    BEGIN {count = 100}
{
	# If this exon contains exon_0
	if ($4 ~ /exon_0/)
	{
		# and if the previous exon contained exon_0 then count is 0 so print the line
		if (0 == count)
			{print $0};
		count = 0
	}
	else
		{count++};
	# If it is the first record print the line (the first record of the reversed file is the last exon of the original file)
	# If 1 == count then the previous line was exon_0 and this line was not exon_0 therefore this line is a last exon print the line
	if (1 == NR || 1 == count)
		{print $0}
}'

# AUTHOR:
# Thomas H. Lipscomb
# New York University 2018 MS thesis

# REPORTING BUGS:
# https://groups.google.com/forum/#!forum/bedtools-discuss
# post titled "How to make distance histograms of RNA-binding protein binding frequency near all types of pre-RNA junctions (intron-exon and exon-intron, etc.)"
# I may not be able to assist in debugging because I will no longer have access to the server I used to generate the example output directory "afterRunningBatchFile".

# COPYRIGHT: freeware copyleft