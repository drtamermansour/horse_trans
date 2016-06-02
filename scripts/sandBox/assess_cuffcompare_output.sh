## Cuffcompre output description of reference annptation
# Kept 1820 ref transcripts out of 1938
# 118 duplicate reference transcripts discarded.
# Reference mRNAs :    1820 in    1786 loci  (928 multi-exon)

## count the super-loci
cat cuffcmp.loci | awk -F "\t" '{ print $3 }' | sort | uniq | wc -l | awk '{print $1-1}'    #1613


#####
## no of transcripts
cat cuffcmp.tracking | awk -F "\t" '{ print $1 }' | wc -l    #351700

## Complete match of intron chain
cat cuffcmp.tracking | awk -F "\t" ' $4 == "=" ' | wc -l    #1606

## Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript
cat cuffcmp.tracking | awk -F "\t" ' $4 == "j" ' | wc -l    #4463

## Contained
cat cuffcmp.tracking | awk -F "\t" ' $4 == "c" ' | wc -l    #0

## Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment.
cat cuffcmp.tracking | awk -F "\t" ' $4 == "e" ' | wc -l    #53

## A transfrag falling entirely within a reference intron
cat cuffcmp.tracking | awk -F "\t" ' $4 == "i" ' | wc -l    #2885

## Generic exonic overlap with a reference transcript
cat cuffcmp.tracking | awk -F "\t" ' $4 == "o" ' | wc -l    #3170

## Possible polymerase run-on fragment (within 2Kbases of a reference transcript)
cat cuffcmp.tracking | awk -F "\t" ' $4 == "p" ' | wc -l    #541

## Unknown, intergenic transcript
cat cuffcmp.tracking | awk -F "\t" ' $4 == "u" ' | wc -l    #264552

## Exonic overlap with reference on the opposite strand
cat cuffcmp.tracking | awk -F "\t" ' $4 == "x" ' | wc -l    #1020

## An intron of the transfrag overlaps a reference intron on the opposite strand (likely due to read mapping errors)
cat cuffcmp.tracking | awk -F "\t" ' $4 == "s" ' | wc -l    #1

## (.tracking file only, indicates multiple classifications)
cat cuffcmp.tracking | awk -F "\t" ' $4 == "." ' | wc -l    #73409




