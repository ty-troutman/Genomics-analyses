######################################

# Author: Ty Dale Troutman, Ph.D.
# Contact: Ty.Troutman@cchmc.org
# Created October, 21 2018

#!/usr/bin/env Rscript

######################################

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("
This script takes a HOMER RNAseq output file and removes genes of less than the specified length.

Usage of the script is as follows, including three arguments:
Rscript /home/ttroutman/scripts/filterRNAseq.R <getDiffExpression.pl output file> <Minimum length> <Desired output path/name.txt>
       
For example: 
Rscript /home/ttroutman/scripts/filterRNAseq.R /directory/file.txt 200 /directory/filter_diff_noadj.txt"
       , call.=FALSE)
}

# Read in arguments from command line
inputPath <- args[1]
length <- as.double(args[2])
outputPath <- args[3]

# Read in the data files
inputData <- read.delim(inputPath, as.is=T, header = T, sep = "\t", check.names = F, quote="", comment.char="")

# Remove genes with less than 'length' criteria
inputData <- inputData[inputData$Length>length,]

# Send filtered data to the output file name
write.table(inputData, outputPath,  sep = "\t", row.names = F, quote = FALSE, col.names = TRUE)

# Finished!
