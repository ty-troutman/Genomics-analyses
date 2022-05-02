######################################

# Author: Ty Dale Troutman, Ph.D.
# Contact: Ty.Troutman@cchmc.org
# Created October, 21 2018

#!/usr/bin/env Rscript

######################################

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4) {
  stop("
This script takes a HOMER RNAseq output file and removes genes of less than the specified length and less than the specified expression in at least one data column.

Usage of the script is as follows, including three arguments:
Rscript /home/ttroutman/scripts/filterRNAseq.R <getDiffExpression.pl output file> <Minimum tags in one data set> <Desired output path/name.txt>
       
For example: 
Rscript /home/ttroutman/scripts/filterRNAseq.R /directory/file.txt 200 2 /directory/filter_diff_noadj.txt"
       , call.=FALSE)
}

# Read in arguments from command line
diffExpressionPath <- args[1]
length <- as.double(args[2])
tags <- as.double(args[3])
outputPath <- args[4]

# Load Libraries
suppressMessages(library("dplyr"))

# Read in the data files
inputData <- read.delim(diffExpressionPath, sep="\t", row.names = 1, header = TRUE, skip = 0)

# Remove genes with less than 'length' criteria
inputData <- filter(inputData,inputData$Length>length)

# Remove genes with less than 'expression' criteria
filterData <- inputData[apply(select(inputData, matches("Tag.Count")), MARGIN=1, function(x) any(x > tags)),]
filterData <- inputData[apply(select(inputData, matches("TPM")), MARGIN=1, function(x) any(x > tags)),]

# Send filtered data to the output file name
write.table(filterData, outputPath,  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Finished!
