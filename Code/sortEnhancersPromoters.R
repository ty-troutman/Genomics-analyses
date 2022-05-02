#!/usr/bin/env Rscript

######################################

# Author: Ty Dale Troutman, Ph.D.
# Contact: Ty.Troutman@cchmc.org
# Created October 2, 2019

######################################

args <- commandArgs(trailingOnly = T)

if (length(args) < 2) {
  stop("
       This script takes a HOMER annotatePeaks.pl input file and outputs two files of specified distance from the TSS
       
       Usage of the script is as follows, including three arguments:
       Rscript /home/ttroutman/scripts/subsetEnhancersPromoters.R <annotatePeaks.pl output file> <Distance from TSS>
       
       For example: 
       Rscript /home/ttroutman/scripts/subsetEnhancersPromoters.R /directory/file.txt 3000 "
       , call.=FALSE)
}

# Read in arguments from command line
inputPath <- args[1]
length <- as.double(args[2])

# Read in the data files
inputData <- read.delim(inputPath,as.is=T,header=T,sep="\t",check.names=F,quote="",comment.char="",row.names=1)
# inputData <- read.table(inputPath,as.is=T,sep="\t",header=T,row.names=1,quote="",comment.char="")

# Remove genes with less than 'length' criteria
enhancer <- inputData[abs(inputData$`Distance to TSS`) >=length,]
promoter <- inputData[abs(inputData$`Distance to TSS`) <length,]

# Send filtered data to the output file name
write.table(enhancer, paste("enhancer_",length,"_",inputPath,sep = ""),  sep = "\t", row.names = T, quote = FALSE, col.names = TRUE, na = "")
write.table(promoter, paste("promoter_",length,"_",inputPath,sep = ""),  sep = "\t", row.names = T, quote = FALSE, col.names = TRUE, na = "")

# Finished!