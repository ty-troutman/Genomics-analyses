# Author: Ty.Troutman@cchmc.org
# Date: April 26, 2022

################################## Description #################################

# This program automates bsub script assembly for SINGLE-READ STAR mapping.
# Data should be fastq.gz, in the same folder, with single fastq.gz files per
# sample, and named in the following convention:

#                 <sampleA>_R1.fastq.gz

################################################################################

# quietly load optparse
suppressPackageStartupMessages(library("optparse"))

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(
    c("-d", "--inputDirectory"),
    type = "character",
    help = "Input directory containing fastq.gz files. [default = %default]",
    metavar = "character"
  ),
  make_option(
    c("-g", "--genome"),
    type = "character",
    default = "mm10",
    help = "Genome assembly to use. Options include mm10 or hg38. [default = %default]",
    metavar = "character"
  ),
  make_option(
    c("-w", "--wall"),
    type = "character",
    default = "2:00",
    help = "Amount of wall time requested. Larger data may require increasing the wall time. [default = %default]",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cores"),
    type = "numeric",
    default = "16",
    help = "Number of cores requested. [default = %default]",
    metavar = "character"
  ),
  make_option(
    c("-m", "--memory"),
    type = "numeric",
    default = "40000",
    help = "Amount of ram requested. Larger data may require more RAM. [default = %default]",
    metavar = "character"
  ),
  make_option(
    c("-o", "--outSAMtype"),
    type = "character",
    default = "BAM SortedByCoordinate",
    help = "Output format for mapped data using STAR command requirements. This can be SAM or BAM SortedByCoordinate. [default = %default]",
    metavar = "character"
  )
)

opt_parser = OptionParser(
  "\n\tRscript Wrapper_STAR_SingleRead.R -d path/to/fastq/data [options]",
  option_list = option_list,
  prog = "Wrapper_STAR_SingleRead.R", 
  epilogue = "\tThis script takes an input directory of fastq.gz files and generates .bat files.
  \tThe .bat file provides instructions for mapping fastq.gz data with STAR/2.7.9 and creation of homer/4.11 tag directories.
  \n\tCommands are appropriate for stranded RNA-seq data (resulting from the standard Troutman lab protocol) and only single-reads [*R1.fastq.gz]! 
  \n\tThe .bat files are then submitted to the LSF scheduler on the Biomedical Informatics Computational Cluster.
  \tUsage of the script should be done in your personal or lab scratch folder.
  \n\tOptional arguments can be included as follows:
  \t\tRscript /data/hpc-troutman/software/scripts/Wrapper_STAR_SingleRead.R -d [fastq directory, character] -g [star genome index, character] -w [wall, hour:minute] -c [cores, integer] -m [memory, integer] -o [STAR output format, character]
  \n\tFor example:
  \t\tRscript /data/hpc-troutman/software/scripts/Wrapper_STAR_SingleRead.R -d /data/troutman-nobackup/username/directory/ -g mm10 -w 2:00 -c 24 -m 40000 -o BAM SortedByCoordinate\n\n"
)

if (length(args)<1){
  print_help(opt_parser)
  stop("path/to/fastq/data/ is required.\n\n\n", call.=FALSE)
}

opt = parse_args(opt_parser,args)

# declare variables
dirPath <- opt$inputDirectory
genome <- opt$genome
genomeDir <- paste0("/data/hpc-troutman/data/genomes/indexes/star/", opt$genome, "_starIndex/")
wall <- opt$wall
cores <- opt$cores
memory <- opt$memory
outSAMtype <- opt$outSAMtype

# quietly load tidyverse
suppressWarnings(suppressMessages(require(tidyverse)))

# create directories for bsub scripts and mapping output
if (dir.exists("bsub_scripts")) {
  
  print("bsub_scripts directory exists")
} else {
  
  # create the "my_new_folder
  dir.create("bsub_scripts")
}

if (dir.exists("star_out")) {
  
  print("out directory exists")
} else {
  
  # create the "my_new_folder
  dir.create("star_out")
}

if (dir.exists("tag_directories")) {
  
  print("out directory exists")
} else {
  
  # create the "my_new_folder
  dir.create("tag_directories")
}

# read in and wrangle fastq.gz file names
# NOTE: wrangling steps are convention specific!
fastqData <-
  tibble(
    read1 = list.files(file.path(dirPath), pattern = "_R1.fastq.gz"),
  ) %>%
  mutate(outPrefix = str_replace(read1, "_R1.fastq.gz", "")) %>%
  select(outPrefix, read1)

# function to print .bat script
fn1 <- function(read1, outPrefix){
  print(str_c(
    "#BSUB -W ", wall,
    "\n#BSUB -n ", cores,
    "\n#BSUB -M ", memory,
    "\n#BSUB -e ./bsub_scripts/", outPrefix, ".err",
    "\n#BSUB -o ./bsub_scripts/", outPrefix, ".out\n",
    "\n# load modules \n",
    "module load STAR/2.7.9\n",
    "module load homer/4.11\n",
    "\n# execute STAR\n",
    
    "STAR \\
    --runThreadN ", cores, " \\
    --genomeDir ", genomeDir, " \\
    --readFilesCommand zcat \\
    --readFilesIn ", file.path(dirPath, read1), " \\
    --outFileNamePrefix ./star_out/", outPrefix, ".", genome, ". \\
    --outSAMtype ", outSAMtype, "\n",
    
    "\n# remove STAR temporary directories\n",
    "rm -r ./star_out/", outPrefix, ".", genome, "._STARtmp\n",
    
    "\n# execute HOMER makeTagDirectory\n",
    "makeTagDirectory ./tag_directories/", outPrefix, ".", genome, "/ \\
    ./star_out/", outPrefix, ".", genome, ".Aligned*am \\
    -genome ", genome, " -checkGC -format sam -flip\n",
    
    "\n# copy STAR log files to the HOMER tag directory\n",
    "cp ./star_out/", outPrefix, ".", "*Log* ./tag_directories/", outPrefix, ".", genome, "/ \\",
    
    "\n\n echo Done!"
  ))
}

# function to write out the printed .bat script
output_bat <- function(data, names){
  write_lines(data, paste0("bsub_scripts/",names, ".bat"))
}

# execute functions:
# build scripts; name list objects; name output files; write scripts to disk
fastqData %>%
    pmap(.f = fn1) %>%
    set_names(., fastqData$outPrefix) %>% 
    list(data = ., names = fastqData$outPrefix) %>% 
    pmap(output_bat) %>% invisible

# submit .bat scripts to the LSF scheduler
system(paste("for i in bsub_scripts/*.bat; do bsub < $i; sleep 0.3; done"))
