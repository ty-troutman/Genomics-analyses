Processing RNA-seq data
================
Ty Troutman
May 2, 2022

### **This example is for stranded RNA-seq data, which is our lab’s [standard protocol](https://docs.google.com/document/d/10gv3azAfEpOEBIOYLgxlRpkMJf1_qg_o-lhHULrTp0g/edit?usp=sharing).**

#### **NOTE: These analyses can’t be run on the cluster head node. You will either need to generate a script to submit to the LSF scheduler, or request an interactive session.** **Quality control, mapping, and tag directories**

1.  Verify data checksums and rename following lab convention:
    \<filenameA\>\_R1.fastq.gz \<filenameA\>\_R2.fastq.gz
    \<filenameB\>\_R1.fastq.gz \<filenameB\>\_R2.fastq.gz

2.  Run [fastqc](https://github.com/s-andrews/FastQC) on the data.

3.  Group data according to mapping needs, e.g. by genome, paired end
    and single read, etc.

4.  Map the data using [STAR](https://github.com/alexdobin/STAR) and
    create [HOMER tag
    directories](http://homer.ucsd.edu/homer/ngs/tagDir.html). We have
    wrapper tools to automate and parallelize this process. There is one
    wrapper script for [single read](../Code/Wrapper_STAR_SingleRead.R)
    chemistry, and one for [paired
    end](../Code/Wrapper_STAR_PairedEnd.R) chemistry. **They depend on a
    loaded R module for use on the BMI HPC cluster.**

    ``` bash
     module load R/4.0.2
    ```

    #### **Caveats!! These tools: depend on file names in the lab convention (\*R1.fastq.gz); require input of only one fastq.gz file per sample per read direction; generate tag directories assuming a stranded library chemistry (e.g. dUTP RNA-seq protocol).**

    Execute the scripts with:

    ``` bash
    Rscript /data/hpc-troutman/software/scripts/Wrapper_STAR_SingleRead.R -d path/to/fastq/dataRscript /data/hpc-troutman/software/scripts/Wrapper_STAR_PairedEnd.R -d path/to/fastq/data
    ```

    Or with use of optional arguments:

    ``` bash
    Rscript /data/hpc-troutman/software/scripts/Wrapper_STAR_PairedEnd.R -d [fastq directory, character] -g [star genome index, character] -w [wall, hour:minute] -c [cores, integer] -m [memory, integer] -o [STAR output format, character]
    ```

    For example:

    ``` bash
    Rscript /data/hpc-troutman/software/scripts/Wrapper_STAR_PairedEnd.R \
    -d /data/troutman-nobackup/username/directory/ \
    -g mm10 -w 2:00 -c 24 -m 40000 -o BAM SortedByCoordinate
    ```

5.  Assess sequencing data using a genome browser.

    1.  The [UCSC Genome Browser](https://genome.ucsc.edu) is a good
        option. We still need to figure out how to host multiwig hubs
        for long term viability. Sharing UCSC browser hubs is a good way
        to share data though.

    2.  [IGV](https://software.broadinstitute.org/software/igv/) is
        another good option. It can be used locally.

#### **(Slightly more) advanced users.**

If RNA-seq data from the dUTP protocol is single read, then use the
-flip option in HOMER makeTagDirectory.

``` bash
makeTagDirectory <tagDirectoryName> <filename.sam> -genome <experiment_genome> \
-checkGC -format sam -flip
```

If RNA-seq data from the dUTP protocol is paired end then use the -sspe
-flip options in HOMER makeTagDirectory.

``` bash
makeTagDirectory <tagDirectoryName> <filename.sam> -genome <experiment_genome> \
-checkGC -sspe -format sam -flip
```

Data type, quality, etc. may require use of altered commands with STAR
and/or HOMER. More information on these applications can be found by
reading the documentation. For reference, below is an example bat script
for mapping paired end RNA-seq data with STAR and generating an
appropriately formatted HOMER tag directory.

``` bash
#BSUB -W 1:00
#BSUB -n 16
#BSUB -M 40000
#BSUB -e ./bsub_scripts/mouse_c57bl6j_Male_BMDM_RNA_KLA2hr_TDT_l20220224_ACGATCTA_GCAATATG_S14_L002.err
#BSUB -o ./bsub_scripts/mouse_c57bl6j_Male_BMDM_RNA_KLA2hr_TDT_l20220224_ACGATCTA_GCAATATG_S14_L002.out# load modulesmodule load STAR/2.7.9module load homer/4.11

# execute STAR
STAR \    
--runThreadN 16 \    
--genomeDir /data/hpc-troutman/data/genomes/indexes/star/mm10_starIndex/ \    --readFilesCommand zcat \    
--readFilesIn .//mouse_c57bl6j_Male_BMDM_RNA_KLA2hr_TDT_l20220224_ACGATCTA_GCAATATG_S14_L002_R1.fastq.gz .//mouse_c57bl6j_Male_BMDM_RNA_KLA2hr_TDT_l20220224_ACGATCTA_GCAATATG_S14_L002_R2.fastq.gz \    
--outFileNamePrefix ./star_out /mouse_c57bl6j_Male_BMDM_RNA_KLA2hr_TDT_l20220224_ACGATCTA_GCAATATG_S14_L002.mm10. \    --outSAMtype SAM

# remove STAR temporary directories
rm -r ./star_out/mouse_c57bl6j_Male_BMDM_RNA_KLA2hr_TDT_l20220224_ACGATCTA_GCAATATG_S14_L002.mm10._STARtmp

# execute HOMER makeTagDirectory
makeTagDirectory ./tag_directories/mouse_c57bl6j_Male_BMDM_RNA_KLA2hr_TDT_l20220224_ACGATCTA_GCAATATG_S14_L002.mm10/ ./star_out/mouse_c57bl6j_Male_BMDM_RNA_KLA2hr_TDT_l20220224_ACGATCTA_GCAATATG_S14_L002.mm10.Aligned*am -genome mm10 -checkGC -format sam -sspe -flip

# copy STAR log files to the HOMER tag directory
cp ./star_out/mouse_c57bl6j_Male_BMDM_RNA_KLA2hr_TDT_l20220224_ACGATCTA_GCAATATG_S14_L002.*Log* ./tag_directories/mouse_c57bl6j_Male_BMDM_RNA_KLA2hr_TDT_l20220224_ACGATCTA_GCAATATG_S14_L002.mm10/

echo Done!
```

This script can be named \<example\>.bat and submitted to the scheduler
on the cluster using:

``` bash
bsub < example.bat
```

Multiple bat scripts can be submitted using a bash for loop:

``` bash
for i in <*wildcardmatchingscriptname.bat>; do bsub < $i; sleep 1; done
```
