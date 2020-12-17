# WorkshopQCBio_trial
Repository created for the QCBio trial workshps

# Workshop: part of NGS workshop: read mapping
# Outline
Expected outcomes
1. Today you will learn how to align your sequencing reads to a genome reference using BWA

Read Mapping exercise
--------------------------------

High-throughput sequencing technologies have in the past few years been 
producing millions of DNA and RNA sequences reads of human genome and
other species. To be useful, this genetic information has to be 'put
together' in a smart way, in the same way as the pieces of a puzzle
(reads) need to be mounted according to a picture (reference genome). 

![Applications](https://github.com/izabelcavassim/WorkshopQCBio_trial/blob/main/Images/applications.png)

Data
--------------------------------
In this exercise section you will be exposed to different softwares used
for mapping. We will use a dataset composed of 30
individuals from 3 different regions: Africa, EastAsia and WestEurasia.

```{bash}
head Sample_meta_data.txt
```

    ##       X...ID    ENA.RUN    population region country latitude longitude
    ## 1 ERS1042176 ERR1019075 Ju_hoan_North Africa Namibia    -18.9      21.5
    ## 2 ERS1042177 ERR1019076 Ju_hoan_North Africa Namibia    -18.9      21.5
    ## 3 ERS1042248 ERR1025622          Esan Africa Nigeria      6.5       6.0
    ## 4 ERS1042265 ERR1025639         Luhya Africa   Kenya      1.3      36.8
    ## 5 ERS1042266 ERR1025640      Mandenka Africa Senegal     12.0     -12.0
    ## 6 ERS1042267 ERR1025641      Mandenka Africa Senegal     12.0     -12.0
    ##      Sex       Illumina.ID
    ## 1   male LP6005441-DNA_B11
    ## 2   male LP6005441-DNA_A11
    ## 3 female LP6005442-DNA_B10
    ## 4   male LP6005442-DNA_E11
    ## 5   male LP6005441-DNA_E07
    ## 6 female LP6005441-DNA_F07

![simons_diversity](https://github.com/izabelcavassim/WorkshopQCBio_trial/blob/main/Images/unnamed-chunk-1-1.png)

This dataset is a subset of the Simons Diversity Project, and as you can
see, it covers a "bit" of the diversity of human population. If you want
to go further in details about this project and their findings, you can read of the publications 
[here](https://www.nature.com/articles/nature18964).

Software
--------------------------------
We will be using the software BWA(https://academic.oup.com/bioinformatics/article/25/14/1754/225615), 
BWA is standard software package for mapping low-divergent sequences (illumina reads) against a large reference genome––such as the human genome.
It contains different algorithms BWA, desgined for short reads (<= 100bp) and  BWA-SW and BWA-MEM for longer reads. Because we have high-quality queries (Simons diversity project) we will use BWA-MEM algorithm. This algorithm turns out to be the latest and fastest.

Log in to the server via terminal
---------------------------------

### For windows users
```{bash}
plink [USERNAME]@hoffman2.idre.ucla.edu
```
### For mac users
```{bash}
ssh [USERNAME]@hoffman2.idre.ucla.edu
```
Data source
-----------
The data is placed in a folder (in my directory) called **data** in the same
directory as users folder. 
In the following tutorial I am using one individual as an example
**ERR1019076**, but we could do the same analyses with every Simon's diversity individual.

Mapping reads against a reference genome
-----------------------------------

We will be using the bwa mapper. If you are interested in understanding
a bit more of the software and its algorithm, you can look it up
[here](http://bio-bwa.sourceforge.net/bwa.shtml). We have thousands of
reads and we want to find out their best location in the genome. We
decided to focus on a 10 MB region of chromosome 2, which can be
downloaded through
[Ensembl](ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/).
This region goes from 135MB to 145MB and it is known to containg the
lactase gene.

Two input files are needed to do genome mapping:

-   Fasta file containing your reference genome
    ([GRCh37](http://grch37.ensembl.org/index.html))
-   The reads in fastq format.

As you learned during the theoretical lecture, fastaq format is a text format that stores 
both the biological sequence and its related quality score. To refresh your memory the explanation for the format is found [here](https://en.wikipedia.org/wiki/FASTQ_format)

We first need to index the reference file for later use. This step is
important for the speed and process of the mapping algorithm. It takes
around 4 minutes. This creates a collecion of files that are used by BWA
to perform the alignment.

Create a soft-link of fasta reference to your home folder (type pwd, to know where you are):

    ln -s /u/home/m/mica20/QCBioWorkshoptrial/data/Homo_sapiens.GRCh37.75.dna.chromosome.2.fa /home/user_name/

Then produce the indexes needed for bwa mapper:

    bwa index -a bwtsw Homo_sapiens.GRCh37.dna_rm.chromosome.2.fa.gz
    
It takes some seconds, so we can start discussing about the next steps. 
Once we have the index file fro the reference sequence created then we can align our reads to that index
 
Multilple options our found including:
* -t INT        number of threads [1]
* -k INT        minimum seed length [19]
* -w INT        band width for banded alignment [100] (gaps longer than 100 will not be found)
* -p            Assume the first input query file is interleaved paired-end FASTA/Q
* -a            output all alignments for SE or unpaired PE

Now you can map the reads back to the reference. This will take around
10 minutes. You can start installing the software that will be used
later in this tutorial (IGV) while you wait for it.

    bwa mem -t 16 -p Homo_sapiens.GRCh37.dna_rm.chromosome.2.fa.gz ERR1019076_reads_135_145.fq > mapped_ERR1019076_reads_135_145.sam

While it is running, let's remember out how the sam format looks like:

![SAM format](https://www.samformat.info/images/sam_format_annotated_example.5108a0cd.jpg)




Have a look at the bam file generated:

    samtools view mapped_ERR1019076.sam | head

Get some useful stats of your mapping:

    samtools flagstat mapped_ERR1019076.sam

Once the map is generated, you can index the bam file to visualize it
using the software IGV. Indexing a genome sorted BAM file allows one to quickly
extract alignments overlapping particular genomic regions. Moreover,
indexing is required by genome viewers such as IGV so that the viewers
can quickly display alignments in each genomic region to which you
navigate.

    samtools index ERR1019076.bam

Dowloading via terminal
-----------------------

You can download the data via terminal by the following:

    scp -P 8922 root@185.45.23.197:/home/Data/ERR1019076.bam Directory/My_computer

IGV software
------------

IGV is an Integrative Genomics viewer and can be very useful to look at
the results of Mapping and SNP calling. We have not installed it in the
cluster, so you can download it to your machine you can go to its
[website](http://software.broadinstitute.org/software/igv/). Three files
are necessary to look at this dataset: a reference sequence and the
**.bam** and **.bai** files, download it from the cluster in a specific
directory. Since we are using a human reference, the sequence is already
available in the software:

Go to Genomes ----&gt; Load Genome from server... ----&gt; Filter by
human and choose the Human hg19 reference (which is the GRCh37).

After it you will the chromosomes and genes. Now you can download the
mapping results by typing: File ----&gt; Load from File... ----&gt;
ERR1019076.bam.

You will see something like this: ![IGV bronwser
illustration](https://github.com/izabelcavassim/PG_2018/blob/master/Exercises/Week02/IGV_example.png)

Try to understand what are the different attributes present in the
viewer. If you zoom in very much you will find single nucleotide
polymorphisms (SNPs), where the reference sequence does not have the
same nucleotide as the data mapped to.

Plotting results
----------------

One of the attributes one could learn from mapping reads back to the
reference is the coverage of reads across the genome. In order to
calculate the coverage depth you can use the command **samtools depth**.

    samtools depth ERR1019076.bam > deduped_ERR1019076.coverage

You can have a look at the resulted file. What do you find in the three
different columns?

    less deduped_ERR1019076.coverage

Now open the subset file in R and plot it. You can do it in the
terminal, you just need to type R.

    library(ggplot2)
    library(dplyr)

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    scaf <- read.table("/Users/PM/Desktop/PHD_incomplete/data/deduped_ERR1019076.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE, col.names = c("Scaffold", "locus", "depth"))
      
    head(scaf)

    ##   Scaffold  locus depth
    ## 1        2 833855     1
    ## 2        2 833856     1
    ## 3        2 833857     1
    ## 4        2 833858     1
    ## 5        2 833859     1
    ## 6        2 833860     1

    # Compressing the dataframe in windows
    scaf %>% 
    mutate(rounded_position = round(locus, -2)) %>%
        group_by(rounded_position) %>% 
            summarize(mean_cov = mean(depth)) -> compressed

    # Plotting the data
    p <- ggplot(data =  compressed, aes(x=rounded_position, y=mean_cov)) + geom_area() + theme_classic() + ylim(0, 400)

    #p

    # Saving your coverage plot
    ggsave(p, device = "pdf")

    ## Saving 7 x 5 in image

What are the conclusions you can extract from these analysis? Does the
coverage match with what you observed with IGV?
