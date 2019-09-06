# Schumacker et al., 2019 Scientific Data manuscript

This repository contains code used for quality control and data analysis presented in: 

> **Schumacker ST, Coppage KR, Enke RA. RNA Sequencing Analysis of the Human Retina and Associated Ocular Tissues. Scientific Data 2019 (in prep).**

----

## Data availability

Data is available in NCBI SRA under the accession number (https:).

## Software used

| Software | Version | URL | 
| --- | --- | --- |
| Trimmomatic | 0.36 | http://www.usadellab.org/cms/?page=trimmomatic  |
| FastQC | 0.11.5 | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| kallisto | 0.42.3 | https://pachterlab.github.io/kallisto/ |
| MultiQC | 1.7 | https://multiqc.info|
| sleuth | - | https://pachterlab.github.io/sleuth/ |

## Data analysis walkthroughs & code

Walkthroughs and code used for all of the quality assessment and data analysis steps are available in each of the below links.

1. [Quality assessment with FastQC](https://github.com/ScottSchumacker/Schumacker2019_Sci_Data/blob/master/Walkthroughs/FastQC-Instructions)
2. [Sequence trimming with Trimmomatic]()
3. [Quality Analysis summary with MultiQC]()
3. [Quantitation with kallisto](https://github.com/ScottSchumacker/Schumacker2019_Sci_Data/blob/master/Walkthroughs/Kallisto-Instructions)
4. [Normalization, visualization, and differential expression analysis with sleuth]()
