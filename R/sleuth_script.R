# R Script for analyzing transcript abundances calculated by Kallisto

# Download Bioconductor if you do not already have it installed
source("https://bioconductor.org/biocLite.R")
biocLite()

# Download BiocUpgrade to update necessary Bioconductor packages if needed
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")

# Download BiomaRt and load biomaRt
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)

# Install and load the following packages for data visualization
install.packages("cowplot")
library(cowplot)

install.packages("ggplot2")
library(ggplot2)

install.packages("ggrepel")
library(ggrepel)

install.packages("ggrepel")
library(ggrepel)

# Install and load Sleuth
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")
library(sleuth)

# Paste the data path to your Kallisto outputs directory
sample_id <- dir(file.path(""))

# Return the names on your samples to check that they all populate
sample_id

# Paste the same data path your Kallisto outputs directory
kal_dirs <- file.path("",sample_id)

# Return the file path to the results containing the kallisto outputs
kal_dirs

# Paste the file path to your design matrix table file
s2c <- read.table(file.path(""),
                  header = TRUE,
                  stringsAsFactors = FALSE,
                  sep = "\t")

# Return the matrix to see if it populates correctly
s2c

# Add the sample directories as a final table column called path
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# Return the table to see if it populates in the console correctly
print(s2c)

# List the Possible marts to view and pick one of your choice
listMarts()

# Choose which mart to use. Insert mart of choice betwwen quotation marks
# For this project, "ENSEMBL_MART_ENSEMBL
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")

# List datasets within mart of choice to see which one you want/need
listDatasets(ensembl)

# Load genes of choice by defining which mart, dataset, and the host our 
# data
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = "www.ensembl.org")

# Get target transcripts from bioMart
t2g <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "transcript_version",
                 "ensembl_gene_id", "external_gene_name", "description",
                 "transcript_biotype"), mart = mart)

# Rename the output from BioMart to make it easier to work with in Sleuth
t2g <- dplyr::mutate(t2g, target_id = paste(ensembl_transcript_id, transcript_version, sep = "."))
t2g <- dplyr::select(t2g, target_id, ens_gene = ensembl_gene_id, 
                     Gene_Name = external_gene_name, description, 
                     transcript_biotype)

# Return a preview of your listing of transcripts to see if they populate
head(t2g)

# Create the "Sleuth Object", a mathmatical
# data structure that contains our results
so <- sleuth_prep(s2c,
                  target_mapping = t2g,
                  extra_bootstrap_summary = TRUE,
                  read_bootstrap_tpm = TRUE)

# Check to see if target mapping worked with data and bioMArt was successfull
View(so$target_mapping)

# Check the structure of the data with a PCA plot

# Plot Basic PCA:
# Create custom bold plot elements
boldTextPlot <- element_text(face = "bold", color = "black", size = 19)
boldTextAxis <- element_text(face = "bold", color = "black", size = 19)

# Plot basic PCA with points:
# "tissue" may be different with your experiment.
# 'tissue" is equivalent to the condition column in the design matrix
plot_pca(so, color_by = 'tissue', point_size = 7)

# Plot with labels
PCA <- plot_pca(so, color_by = 'tissue', text_labels = TRUE)

# Assign the PCA plot to the variable 'PCA'
PCA <- plot_pca(so, color_by = 'tissue', point_size = 7)

# Return the plot to view
PCA

# Creating plot edits (Bold Elements)
PCA <- PCA + theme(text = boldTextPlot)
PCA <- PCA + theme(axis.text = boldTextAxis)

# Repeling the data labels with ggrepel
PCA <- PCA + 
  geom_label_repel(aes(label = sample_id, size = 19),
                   box.padding   = 0.7, 
                   point.padding = 0.7,
                   segment.color = 'black', label.size = 0.40)

# Return PCA to check final plot
PCA

# Test for differential gene expression

# Create two models to compare for a difference
so <- sleuth_fit(so, ~tissue, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
models(so)

# Use the following function to perform Likelihood Ratio Test
# so <- sleuth_lrt(so, 'reduced', 'full')

# A Wald test was used for this specific project
so <- sleuth_wt(so, which_beta = 'tissueretina', which_model = 'full')

# Assign results
results <- sleuth_results(so, 'tissueretina', 'wt')

# the_results_lrt <- sleuth_results(so, 'reduced:full', 'lrt',
#                                  show_all = FALSE)

# View Results for the Wald test
View(results)

# Delete any data that has value NA
results_noNA <- na.omit(results)
View(results_noNA)

# View the results of lrt
# View(the_results_lrt)

# Filter out the significant genes according to a set qvalue
sleuth_significant <- dplyr::filter(results_noNA, qval <= 0.05)

# Filter for up and down regulation based on "b" value
SigRetinaUpregulated <- dplyr::filter(sleuth_significant, b >= 0)
SigMaculaUpregulated <- dplyr::filter(sleuth_significant, b <= 0)

# Removal of NA values
sleuth_significant = na.omit(sleuth_significant)

# View results
View(sleuth_significant)
View(SigRetinaUpregulated)
View(SigMaculaUpregulated)

# Return a preview of the first 20 genes in this list
# The number 20 is not a set number and can be changed
head(sleuth_significant, 20)

# Write the etire gene results table out into a 
# .csv file. Place the name of file into quotation marks
write.csv(results,
          file = "Results.csv")

# Use Rshiny to do some interactive visualizations
sleuth_live(so)