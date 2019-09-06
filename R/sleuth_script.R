# This script is for Sleuth
# Scott Schumacker
# Script for Scientific Data paper
# Schumacker et al.

# Download Bioconductor if you do not already have it installed
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()

# Download BiocUpgrade to update necessary Bioconductor packages if needed
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")     ## R version 2.15 or later

# Download BiomaRt
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")

# Install cowplot and load for plots and visuals
install.packages("cowplot")
library(cowplot)

install.packages("ggplot2")
library(ggplot2)

install.packages("ggrepel")
library(ggrepel)

# Download Sleuth
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")

# Load Sleuth into the project
library(sleuth)

# Modify the following line so the location of your kallisto results (each
# in their own folder) 
sample_id <- dir(file.path("/Users/ScottSchumacker/Desktop/SciDataPaper/MacRetData"))

# Return the names on your samples to check that they all populate
# in the console
sample_id

# Access each file path for the each sample. Again, paste the same
# file path as before and also write the sample_id
kal_dirs <- file.path("/Users/ScottSchumacker/Desktop/SciDataPaper/MacRetData",sample_id)

# Return the file path to the results containing the kallisto outputs
kal_dirs

# Create a table (.txt) file and save
# Load the table for the experimental design. Paste file path
# of the experimental design table in quotations
s2c <- read.table(file.path("/Users/ScottSchumacker/Desktop/SciDataPaper/MacRetTable.txt"),
                  header = TRUE,
                  stringsAsFactors = FALSE,
                  sep = "\t")

# Check the table to see if it appears in the console
s2c

# Add the directories as a final table column called path
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# check the table to see if it appears in the console
print(s2c)

# Load biomaRt into the project
library(biomaRt)

# List the Possible Marts to view and pick one of your choice
listMarts()

# Choose which mart to use. Insert mart of choice betwwen quotation marks
ensembl = useMart("ENSEMBL_MART_ENSEMBL")

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

# Renaming the output from BioMart to make it easier to work
# with in sleuth
t2g <- dplyr::mutate(t2g, target_id = paste(ensembl_transcript_id, transcript_version, sep = "."))
t2g <- dplyr::select(t2g, target_id, ens_gene = ensembl_gene_id, Gene_Name = external_gene_name, description, transcript_biotype)

# Check your listing of transcripts to see if they populate in 
# console. Not all transcripts will populate
head(t2g)

# Create the "Sleuth Object", a mathmatical
# data structure that contains our results
so <- sleuth_prep(s2c,
                  target_mapping = t2g,
                  extra_bootstrap_summary = TRUE,
                  read_bootstrap_tpm = TRUE)

# Check to see if target mapping worked with data and bioMArt was
# successful
View(so$target_mapping)

# Check the structure of the data with a PCA plot
# You will get some expected warnings

# Plot Basic PCA:
# Create custom bold plot elements
boldTextPlot <- element_text(face = "bold", color = "black", size = 19)
boldTextAxis <- element_text(face = "bold", color = "black", size = 19)

# Plot Basic PCA:
# With points
plot_pca(so, color_by = 'tissue', point_size = 7)

# With labels
PCA = plot_pca(so, color_by = 'tissue', text_labels = TRUE)

# Assigning the pca plot to the variable 'PCA'
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

# Create a heatmap for your samples
plot_sample_heatmap(so, use_filtered = TRUE, color_high = "white", color_low = "dodgerblue",
                    x_axis_angle = 75)

# Assign the heatmap to the variable 'Heat'
Heat <- plot_sample_heatmap(so, use_filtered = TRUE, color_high = "white", color_low = "dodgerblue",
                            x_axis_angle = 75)
dev.off()
# Return the plot check that it works
Heat

# Creating the plot edits (Bold Elements)
Heat <- Heat + theme(axis.text = boldTextAxis, text = boldTextPlot)
Heat <- Heat + aes(size = 19, face = "bold")

# Return 'Heat' to check final plot
Heat

# This Next Section is to Check for Possible outliers

# We can see genes involved in the the 1st PC by looking
# at the loadings (primary genes whose linear combinations define
# the principal components)

plot_loadings(so, pc_input = 2)

# Investigate if there are any outliers
# There may be a possible outliers
# Investigate further to see if they should be removed
plot_bootstrap(so, 'ENST00000361789.2', color_by = 'tissue')


#*******REMOVAL OF OUTLIERS*****************
# If there is an outlier, use the folling function to delete it

s2c <- dplyr::filter(s2c, sample != '')

# Set up the new Sleuth Object with the removed outlier.
# Use this command only if you deleted an outlier

so <- sleuth_prep(s2c,
                  target_mapping = t2g,
                  extra_bootstrap_summary = TRUE)

# After the new sleuth object is created, re-analyze the PCA plot to see
# if the data is distributed any better
plot_pca(so, color_by = 'condition', point_size = 7)
plot_pca(so, color_by = 'condition', text_labels = TRUE)
#*******END REMOVAL OF OUTLIERS**************


# Testing for differential genes

# Create the full model (i.e. a model that contains all covariates)
# then create an additional model with respect to one of the covariates
# finally compare both models to identify the differences

# In this scenario, the full model will contain just the covariate
# condition because all the sample under tissue are the same
so <- sleuth_fit(so, ~tissue, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
models(so)

# Create a likelihood ratio test
so <- sleuth_lrt(so, 'reduced', 'full')

# Create a Wald Test
so <- sleuth_wt(so, which_beta = 'tissueretina', which_model = 'full')

# Assign results
retMacSheet <- sleuth_results(so, 'tissueretina', 'wt')

the_results_lrt <- sleuth_results(so, 'reduced:full', 'lrt',
                                  show_all = FALSE)
# View Results
View(retMacSheet) # 188,753 rows of data before cleaning
retMacSheet_noNA <- na.omit(retMacSheet)
View(retMacSheet_noNA) #72,794 rows of data after cleaning


View(the_results_lrt)

# Filter out the significant genes according to a set qvalue
sleuth_significant <- dplyr::filter(the_results_wt, qval <= 0.05)
SigRetinaUpregulated <- dplyr::filter(sleuth_significant, b >= 0)
SigMaculaUpregulated <- dplyr::filter(sleuth_significant, b <= 0)

# Removal of NA values
sleuth_significant = na.omit(sleuth_significant)

# View results
View(sleuth_significant)
View(SigRetinaUpregulated)
View(SigMaculaUpregulated)

# Create plot of bootstraps representing differential transcript expression
plot_bootstrap(obj = so, target_id = 'ENST00000527717.5', 
               units = 'est_counts')

# View the first 20 genes in this list
# The number 20 is not a set number and can be changed
head(sleuth_significant, 20)

# Write the etire gene results table out into a 
# .csv file. Place the name of file into quotation marks
write.csv(retMacSheet_noNA,
          file = "totalData.csv")

# Use Rshiny to do some interactive visualizations
# Load Rshiny
sleuth_live(so)