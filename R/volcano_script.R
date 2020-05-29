# R script for creating volcano plot from Kallisto/Sleuth outputs

# Install and load data visualization packages
install.packages("dplyr")
library(dplyr)

install.packages("ggplot2")
library(ggplot2)

install.packages("cowplot")
library(cowplot)

install.packages("ggrepel")
library(ggrepel)

# Load in sleuth .csv data results
inputData <- read.csv(file="/Users/ScottSchumacker/Desktop/SciDataPaper/ScientificData/totalData.csv", 
                 header=TRUE, sep=",")

# Return the results
View(inputData)

# Mutate the data to identify which points are significant with q-val < 0.05
inputData <- mutate(inputData, sig=ifelse(inputData$qval<0.05, "FDR < 0.05", "Not Significant"))

# Return the result
View(inputData)

# Search for transcripts of interest and retrieve their index number
View(inputData$target_id)

# For gene labels, 'ENST' is transcript_id name and 
# target_id[] is it's corresponding index number

# Genes used in this example:
# NRL, RCVRN, PDE6B, RHO, GUCY2F, GNAT1

# RCVRN
ENST00000226193.5 <- inputData$target_id[] # Enter index of transcript of choice

# PDE6B
ENST00000255622.10 <- inputData$target_id[]

# NRL
ENST00000560550.1 <- inputData$target_id[]

# GNAT
ENST00000467787.1 <- inputData$target_id[]

# RHO
ENST00000296271.3 <- inputData$target_id[]

# Create bold text elements
boldVolcTextPlot <- element_text(face = "bold", color = "black", size = 19)
boldVolcTextAxis <- element_text(face = "bold", color = "black", size = 13)

# Addeing bold text theme
customPlot <- customPlot + theme(text = boldVolcTextPlot)
customPlot <- customPlot + theme(axis.text = boldVolcTextAxis)

# Creating custom plot
customPlot <- ggplot(inputData, aes(b, -log10(qval))) + geom_point(aes(x=b, y=-log10(qval), 
                                                        color = ifelse(-log10(qval)>-log10(0.05), "Statistically Significant", "no")), size=1.5)

customPlot <- customPlot + geom_point(aes(x=b, y = -log10(qval), color = ifelse(b < -1.5 & -log10(qval) > -log10(0.05) | b > 1.5 & -log10(qval) > -log10(0.05), "Biologically and Statistically Significant" , NA)))
customPlot <- customPlot + geom_point(aes(x=b, y = -log10(qval), color = ifelse(-log10(qval) < -log10(0.05) & (b) < -1.5 | (-log10(qval) < -log10(0.05) & (b) > 1.5), "Biologically Significant", NA)))
customPlot <- customPlot + scale_color_manual(values = c("#44AA99", "orange", "black","#332288"), label = c("no" = "Not Signficant"))
customPlot <- customPlot + theme(axis.text = boldVolcTextAxis)
customPlot <- customPlot + theme(axis.title = boldVolcTextAxis)

customPlot <- customPlot + geom_vline(xintercept=c(-1.5, 1.5), linetype="dashed")
customPlot <- customPlot + geom_hline(yintercept=c(1.30), linetype="dashed")
customPlot

# Label the gene/transcript names
customPlot <- customPlot + geom_label_repel(aes(x=b, y = -log10(qval), 
                              label = ifelse(
                                  
                                  inputData$target_id == ENST00000226193.5, "RCVRN", ifelse(
                                    
                                    inputData$target_id == ENST00000255622.10, "PDE6B", ifelse(
                                      
                                      inputData$Gene_Name == GUCY2F, "GUCY2F", ifelse(
                                        
                                        inputData$target_id == ENST00000467787.1, "GNAT1", ifelse(
                                          
                                        inputData$target_id == ENST00000560550.1, "NRL", ifelse(
                                          
                                                    inputData$target_id == ENST00000296271.3, "RHO", NA))))))), 
                          color = "black", fontface = "bold", size=3, segment.size = 0.1, box.padding = unit(0.5, "lines"))


# Return the custom volcano plot
customPlot