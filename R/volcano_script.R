# Script for creating volcano plot
# Create Cutom Volcano:
install.packages("dplyr")
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggrepel)

# Load in sleuth .csv data
inputData <- read.csv(file="/Users/ScottSchumacker/Desktop/SciDataPaper/ScientificData/totalData.csv", 
                 header=TRUE, sep=",")

View(inputData)
# mutating the data to tell which points are significant
inputData <- mutate(inputData, sig=ifelse(inputData$qval<0.05, "FDR < 0.05", "Not Significant"))

View(inputData)
# search for genes of interest and retrieve their index number
View(inputData$target_id)

# Assigning gene variable names
RHO = inputData$Gene_Name[4989]
PDE6A = inputData$Gene_Name[4041]
CNGA1 = inputData$Gene_Name[3667]
CNGB3 = inputData$Gene_Name[29401]
GNGT1 = inputData$Gene_Name[2211]
GUCA1A = inputData$Gene_Name[3131]
GUCA1B = inputData$Gene_Name[1845]
GUCY2D = inputData$Gene_Name[1824]
GUCY2F = inputData$Gene_Name[6613]
GNB1 = inputData$Gene_Name[1407]
SAG = inputData$Gene_Name[3522]
SLC24A1 = inputData$Gene_Name[2123]
NR2E3 = inputData$Gene_Name[3270]
CRX = inputData$Gene_Name[1064]
ENST00000226193.5 = inputData$target_id[1216]

ENST00000255622.10 = inputData$target_id[1143]

ENST00000560550.1 = inputData$target_id[566]

ENST00000467787.1 = inputData$target_id[2189]





# Creating the volcano plot variable
boldVolcTextPlot <- element_text(face = "bold", color = "black", size = 19)
boldVolcTextAxis <- element_text(face = "bold", color = "black", size = 13)


customPlot <- customPlot + theme(text = boldVolcTextPlot)
customPlot <- customPlot + theme(axis.text = boldVolcTextAxis)

customPlot = ggplot(inputData, aes(b, -log10(qval))) + geom_point(aes(x=b, y=-log10(qval), 
                                                        color = ifelse(-log10(qval)>-log10(0.05), "Statistically Significant", "no")), size=1.5)

customPlot <- customPlot + geom_point(aes(x=b, y = -log10(qval), color = ifelse(b < -1.5 & -log10(qval) > -log10(0.05) | b > 1.5 & -log10(qval) > -log10(0.05), "Biologically and Statistically Significant" , NA)))
customPlot <- customPlot + geom_point(aes(x=b, y = -log10(qval), color = ifelse(-log10(qval) < -log10(0.05) & (b) < -1.5 | (-log10(qval) < -log10(0.05) & (b) > 1.5), "Biologically Significant", NA)))
customPlot <- customPlot + scale_color_manual(values = c("#44AA99", "orange", "black","#332288"), label = c("no" = "Not Signficant"))
customPlot <- customPlot + theme(axis.text = boldVolcTextAxis)
customPlot <- customPlot + theme(axis.title = boldVolcTextAxis)

customPlot <- customPlot + geom_vline(xintercept=c(-1.5, 1.5), linetype="dashed")
customPlot <- customPlot + geom_hline(yintercept=c(1.30), linetype="dashed")
customPlot

# Labeling the gene names
customPlot <- customPlot + geom_label_repel(aes(x=b, y = -log10(qval), 
                              label = ifelse(
                                  
                                  inputData$target_id == ENST00000226193.5, "RCVRN", ifelse(
                                    
                                    inputData$target_id == ENST00000255622.10, "PDE6B", ifelse(
                                      
                                      inputData$Gene_Name == GUCY2F, "GUCY2F", ifelse(
                                        
                                        inputData$target_id == ENST00000467787.1, "GNAT1", ifelse(
                                          
                                        inputData$target_id == ENST00000560550.1, "NRL", ifelse(
                                          
                                                    inputData$Gene_Name == RHO, "RHO", NA))))))), 
                          color = "black", fontface = "bold", size=3, segment.size = 0.1, box.padding = unit(0.5, "lines"))


# Plot
customPlot
