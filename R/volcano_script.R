# Script for creating volcano plot
# Create Cutom Volcano:
install.packages("dplyr")
library(dplyr)

# Load in sleuth .csv data
inputData <- read.csv(file="/Users/ScottSchumacker/Desktop/SciDataPaper/ScientificData/totalData.csv", 
                 header=TRUE, sep=",")

# mutating the data to tell which points are significant
inputData <- mutate(inputData, sig=ifelse(inputData$qval<0.05, "FDR < 0.05", "Not Significant"))

View(inputData)
# search for genes of interest and retrieve their index number
View(inputData$Gene_Name)

# Assigning gene variable names
RHO = inputData$Gene_Name[4989]
PDE6B = inputData$Gene_Name[1143]
PDE6A = inputData$Gene_Name[4041]
CNGA1 = inputData$Gene_Name[3667]
CNGB3 = inputData$Gene_Name[29401]
GNAT1 = inputData$Gene_Name[2189]
GNGT1 = inputData$Gene_Name[2211]
GUCA1A = inputData$Gene_Name[3131]
GUCA1B = inputData$Gene_Name[1845]
GUCY2D = inputData$Gene_Name[1824]
GUCY2F = inputData$Gene_Name[6613]
RCVRN = inputData$Gene_Name[1216]
GNB1 = inputData$Gene_Name[1407]
SAG = inputData$Gene_Name[3522]
SLC24A1 = inputData$Gene_Name[2123]
NRL = inputData$Gene_Name[566]
NR2E3 = inputData$Gene_Name[3270]
CRX = inputData$Gene_Name[1064]

# Creating the volcano plot variable
customPlot = ggplot(inputData, aes(b, -log10(qval))) + geom_point(aes(x=b, y=-log10(qval), 
                                                        color = ifelse(-log10(qval)>-log10(0.05), "Statistically Significant", "no")), size=1.5)

customPlot <- customPlot + geom_point(aes(x=b, y = -log10(qval), color = ifelse(b < -1.5 & -log10(qval) > -log10(0.05) | b > 1.5 & -log10(qval) > -log10(0.05), "Biologically and Statistically Significant" , NA)))
customPlot <- customPlot + geom_point(aes(x=b, y = -log10(qval), color = ifelse(-log10(qval) < -log10(0.05) & (b) <= -1.5 | (-log10(qval) < -log10(0.05) & (b) > 1.5), "Biologically Significant", NA)))
customPlot <- customPlot + scale_color_manual(values = c("green", "orange", "black","red"), label = c("no" = "Not Signficant"))
customPlot

customPlot <- customPlot + geom_vline(xintercept=c(-1.5, 1.5), linetype="dashed")
customPlot <- customPlot + geom_hline(yintercept=c(1.30), linetype="dashed")
customPlot

# Labeling the gene names
customPlot <- customPlot + geom_label_repel(aes(x=b, y = -log10(qval), 
                              label = ifelse(
                                  
                                  inputData$Gene_Name == GUCA1A, "GUCA1A", ifelse(
                                    
                                    inputData$Gene_Name == GUCA1B, "GUCA1B", ifelse(
                                      
                                      inputData$Gene_Name == GUCY2F, "GUCY2F", ifelse(
                                                    
                                                    inputData$Gene_Name == RHO, "RHO", NA))))), 
                          color = "black", fontface = "bold", size=3, segment.size = 0.1, box.padding = unit(0.5, "lines"))
# Plot
customPlot
