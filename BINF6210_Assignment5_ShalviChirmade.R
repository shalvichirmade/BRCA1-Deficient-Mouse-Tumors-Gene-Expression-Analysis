#### BINF 6210 - Assignment 5 - Due Friday December 17, 2021 by 5 pm ----
# By Shalvi Chirmade


### TITLE

#### 1- INTRODUCTION ----








#### 2- DATA SET ----









#### 3- DATA ACQUISITION, EXPLORATION, FILTERING AND QUALITY CONTROL ----

##If needed, install the packages below by removing the comment on the specific line. Otherwise, please load in these required packages.

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(c("limma", "Glimma", "edgeR", "Mus.musculus"))
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

#install.packages("R.utils")
library(R.utils)
#install.packages("gplots")
library(gplots)
#install.packages("tidyverse")
library(tidyverse)


##Reading in the data from the GEO database. As I am working on a Mac, it automatically unzips the file hence the file used here is only .txt. This file is provided to you.

#Load the complete data set into the script.
dfComplete_Data <- read.delim("GSE148569_scRNAseq_C1_count.txt")

#Inspect the data.
head(complete_data, n = c(10, 5))
#There are 477 samples in this data set and 21004 genes analyzed.

#This data set has an association to a published paper by Sun et al. This paper is not explicit with the details of this data file so I found it quite difficult to come up with a conclusion for what the data is showing me. I speculate that each of the different SC samples correspond to a different mouse; so SC1, SC2 and so on depict the names of the particular mice. Hence, each sample associated with a particular SC is a sample from that mouse. For example, if we take the first sample from the data frame we visualized above, SC1_01 refers to the 01 sample from mouse SC1. So, sample SC5_29 would mean sample 29 from mouse SC5. Four of these mice have been used for tumor samples and the other four mice have been used for the luminal samples.

#After talking with Sally about this data set, we came to a conclusion that the best approach for this assignment would be to randomly select three samples from each mouse for my further analysis. The next section of the code depicts this.

#I am going to duplicate the data frame to subset my data.
dfData <- dfComplete_Data
#Check dimensions.
dim(dfData)

#The first step is to make the columns displaying the gene names into the row names.
rownames(dfData) <- dfData$X
#Check to see if it worked
rownames(dfData)

#Delete the gene name column.
dfData <- dfData[-1]
#Check to see if it worked.
dim(dfData)

sample(dfData, 10)
#for loop for column names








#### 4- MAIN SOFTWARE TOOLS ----










#### 5- MAIN ANALYSIS ----













#### 6- RESULTS AND DISCUSSION ----












#### 7- ACKNOWLEDGEMENTS ----








#### 8- REFERENCES ----

#https://master.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
#RNASeq Bioconductor vignette

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148569
#Data

#Sun H, Zeng J, Miao Z, Lei KC, Huang C, Hu L, Su SM, Chan UI, Miao K, Zhang X, Zhang A, Guo S, Chen S, Meng Y, Deng M, Hao W, Lei H, Lin Y, Yang Z, Tang D, Wong KH, Zhang XD, Xu X, Deng CX. Dissecting the heterogeneity and tumorigenesis of BRCA1 deficient mammary tumors via single cell RNA sequencing. Theranostics 2021; 11(20):9967-9987.
#Paper of data








### QUESTIONS TO ASK

#1- There are 477 samples and 21004 genes. Should I choose 1-2 samples from each SC?
#Categories are tumor and luminal (4 SCs from each)