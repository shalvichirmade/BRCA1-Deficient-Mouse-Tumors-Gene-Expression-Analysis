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
complete_data <- read.delim("GSE148569_scRNAseq_C1_count.txt")

#Inspect the data.
head(complete_data, n = c(10, 5))
#There are 478 samples in this data set and 21003 genes analyzed.





#### 4- MAIN SOFTWARE TOOLS ----










#### 5- MAIN ANALYSIS ----













#### 6- RESULTS AND DISCUSSION ----












#### 7- ACKNOWLEDGEMENTS ----








#### 8- REFERENCES ----