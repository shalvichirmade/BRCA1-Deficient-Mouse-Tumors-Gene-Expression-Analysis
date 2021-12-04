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
head(dfComplete_Data, n = c(10, 5))
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


#Create a variable with the different SC names.
SC_Names <- c("SC1_", "SC2_", "SC4_", "SC5_", "SC6_", "SC7_", "SC9_", "SC10_")

#Randomly sample three columns for each SC mouse.
#sample(dfData, 10) #randomly samples 10 columns

#Create a function to randomly select three samples from each mouse. As I have to do this eight times, I thought creating a function would be beneficial. This function can be used for any data frame, for any string to be matched and any number of samples. It will output a new data frame with the name of the string entered in the function.
Randomly_Sample_Columns <- function(df, str, n){
  
  #df - the name of the data frame from which you would like to sample columns
  
  #str - the string that you want to search for in the column names. It must be entered with ""
  
  #n - number of columns with the string you would like to randomly sample
  
  #Extract the string to use as the name of the data frame to be created
  name <- paste("df", str, sep = "")
  
  dfOutput <- df %>%
    select(matches(str)) %>%
    sample(n)
  
  #Add dfOutput to the global enviornment using the str as the name of the data frame
  assign(name, data.frame(dfOutput), envir = .GlobalEnv)
  
}

#Using a for loop, create new data frames for each mouse with the help of the new function.
set.seed(6210)

for (name in SC_Names) {
  
  Randomly_Sample_Columns(dfData, name, 3)
  
}

#Combine these new data frames to use for the remainder of the analysis.
dfData <- cbind(dfSC1_, dfSC2_, dfSC4_, dfSC6_, dfSC5_, dfSC7_, dfSC9_, dfSC10_)
#This order is dependent on the type of sample; the first four are tumor and the last four are luminal.

#Remove the data frames no longer needed.
rm(dfSC1_, dfSC2_, dfSC4_, dfSC5_, dfSC6_, dfSC7_, dfSC9_, dfSC10_, SC_Names, name)

#Convert data frame into a DGE object for gene expression analysis. This cn be done using the DGEList function from edgeR. The input data frame is the one we just created, dfData and the groupings are the cell type, so either luminal or tumor.

#This vector is to label each mouse sample type, luminal or tumor.
group <- rep(c("tumor", "luminal"), each = 12)

#Create the DGEList object.
DGE_Data <- DGEList(dfData, group = group)
#Check to see if this worked.
DGE_Data
#Can view each individual element of this list as well.
View(DGE_Data$counts)
View(DGE_Data$samples)

#Remove variables that are no longer needed.
rm(dfData, group)



#### 4- MAIN SOFTWARE TOOLS ----

#Check to see if there any any duplicated gene names present in our data set.
anyDuplicated(rownames(DGE_Data$counts)) #0

#As the data set I acquired already had gene names for each row, I did not need to extract these values using the Mus.musculus package.

#There are also various gene names in the data set that consist of numbers and then the name Riken; for example, 4931408C20Rik. I do not have much experience working with the mouse genome but according to my research, these belong to names of genes that have not been officially annotated (Hayashizaki, 2003). Riken seems to be a large genomic institute in Japan that collaborates with various other institutes across the world to form comprehensive research on multiple different organisms. When I enter these names in UniProt, the associated gene name are the Rik ones we see here (https://www.uniprot.org/uniprot/E9PWP9).

#In general for gene expression analysis, raw counts are not considered accurate values for comparison. Instead, counts per million (cpm) or log2-counts per million (lcpm) are utilized. This allows for the library size to be considered in the creation of these values. Library size is the total number of counts from each sample.
#The total library size for our subset data set is:
sum(DGE_Data$samples$lib.size) #39624581
#The average library size for each sample is:
mean(DGE_Data$samples$lib.size) #1651024

#Converting the counts to cpm and lcpm.
dfCPM <- cpm(DGE_Data)
dfLCPM <- cpm(DGE_Data, log = T)
#You can see that in the cpm matrix, all the 0 gene expression values remained 0 but in the lcpm matrix, this has been transformed to a more relational data set.







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

#Hayashizaki Y. (2003). RIKEN mouse genome encyclopedia. Mechanisms of ageing and development, 124(1), 93–102.






### QUESTIONS TO ASK

#1- There are 477 samples and 21004 genes. Should I choose 1-2 samples from each SC?
#Categories are tumor and luminal (4 SCs from each)