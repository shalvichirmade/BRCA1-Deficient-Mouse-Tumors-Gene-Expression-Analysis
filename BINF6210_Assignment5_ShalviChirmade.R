#### BINF 6210 - Assignment 5 - Due Friday December 17, 2021 by 5 pm ----
# By Shalvi Chirmade


### TITLE

#### 1- INTRODUCTION ----








#### 2- DATA SET ----









#### 3- DATA ACQUISITION, EXPLORATION, FILTERING AND QUALITY CONTROL ----

##If needed, install the packages below by removing the comment on the specific line. Otherwise, please load in these required packages.

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(c("limma", "Glimma", "edgeR"))
library(edgeR)
library(Glimma)
library(limma)


#install.packages("gplots")
library(gplots)
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("R.utils")
library(R.utils)
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
rm(dfSC1_, dfSC2_, dfSC4_, dfSC5_, dfSC6_, dfSC7_, dfSC9_, dfSC10_,name)

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

#Extract column names for later analysis. 
samplenames <- colnames(DGE_Data$counts)

#Remove variables that are no longer needed.
rm(dfData)


#Check to see if there any any duplicated gene names present in our data set.
anyDuplicated(rownames(DGE_Data$counts)) #0

#As the data set I acquired already had gene names for each row, I did not need to extract these values using the Mus.musculus package.

#There are also various gene names in the data set that consist of numbers and then the name Riken; for example, 4931408C20Rik. I do not have much experience working with the mouse genome but according to my research, these belong to names of genes that have not been officially annotated (Hayashizaki, 2003). Riken seems to be a large genomic institute in Japan that collaborates with various other institutes across the world to form comprehensive research on multiple different organisms. When I enter these names in UniProt, the associated gene name are the Rik ones we see here (https://www.uniprot.org/uniprot/E9PWP9).

#In general for gene expression analysis, raw counts are not considered accurate values for comparison. Instead, counts per million (cpm) or log2-counts per million (lcpm) are utilized. This allows for the library size to be considered in the creation of these values. Library size is the total number of counts from each sample.
#The total library size for our subset data set is:
sum(DGE_Data$samples$lib.size) #39624581
#The average library size for each sample is:
mean(DGE_Data$samples$lib.size) #1651024
#Will display library size as a bar plot a little later on.

#Converting the counts to cpm and lcpm.
dfCPM <- cpm(DGE_Data)
dfLCPM <- cpm(DGE_Data, log = T)
#You can see that in the cpm matrix, all the 0 gene expression values remained 0 but in the lcpm matrix, this has been transformed to a more relational data set; it is useful when creating exploratory plots. It reduces the inter-sample changes and prevents a large separation between the lowest and highest count values. 

#Calculating the L and M parameters used in the lcpm calculations; it will be used for generating read density figures later on.
L <- mean(DGE_Data$samples$lib.size) * 1e-6
M <- median(DGE_Data$samples$lib.size) * 1e-6
#The library size for this data set is very very low in comparison to the vignette data set. This has about 1.6 million on average while the vignette data set was at 45.5 million. We will see how this is affected when we compare the genes lost during filtration in the next step.

#Create a bar plot showing the difference in library size.
#The average library size for each mouse.
SC1_lib <- mean(DGE_Data$samples$lib.size[1:3])
SC2_lib <- mean(DGE_Data$samples$lib.size[4:6])
SC4_lib <- mean(DGE_Data$samples$lib.size[7:9])
SC6_lib <- mean(DGE_Data$samples$lib.size[10:12])
SC5_lib <- mean(DGE_Data$samples$lib.size[13:15])
SC7_lib <- mean(DGE_Data$samples$lib.size[16:18])
SC9_lib <- mean(DGE_Data$samples$lib.size[19:21])
SC10_lib <- mean(DGE_Data$samples$lib.size[22:24])

#Create a data frame containing the information needed to create a visualization.
dfLib <- data.frame(name = SC_Names, size = rbind(SC1_lib, SC2_lib, SC4_lib, SC6_lib, SC5_lib, SC7_lib, SC9_lib, SC10_lib), cell_type = c("tumor","tumor", "tumor", "tumor", "luminal", "luminal", "luminal", "luminal"))

ggplot(dfLib, aes(x = factor(name, levels = SC_Names), y = size, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Cell-type", values = c("#47E0E5", "#62E547")) +
  theme_bw() +
  labs(x = "Mouse Name", y = "Library Size", title = "Library Size in Samples") 
#The distribution of these library sizes on average per cell-type is not noticeably different. There is one sample from each that has a lower library size; this can be due to the fact that I randomly selected three samples from each mouse. If we had used the full data set, we may not see this distribution.


#Let's determine if there are any genes that have zero expression in all samples.
table(rowSums(DGE_Data$counts == 0) == 24)
#There are 7429 genes with zero expression. These will be deleted from our data set as they can hinder our downstream analysis.

#This function, filterByExpr from the edgeR package, allows us to filter out the lowly expressed genes while still maintaining a large amount of data for analysis. It creates a logical vector displaying the rows to be kept and removed.
KeepExprs <- filterByExpr(DGE_Data, group = group)
DGE_Data <- DGE_Data[KeepExprs,, keep.lib.sizes = F]

#Let's compare the amount of genes we lost to the number that had zero counts in all samples.
dim(DGE_Data) #3382 24
(21004-3382)/21004 * 100 #Removed 83.9% of the genes..
#It has removed a SUBSTANTIAL amount of genes! After re-analyzing the original data set, I came to realize that most of the rows had values of 0 throughout each gene. This could be because I have only taken a small subset of the samples actually analyzed by the authors of the paper or it could be the due to the quality of reads attained by the authors while carrying out the experiment. I am going to continue my analysis using the 3400 genes I do have in this data as it should still be able to yield functional plots for my interpretation. If I am obstructed with errors, I will have to re-evaluate this data set or the filtration step itself. 

#I will re-try this step and reduce the minimum count to 5 as the default is set to 10. I hope we notice a fewer amount of genes being discarded.
# KeepExprs <- filterByExpr(DGE_Data, group = group, min.count = 5)
# DGE_Data <- DGE_Data[KeepExprs,, keep.lib.sizes = F]
# dim(DGE_Data) #3914 24
# (21004-3914)/21004 * 100 #81.4%
#As the number of samples barely reduced. I will keep the default minimum count as it is recommended by the vignette for a more accurate downstream analysis. The vignette also states that the genes retained have counts in mostly all samples for the same grouping. In our case, these groupings would be luminal and tumor. So if a gene has multiple zero counts for all the luminal or tumor samples, this gene would be disregarded. If the gene is of interest to the study, most of the samples from the same groupings should have a count associated. If not, then the count for that gene could be a rogue value. After reading this, I re-checked to make sure I made the correct groupings based on the paper and the GEO database; as the authors did not specifically specify what each sample corresponds to, my speculations could be inaccurate. Furthermore, the vignette data set lost about 60% of their data during this filtration step so my assumptions could be the right choice. They also mention that a lower library size can be a factor of losing more data as there is less information to evaluate. Let's move on and see what our current data set can provide us.

#I will now produce a figure comparing the density of reads from the raw unfiltered data and the filtered data. This code comes from the vignette for RNA-Seq analysis by Bioconductor.

lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(DGE_Data)
col <- colorRampPalette(brewer.pal(12, "Paired")) (nsamples) #Add more colors to the palette.
par(mfrow=c(1,2))
plot(density(dfLCPM[,1]), col=col[1], lwd=2, ylim=c(0,2), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(dfLCPM[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topleft", samplenames, text.col=col, bty="n", cex = 0.75)

dfLCPM_filtered <- cpm(DGE_Data, log=TRUE)
plot(density(dfLCPM_filtered[,1]), col=col[1], lwd=2, ylim=c(0,2), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(dfLCPM_filtered[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topleft", samplenames, text.col=col, bty="n", cex = 0.75) #legends don't show in zoom

#We can see a drastic difference between the raw data and the filtered. A large amount of data were at zero to low values. We can also see a varying expression for each sample; this can account to the quality of the data from this study.



#The next step in our quality control is normalization.

#Normalization is carried out to allow for a reasonably similar expression pattern among all the samples being analyzed. As we saw in the bar plot, there were a few samples with lower library size and one with almost half of the largest library. This can occur from external factors during the experiment such as sequencing depth and even GC count (Evans et al., 2017). The main reason for normalization is to showcase the true differences in samples instead of using the raw values given by the instrument (Evans et al., 2017). Another reason for the use normalization is that the variance in gene counts can be from the outcome of differential coverage instead of differential gene expression (Evans et al., 2017). THe vignette uses the method of Trimmed Mean of the M-values (TMM) where a sample is chosen as a reference and fold-changes are calculated relative to this sample (Evans et al., 2017). Apparently this method of usage is known only from being used in the edgeR package, which is the one we are using today. This method did not give a good representation of normalized data so I chose TMMwsp instead which is 'TMM with singleton pairing' method. This performs better with data consisting of multiple 0 expression values; as we have seen, this data set is riddled with 0 expression, hence it works better for my further analysis.

DGE_Data <- calcNormFactors(DGE_Data, method = "TMMwsp")
#We see a new element in our DGE_Data object called norm.factors.
DGE_Data$samples$norm.factors

#To better visualize the impact of normalization, we will create a boxplot showcasing log-CPM values of expression distribution for both unnormalized and normalized data. We will see that in the unnormalized data, there will be a variance in expression levels but this will be standardized in the normalized image. The DGE object is duplicated to manipulate the already normalized data into showcasing what unnormalized data looks like.
DGE_Data_2 <- DGE_Data
DGE_Data_2$samples$norm.factors <- 1 #Setting all the norm factors to 1
DGE_Data_2$counts[,1] <- ceiling(DGE_Data_2$counts[,1] * 0.05)
DGE_Data_2$counts[,2] <- DGE_Data_2$counts[,2] * 5

par(mfrow = c(1,2))
dfLCPM <- cpm(DGE_Data_2, log = T)
boxplot(dfLCPM, las = 2, col = col, main = "A. Unnormalizaed Data Example", ylab = "Log-cpm") #Using the same colors as the plot above
DGE_Data_2 <- calcNormFactors(DGE_Data_2, method = "TMMwsp") #Normalizing data in this duplicted DGE object
dfLCPM <- cpm(DGE_Data_2, log = T)
boxplot(dfLCPM, las = 2, col = col, main = "A. Normalizaed Data Example", ylab = "Log-cpm")
par(mfrow = c(1,1))

#We can see that the mean values of all the samples still do not line up as well as the example data set in this vignette, however, this could be a result of the quality of data collection done by the authors of the paper or even the mice themselves. We can see a noticeable difference in mouse 10 compared to the others. As the authors did not lay out the exact differences between the mice tested, my speculation is that this mouse was affected by some other external factor compared to the other mice in its category.


# #Unsupervised clustering of samples - DELETE
# 
# #A multi-dimensional scaling (MDS) plot can be created to visualize the similarities and differences in the samples being used. This can be done using the plotMDS function in limma.
#dfLCPM <- cpm(DGE_Data, log = T)
# Color_By_Group <- group
# levels(Color_By_Group) <- brewer.pal(2, "Spectral") #Colorblind friendly
# Color_By_Group <- as.character(Color_By_Group)
# plotMDS(dfLCPM, labels = group, col = levels(Color_By_Group), main = "Sample Groups")



#### 4- MAIN SOFTWARE TOOLS ----




#### 5- MAIN ANALYSIS ----

#The main analysis of this script is showcasing the genes that are differentially expressed among the two groupings of our data, tumor and luminal. For this, we are in the assumption that our data is normally distributed. First we create a design matrix for the cell-type groupings. I used the vignette (https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html) to help me determine if I need to use an intercept term. I will not include an intercept term as the grouping are not covariates and the models are considered equivalent.
mtDesign <- model.matrix(~0+group)
colnames(mtDesign) <- gsub("group", "", colnames(mtDesign))
mtDesign

#We have to set up model contrasts for the our groups as they do not have an associated intercept. The makeContrasts function from limma is used to make pairwise comparison between the two cell populations.
mtContrasts <- makeContrasts(TvsL = tumor - luminal, levels = colnames(mtDesign))
mtContrasts

#Honestly, I found going along the vignette quite difficult for my data set. As I am not completely sure if my groupings are accurate, I am always worried that my results are incorrect. However, as this assignment is about conducting an analysis along with learning the ups and downs about data sets and packages, I believe that I am learning a lot along the way.


#Removing heteroscedasticity from count data

#When using raw counts for RNA-seq data, the variance is not considered independent from the mean. For this script, we use the lcpm values where we assume our data is normally distributed. The function, voom, calculates precision weights on the mean-variance dependency by using library size and the normalization factors.
par(mfrow = c(1,2))
V_EList <- voom(DGE_Data, mtDesign, plot = T)
VFit <- lmFit(V_EList, mtDesign)
VFit <- contrasts.fit(VFit, contrasts = mtContrasts)
EFit <- eBayes(VFit)
plotSA(EFit, main = "Final Model: Mean-variance trend")

#In these plots, the means are plotted on the x axis and the variances are plotted on the y axis. The comparison image shows the distribution before and after voom is applied to the data set. After briefly reading the paper introducing voom (Law et al., 2014), I was unable to understand why the data set before voom has a slight increase trend before decreasing. My understanding of linear models is not adequate enough to have an explanation; if this can be explained, I would really like to understand why. The black dots in each of these plots correlates to a gene in our data set. We can see that this data set is smaller in comparison to the data set used in the vignette. Voom creates an EList object containing various information we have already seen in the DGEList object as well as additional information such as expression values and precision weights.


#Examining the number of differentially expressed genes

#A simple table can be created using the decideTests function from limma to summarize the number of significantly up- and down-regulated genes. As per usual, significance is cut off at a p-value of 0.05 by default.
dt <- decideTests(EFit)
summary(dt)
#We can see that most of the genes are not significantly different from the two different groupings. This tells us that there are 28 genes that are up-regulated in the tumor samples and 26 down-regulated genes in the tumor samples relative to the luminal samples. This can mean that there only a few gene expressions that are affected in our data set when comparing tumor and luminal samples of mice.

#This step allows us to extract the differentially expressed gene names. As 0 represents genes that are not differentially expressed, we are looking for the rest.
DE_Common <- which(dt[,1] != 0)
length(DE_Common) #54 genes
DE_Genes <- rownames(EFit$coefficients)[DE_Common] #These are the significantly expressed genes. Storing this vector to use for finding GO terms later on.
DE_Genes

# TFit <- treat(VFit, lfc = 1)
# dt2 <- decideTests(TFit)
# summary(dt2) #7 up and 1 down regulated genes


#Displaying the differential expression results graphically

#Limma has a function called plotMD which allows for the genes to be displayed visually by mean-difference plots. This uses log-FC values from our lineal model analysis against the mean lcpm value.
par(mfrow = c(1,1))
plotMD(EFit, column = 1, status = dt[,1], main = colnames(EFit)[1], xlim = c(0, 13))

#A MA plot is an interpretation of the Bland-Altman plot which calculates the agreement between quantitative measures (Giavarina, 2015). The mean and sd of the differences between the measurements are used to calculate the statistical limits (Giavarina, 2015). It is a scatter plot where the x axis displays average of the means and the y axis measures the difference of the paired measurements we calculated earlier when creating our contrast matrix (Giavarina, 2015). As we can see in the colored dots, the red ones display the up-regulated genes and the blue dots show the down-regulated genes.

#A heatmap can also be used to display the differentially expressed genes. We will be using all 54 DE genes. We can utilize the heatmap.s function from the gplots package to create this.
dfLCPM <- cpm(DGE_Data, log = T)
Color_HeatMap <- colorpanel(1000, "blue", "white", "red")

heatmap.2(dfLCPM[DE_Common,], 
          scale = "row",
          labRow = DE_Genes,
          labCol = group,
          col = Color_HeatMap,
          trace = "none",
          density.info = "none",
          lhei = c(2,10),
          dendrogram = "column")

#As we can see, this data set is quite a mess. There are two tumor samples that have clustered with the luminal samples. As I talked about in the Data Acquisition section, both these cell-types are from a BRCA1-deficient mouse; the paper was studying the intratumor variations so a clear distinction of gene expression is not going to be seen. However, we can mildly see a distinction between the up-regulated and down-regulated genes in the tumor and luminal samples. Most of the genes at the top of the plot show an up-regulation in tumor cells and the most of the genes at the bottom of the plot show an up-regulation in the luminal cells. The next step for my analysis is to correlate these differentially expressed genes with their associated GO terms; hopefully we can figure out why the particular genes are either up- or down-regulated in our two groups.



#### 6- RESULTS AND DISCUSSION ----












#### 7- ACKNOWLEDGEMENTS ----








#### 8- REFERENCES ----

#https://master.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
#RNASeq Bioconductor vignette

#https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html
#Design matrices vignette

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148569
#Data

#Sun H, Zeng J, Miao Z, Lei KC, Huang C, Hu L, Su SM, Chan UI, Miao K, Zhang X, Zhang A, Guo S, Chen S, Meng Y, Deng M, Hao W, Lei H, Lin Y, Yang Z, Tang D, Wong KH, Zhang XD, Xu X, Deng CX. Dissecting the heterogeneity and tumorigenesis of BRCA1 deficient mammary tumors via single cell RNA sequencing. Theranostics 2021; 11(20):9967-9987.
#Paper of data

#Hayashizaki Y. (2003). RIKEN mouse genome encyclopedia. Mechanisms of ageing and development, 124(1), 93–102.

#Evans, C., Hardin, J., & Stoebel, D. M. (2018). Selecting between-sample RNA-Seq normalization methods from the perspective of their assumptions. Briefings in bioinformatics, 19(5), 776–792.

#Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome biology, 15(2), R29.

#Giavarina D. (2015). Understanding Bland Altman analysis. Biochemia medica, 25(2), 141–151.







### QUESTIONS TO ASK
