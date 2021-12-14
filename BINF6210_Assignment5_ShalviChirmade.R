#### BINF 6210 - Assignment 5 - Due Friday December 17, 2021 by 5 pm ----
# By Shalvi Chirmade

### Gene Expression Analysis on BRCA1-Deficient Mouse Cells

#GitHub link: https://github.com/shalvichirmade/BRCA1-Deficient-Mouse-Tumors-Gene-Expression-Analysis

#### 1- INTRODUCTION ----

#I began my interest in genetic counselling when I was in high school. This came from a small paragraph in my biology textbook and evolved into a career passion. As you may know, genetic counselors are health professionals that help guide patients with understanding and providing information on their genetic reports and/or disorders. This can range from fields of early developmental and birth defects to inherited disorders and familial cancer. I personally have volunteered and worked with genetic counselors after graduating from my BSc and continued to be entranced by the work and research they conduct. My goal was to one-day become a genetic counselor myself until I found the field of bioinformatics and realized I could turn my fascination of human genetic disease into computational analysis. This is the main reason I decided to search for a data set analyzing breast cancer tumors. A comprehensive understanding of how BRCA1 affects the formation of breast and/or ovarian tumors is not complete so research of this kind is important for society. It aids in acquiring knowledge of hereditary and sporadic carcinogenics (Yoshida & Miki, 2004). 
#BRCA1 is an important gene responsible for DNA repair, cell cycle checkpoints and response to DNA damage (Yoshida & Miki, 2004). It also plays a crucial role in genomic stability; the protein synthesized interacts with its proteome to allow for recognition and repairing of damaged DNA (Jhanwar-Uniyal, 2003). As it interacts with multiple different proteins to carry out this process, it plays an very important role in suppressing tumorigenesis (Jhanwar-Uniyal, 2003). BRCA1 is a well-conserved gene among various species, research in hope to understand its biological process is crucial for humans and other related organisms (Jhanwar-Uniyal, 2003). Some mutations in BRCA1 mask the functional ability of the wild-type allele which causes a higher probability of the individual developing a breast cancer (Jhanwar-Uniyal, 2003). Mutations of tumor-inducing genes, such as BRCA1, allow for an early-onset of breast and/or ovarian cancer in women (Sun et al., 2021). The research conducted by these authors investigate the origin of mammary tumors from the luminal epithelia and compare the gene expression to cells derived from the mammary tumors themselves. The data set derived from this study is what I will be using for this script; it will be described in greater detail in the next section. As for the gene expression analysis being conducted in this script, choosing the correct software tools for statistical analysis is of great importance. I will go into deeper detail about this later on when explaining which software tool I decided to use.
#The question I would like evaluate for this script is, what are the biological processes of the differentially expressed genes in this data set? Are the DE genes related to the biological function of BRCA1, i.e., are they part of the individual networks for DNA repair and damage or are their functions on a broader spectrum like immune response or apoptosis? I plan to carry out gene expression analysis using the data set obtained by Sun et al. by utilizing the packages limma and goseq by Bioconductor.




#Genetic counselling
#Why choose BRCA1 study
#Difficulty with data - brief
#Difference in luminal and tumor cells
#DE question - BRCA1 regulatory proteins DE?
#What you expect based on your knowledge






#### 2- DATA SET ----

#I was searching through the GEO database when I came across this data set (November 22, 2021). It comes from a very interesting paper written by Sun et al. in 2021 titled "Dissecting the heterogeneity and tumorigenesis of BRCA1-deficient mammary tumors via single cell RNA sequencing". This data set was obtained and used by the authors for expression profiling and the understanding of the BRACA1-deficient mouse tumors via RNA sequencing. BRCA1 is a gene known to be tumor-inducing when mutated (Sun et. al., 2021); this can occur either sporadically or be inherited from a parent. This mutation can predispose breast and/or ovarian cancers (Yoshida & Miki, 2004). As the normal function of this gene is involved with DNA repair and apoptosis (Yoshida & Miki, 2004), we can understand how dysfunctional and detrimental a mutation in these regulatory processes can be. BRCA1 is aided by many regulatory proteins and it will be interesting to see if this data set shows differential expression in genes that interact with BRCA1. As all the mice utilized in this data set were BRCA1-deficient (the same cross, strain: C57BL/6 and 129/Sv mixed), the analysis carried out by the authors were on the differential expression between the luminal and tumor cells from their mammary glands. The authors main analysis was to showcase the heterogeneity of these cells in both intra- and inter-tumor levels. They wanted to assess the initiation and progression of tumor formation hence isolating samples from the luminal cells of the mammary gland and the tumor cells from the carcinoma of the same region. The data available on GEO has 477 samples which encompasses samples from eight different mice. Each of these mice are replicates from the same cross I mentioned earlier. There are between 30-60 samples from each mouse which I will discuss in the next section when carrying out my data exploration. As you will see, I found understanding this data set quite difficult. In the published paper, the authors do not correspond the names of the mice chosen to their publicly available data set hence the exact difference between each mouse was very hard to comprehend. The only difference I was able to find, once I went through about 100 of the 477 samples, was only if the sample came from luminal or tumor regions. As the authors also take into account gestation of the mouse, this information was not differentiated in the public data set. For each sample, there is a separate accession display but the description are all the same, which states "Three-to-four month-old female virgin or pregnant at day 12.5 mice were sacrificed". Hence, I was unaware if the mouse was either virgin or 12.5 days pregnant. This confusion carries on during the rest of my analysis as we will see that the differentially expressed genes are very minimal and not highly categorized. My speculation is, due to the fact that I am not explicitly aware of the difference between each sample, I may have grouped the samples in a less efficient way than they were intended. I have only differentiated the samples by luminal or tumor and have not taken into account gestation. I will discuss this further in my final results based on the analysis I accomplish.

#The data set can be found at:
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148569


#### 3- DATA ACQUISITION, EXPLORATION, FILTERING AND QUALITY CONTROL ----

##If needed, install the packages below by removing the comment on the specific line. Otherwise, please load in these required packages.

#BIOCONDUCTOR
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(c("limma", "Glimma", "edgeR", "goseq", "Mus.musculus","TxDb.Mmusculus.UCSC.mm39.refGene"))
library(edgeR)
library(Glimma)
library(goseq)
library(limma)
library(Mus.musculus)
library(TxDb.Mmusculus.UCSC.mm39.refGene)

#CRAN
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
head(rownames(dfData),10)

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
  
  #Add dfOutput to the global environment using the str as the name of the data frame
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

#Convert data frame into a DGE object for gene expression analysis. This cn be done using the DGEList function from edgeR. The input data frame is the one we just created, dfData, and the groupings are the cell type, so either luminal or tumor.

#This vector is to label each mouse sample type, luminal or tumor.
group <- rep(c("tumor", "luminal"), each = 12)

#Create the DGEList object.
DGE_Data <- DGEList(dfData, group = group)
#Check to see if this worked.
DGE_Data
#Can view each individual element of this list as well. 
#TODO - view small amount for RM
head(DGE_Data$counts, 5)
head(DGE_Data$samples, 5)

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
#Let's look at them.
head(dfCPM, 5)
head(dfLCPM, 5)
#You can see that in the cpm matrix, all the 0 gene expression values remained 0 but in the lcpm matrix, this has been transformed to a more relational data set; it is useful when creating exploratory plots. It reduces the inter-sample changes and prevents a large separation between the lowest and highest count values. 

#Calculating the L and M parameters used in the lcpm calculations; it will be used for generating read density figures later on.
L <- mean(DGE_Data$samples$lib.size) * 1e-6
M <- median(DGE_Data$samples$lib.size) * 1e-6
#The library size for this data set is very very low in comparison to the vignette data set. This has about 1.6 million on average while the vignette data set was at 45.5 million. We will see how this is affected when we compare the genes lost during filtration in the next step.

#Create a bar plot showing the difference in library size.
#The average library size for each mouse. I was going to make this via lapply but I was confused as to how to change the ranges with each iteration without using a for loop.
SC1_lib <- mean(DGE_Data$samples$lib.size[1:3])
SC2_lib <- mean(DGE_Data$samples$lib.size[4:6])
SC4_lib <- mean(DGE_Data$samples$lib.size[7:9])
SC6_lib <- mean(DGE_Data$samples$lib.size[10:12])
SC5_lib <- mean(DGE_Data$samples$lib.size[13:15])
SC7_lib <- mean(DGE_Data$samples$lib.size[16:18])
SC9_lib <- mean(DGE_Data$samples$lib.size[19:21])
SC10_lib <- mean(DGE_Data$samples$lib.size[22:24])

#Create a data frame containing the information needed to create a visualization.
dfLib <- data.frame(name = SC_Names, 
                    size = rbind(SC1_lib, SC2_lib, SC4_lib, SC5_lib, SC6_lib, SC7_lib, SC9_lib, SC10_lib), 
                    cell_type = c("tumor","tumor", "tumor", "luminal", "tumor", "luminal", "luminal", "luminal"))

dfLib %>% ggplot(aes(x = factor(name, levels = SC_Names), y = size, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Cell-type", values = c("#47E0E5", "#62E547")) +
  theme_bw() +
  labs(x = "Mouse Name", y = "Library Size", title = "Library Size in Samples") 
#The distribution of these library sizes on average per cell-type is not noticeably different. There is one sample from each that has a lower library size; this can be due to the fact that I randomly selected three samples from each mouse. If we had used the full data set, we may not see this distribution.

#I want to see the average library size for each grouping and see if we have  similar numbers.
mean(SC5_lib, SC7_lib, SC9_lib, SC10_lib) #luminal - 1504943
mean(SC1_lib, SC2_lib, SC4_lib, SC6_lib) #tumor - 1905406
#Their average is very similar so we won't worry about the individual variation.


#Let's determine if there are any genes that have zero expression in all samples.
table(rowSums(DGE_Data$counts == 0) == 24)
#There are 7429 genes with zero expression. These will be deleted from our data set as they can hinder our downstream analysis.

#The function, filterByExpr from the edgeR package, allows us to filter out the lowly expressed genes while still maintaining a large amount of data for analysis. It creates a logical vector displaying the rows to be kept and removed.
KeepExprs <- filterByExpr(DGE_Data, group = group)
DGE_Data <- DGE_Data[KeepExprs,, keep.lib.sizes = F]

#Let's compare the amount of genes we lost to the number that had zero counts in all samples.
dim(DGE_Data) #3382 24
(21004-3382)/21004 * 100 #Removed 83.9% of the genes..
#It has removed a SUBSTANTIAL amount of genes! After re-analyzing the original data set, I came to realize that most of the rows had values of 0 throughout each gene. This could be because I have only taken a small subset of the samples actually analyzed by the authors of the paper or it could be the due to the quality of reads attained by the authors while carrying out the experiment. I am going to continue my analysis using the 3382 genes I do have in this data as it should still be able to yield functional plots for my interpretation. If I am obstructed with errors, I will have to re-evaluate this data set or the filtration step itself. 

#I will re-try this step and reduce the minimum count to 5 as the default is set to 10. I hope we notice a fewer amount of genes being discarded.
# KeepExprs <- filterByExpr(DGE_Data, group = group, min.count = 5)
# DGE_Data <- DGE_Data[KeepExprs,, keep.lib.sizes = F]
# dim(DGE_Data) #3914 24
# (21004-3914)/21004 * 100 #81.4%
#As the number of samples barely increased, I will keep the default minimum count as it is recommended by the vignette for a more accurate downstream analysis. The vignette also states that the genes retained have counts in mostly all samples for the same grouping. In our case, these groupings would be luminal and tumor. So if a gene has multiple zero counts for all the luminal and tumor samples, this gene would be disregarded. If the gene is of interest to the study, most of the samples from the same groupings should have a count associated. If not, then the count for that gene could be a rogue value. This would mean that the expression of the particular gene cannot depict a reliable conclusion as it may not play a major role in our analysis. After reading this, I re-checked to make sure I made the correct groupings based on the paper and the GEO database; as the authors did not specifically specify what each sample corresponds to, my speculations could be inaccurate. Furthermore, the vignette data set lost more than 60% of their data during this filtration step so my assumptions could be the right choice. They also mention that a lower library size can be a factor of losing more data as there is less information to evaluate. Let's move on and see what our current data set can provide us.

#Remove variables that are not required.
rm(SC1_lib, SC2_lib, SC4_lib, SC5_lib, SC6_lib, SC7_lib, SC9_lib, SC10_lib, KeepExprs)

#I will now produce a figure comparing the density of reads from the raw unfiltered data and the filtered data. This code is adapted from the vignette for RNA-Seq analysis by Bioconductor.
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

#We can see a drastic difference between the raw data and the filtered. A large amount of data were at zero to low values. We can also see a varying expression for each sample; this can account to the quality of the data from this study and the difference in their individual library size.



#The next step in our quality control is normalization.

#Normalization is carried out to allow for a reasonably similar expression pattern among all the samples being analyzed. As we saw in the bar plot, there were a few samples with lower library size and one with almost half of the largest library. This can occur from external factors during the experiment such as sequencing depth and even GC count (Evans et al., 2017). The main reason for normalization is to showcase the true differences in samples instead of using the raw values given by the instrument (Evans et al., 2017). Another reason for the use of normalization is that the variance in gene counts can be from the outcome of differential coverage instead of differential gene expression (Evans et al., 2017). The vignette uses the method of Trimmed Mean of the M-values (TMM) where a sample is chosen as a reference and fold-changes are calculated relative to this sample (Evans et al., 2017). Apparently this method of usage is known only from being used in the limma package, which is the one we are using today. This method did not give a good representation of normalized data so I chose TMMwsp instead which is 'TMM with singleton pairing' method. This performs better with data consisting of multiple 0 expression values; as we have seen, this data set is riddled with 0 expression, hence it works better for my further analysis.

DGE_Data <- calcNormFactors(DGE_Data, method = "TMMwsp")
#We see a new element in our DGE_Data object called norm.factors.
DGE_Data$samples$norm.factors

#TODO comment out image
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

#When plotted, we can see that the mean values of all the samples still do not line up as well as the example data set in this vignette, however, this could be a result of the quality of data collection done by the authors of the paper or even the mice themselves. We can see a noticeable difference in mouse 10 compared to the others. As the authors did not lay out the exact differences between the mice tested (gestation period), my speculation is that this mouse was affected by some other external factor compared to the other mice in its category.

#I had also created a multi-dimensional scaling(MDS) visualization to display the similarities and differences in the samples being used. This was done using the plotMDS function in limma. However, I decided to delete the plot from the assignment as there was no clear distinction between the two groupings of the samples. Two of the tumor samples were grouped together with the luminal samples; as I mentioned before, I am unable to differentiate the gestation period of each mouse sample, so this could be a reason as to why two samples grouped irregularly.

#Remove variables that are not required for the rest of the analysis.
rm(DGE_Data_2, dfLCPM_filtered, dfLib, L, M, SC_Names, samplenames)


#### 4- MAIN SOFTWARE TOOLS ----

#The main vignettes I decided to adapt for this script are both from Bioconductor: RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR, and goseq: Gene Ontology testing for RNA-seq datasets.
#https://master.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
#https://bioconductor.org/packages/devel/bioc/vignettes/goseq/inst/doc/goseq.pdf
#TODO make into bullet points
#To begin, I researched different RNA sequencing packages available and used the paper by Seyednasrollah et al. as a guideline for selecting a specific software tool. To explain briefly, the paper outlines the most popular software tools currently used for analyzing RNA-sequencing data. The authors used real data to conduct their analysis instead of using simulated data; this accounts for natural variability seen in real-life data and allows the comparison to be more realistic. As the human understanding of variability is imperfect, using this type of real data becomes a large advantage while comparing computational tools. Through their analysis, I decided to use the tool limma by Bioconductor as this package had a high precision level which corresponds to the accuracy of analyzing replicates. The authors also display a proportion of false discoveries for the software tools analyzed and limma displays a low proportion. As I also considered using edgeR by Bioconductor, I decided not to choose this as it shows a varying degree of false discovery proportion between sample types as well as lower precision levels when analyzing replicates. Both limma and edgeR are available packages by Bioconductor and I utilize both of these when conducting my analysis. Sun et al. use Cufflinks to carry out their RNA-seq analysis and TopHat for alignment. My main analysis is using the function voom by limma which uses linear regression modeling  while analyzing the read counts for the RNA sequence data (Law et. al., 2014). The functions I use from edgeR are DGEList which helps me covert my data frame downloaded from GEO into a DGEList object and filterByExpression which allows the data to filter out the lowly expressed genes. I also use the package geoseq by Bioconductor to carry out my Gene Ontology analysis which accounts for the gene length bias of our differentially expressed genes. Understanding this vignette and carrying out the analysis was a very new experience for me. It was thoroughly an educational experience and helped me understand the difficulties of using new data objects. I consulted many vignettes and tutorials to be able to carry out the analysis seen in the last figure. of this script. However, this journey has encouraged me to develop new investigative skills and has me excited for my future summer project!



#### 5- MAIN ANALYSIS ----

#The main analysis of this script is showcasing the genes that are differentially expressed among the two groupings of our data, tumor and luminal. For this, we are in the assumption that our data is normally distributed. First we create a design matrix for the cell-type groupings. I used the vignette (https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html) to help me determine if I need to use an intercept term. I will not include an intercept term as the grouping are not covariates and the models are considered equivalent.
mtDesign <- model.matrix(~0+group)
colnames(mtDesign) <- gsub("group", "", colnames(mtDesign))
head(mtDesign)

#We have to set up model contrasts for the our groups as they do not have an associated intercept. The makeContrasts function from limma is used to make pairwise comparison between the two cell populations.
mtContrasts <- makeContrasts(TvsL = tumor - luminal, levels = colnames(mtDesign))
mtContrasts

#Honestly, I found going along the vignette quite difficult for my data set. Each step took a few hours to determine if my arguments are entered correctly or if I have to manipulate my data to input the correct class or even if my results make sense. As I am not completely sure if my groupings are accurate, I am always worried that my results are incorrect. However, as this assignment is about conducting an analysis along with learning the ups and downs about data sets and packages, I believe that I am learning a lot along the way.


#Removing heteroscedasticity from count data

#When using raw counts for RNA-seq data, the variance is not considered independent from the mean. For this script, we use the lcpm values where we assume our data is normally distributed. The function, voom, calculates precision weights on the mean-variance dependency by using library size and the normalization factors.
par(mfrow = c(1,2))
V_EList <- voom(DGE_Data, mtDesign, plot = F, save.plot = T)
plot(V_EList$voom.xy, 
     xlab = "Average log of counts", 
     ylab = "Sqrt (standard deviation)", 
     main = "voom: Mean-variance trend", 
     pch = 20, cex = 0.1)
lines(V_EList$voom.line, col = "red")
VFit <- lmFit(V_EList, mtDesign)
VFit <- contrasts.fit(VFit, contrasts = mtContrasts)
EFit <- eBayes(VFit)
plotSA(EFit, main = "Final Model: Mean-variance trend", 
       ylab = "Sqrt (sigma)", 
       xlab = "Average log of expression")

#In these plots, the means are plotted on the x axis and the variances are plotted on the y axis. After briefly reading the paper introducing voom (Law et al., 2014), I was unable to understand why the voom plot has a slight increase trend before decreasing. My understanding of linear models is not adequate enough to have an explanation; if this can be explained, I would really like to understand why. In my research, I have come across the downwards trend correlating to inaccurate groupings or the need to filter the data higher.The black dots in each of these plots correlates to a gene in our data set. We can see that this data set is smaller in comparison to the data set used in the vignette. Voom creates an EList object containing various information we have already seen in the DGEList object as well as additional information such as expression values and precision weights.


#Examining the number of differentially expressed genes

#A simple table can be created using the decideTests function from limma to summarize the number of significantly up- and down-regulated genes. As per usual, significance is cut off at a p-value of 0.05 by default.
dt <- decideTests(EFit)
summary(dt)
#We can see that most of the genes are not significantly different in the two groupings. This tells us that there are 28 genes that are up-regulated in the tumor samples and 26 down-regulated genes in the tumor samples relative to the luminal samples. This can mean that there only a few gene expressions that are affected in our data set when comparing tumor and luminal mammary samples of mice.

#This step allows us to extract the differentially expressed gene names. As 0 represents genes that are not differentially expressed, we are looking for the rest.
DE_Common <- which(dt[,1] != 0)
length(DE_Common) #54 genes
DE_Genes <- rownames(EFit$coefficients)[DE_Common] #These are the significantly expressed genes. Storing this vector to use for finding GO terms later on.
DE_Genes



#Displaying the differential expression results graphically

#Limma has a function called plotMD which allows for the genes to be displayed visually by mean-difference plots. This uses log-FC values from our lineal model analysis against the mean lcpm value.
par(mfrow = c(1,1))
plotMD(EFit, column = 1, status = dt[,1],
       main = "Tumor vs Luminal", 
       xlim = c(0, 13),
       xlab = "Average log of expression",
       ylab = "logFC")


#A MA plot is an interpretation of the Bland-Altman plot which calculates the agreement between quantitative measures (Giavarina, 2015). The mean and sd of the differences between the measurements are used to calculate the statistical limits (Giavarina, 2015). It is a scatter plot where the x axis displays average of the means and the y axis measures the difference of the paired measurements we calculated earlier when creating our contrast matrix (Giavarina, 2015). As we can see in the colored dots, the red ones display the up-regulated genes and the blue dots show the down-regulated genes.

#A heatmap can also be used to display the differentially expressed genes. We will be using all 54 DE genes. We can utilize the heatmap.2 function from the gplots package to create this.
dfLCPM <- cpm(DGE_Data, log = T)
Color_HeatMap <- colorpanel(1000, "turquoise4", "white", "tomato3")

heatmap.2(dfLCPM[DE_Common,], 
          scale = "row",
          labRow = DE_Genes,
          labCol = group,
          key = T,
          key.par = list(mar = c(0,1,5,1)),
          col = Color_HeatMap,
          trace = "none",
          density.info = "none",
          lhei = c(2,10),
          dendrogram = "column",
          ColSideColors = rep(c("lightgoldenrod2", "palegreen2"), each = 12))

legend(x = -4, y = 0.6, 
       title = as.expression(bquote(bold("Groups"))), 
       legend = c("tumor", "luminal"), 
       fill = c("lightgoldenrod2", "palegreen2"),
       cex = 0.7,
       box.lty = 0, 
       title.adj = 0.2) 
#Side note, it took me over two hours to get this right..why are legends so difficult?

#As we can see, this data set is quite a mess. There are two tumor samples that have clustered with the luminal samples. As I talked about in the data acquisition section, both these cell-types are from a BRCA1-deficient mouse; the paper was studying the intratumor variations so a very clear distinction of gene expression is not going to be seen. However, we can mildly see a distinction between the up-regulated and down-regulated genes in the tumor and luminal samples. Most of the genes at the top of the plot show an up-regulation in tumor cells and the most of the genes at the bottom of the plot show an up-regulation in the luminal cells. The next step for my analysis is to correlate these differentially expressed genes with their associated GO terms.


##Gene Ontology Analysis

#To perform this analysis, I have to create a named vector containing all the genes analyzed in this study and if there are differentially expressed. I will start by creating a vector containing the names of the 3382 genes we used.
Assayed_Genes <- rownames(dfLCPM)

#DE_Genes is a vector we already created to showcase the differentially expressed genes. It consists of 54 genes.

#After looking through the supportedOrganisms() list provided by goseq, I realized that the newest genome for Mus musculus (mm39) is not available through the offline database of goseq, hence I needed to install the Bioconductor package for this specific genome and all its gene annotations (TxDb.Mmusculus.UCSC.mm39.refGene). I also had to go through the vignette for Genomic Feautures (https://bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf) as it helped with extracting the data I required from the TxDb object.

#As the information I have is the gene names (SYMBOL), I have to use the Mus.musculus package to incorporate the gene names into the TxDb object. I am switching the mm10 TxDb from the Mus.musculus package to the newer mm39 TxDb.
TxDb(Mus.musculus) <- TxDb.Mmusculus.UCSC.mm39.refGene
Mus.musculus
#Can see that the TxDb object has been replaced.

#Obtaining the Gene ID's for each gene vector.
dbGeneID <- AnnotationDbi::select(Mus.musculus, keys = Assayed_Genes, columns = "GENEID", keytype = "SYMBOL")
dbDE_GENEID <- AnnotationDbi::select(Mus.musculus, keys = DE_Genes, columns = "GENEID", keytype = "SYMBOL")

#Create vectors with the GENEID for both the full list of genes and DE genes.
Assayed_Genes_ID <- as.vector(dbGeneID$GENEID)
DE_Genes_ID <- as.vector(dbDE_GENEID$GENEID)

#Now we will incorporate these two vectors to output a numerical vector stating if the original list of genes were differentially expressed (1) or not differentially expressed (0).
Gene_Vector <- as.integer(Assayed_Genes_ID %in% DE_Genes_ID)
names(Gene_Vector) <- Assayed_Genes_ID

#Let's look to see if we did what we intended.
head(Gene_Vector, 20)


#To carry out the rest of our GO analysis, we have to manually extract the gene lengths for each of these genes before fitting the probability weighting function (PWF).
dfTxDb_Gene_Length <- transcriptLengths(TxDb.Mmusculus.UCSC.mm39.refGene)
head(dfTxDb_Gene_Length, 10)

#Extract the gene lengths for each of the genes we require.
dfTxDb_Gene_Length_Subset <- dfTxDb_Gene_Length %>%
  filter(gene_id %in% Assayed_Genes_ID)
#We can see that many genes have multiple entries. We will have to average the gene lengths for the repeated entries.
dfTxDb_Gene_Length_Subset %<>%
  dplyr::select(gene_id, tx_len) %>%
  group_by(gene_id) %>%
  summarise(mean_length = ceiling(mean(tx_len)))

head(dfTxDb_Gene_Length_Subset, 10)
dim(dfTxDb_Gene_Length_Subset)

#We can see that not every gene in our Assayed_Genes have an entry in this data frame. As those genes cannot be analyzed for the rest of out script, I will delete those from Gene_Vector.
Assayed_Genes_ID <- as.vector(dfTxDb_Gene_Length_Subset$gene_id)
Gene_Vector <- as.integer(Assayed_Genes_ID %in% DE_Genes_ID)
names(Gene_Vector) <- Assayed_Genes_ID

pwf <- nullp(Gene_Vector, bias.data = dfTxDb_Gene_Length_Subset$mean_length, plot.fit = F)
summary(pwf$pwf)

#All the PWF values are lower than 0.018. PWF shows the probability of each gene being differentially expressed based on their length alone. My understanding from the plot created and the summary of values is that due to the small size of the data set and the very few DE genes, the calculation of PWF here cannot give us a lot of accurate information. When plotted, the points are so far apart and do not show any comprehensive trend for gene length. The plot somewhat shows when the length of the genes are longer, the genes are less likely to be differentially expressed (the line curves downwards, opposite from the vignette).

#Let's conduct the GO enrichment analysis and see what happens.
#I searched and tried for many hours to use use the mm39 mouse genome for this next step. I was unable to find a way to incorporate the new genome into this function. I used mm9 which is the latest mouse reference genome available on the goseq offline database. Hopefully over the next week or the new months (before my project), I can figure out how to utilize the newest reference genome.
#GO_Results_1 <- goseq(pwf = pwf, genome = "mm39", id = "GENEID", test.cats = "GO:BP")
#Error message: Couldn't grab GO categories automatically.  Please manually specify.
#As mm39 is not in the offline database of goseq, I would have needed to find the GO categories for each gene manually and create a vector to input into this function. Due to time restraints for this assignment, I decided to use the most recent genome available in the offline goseq database, mm9 to finish my analysis.

Gene_Vector2 <- as.integer(Assayed_Genes %in% DE_Genes)
names(Gene_Vector2) <- Assayed_Genes

pwf2 <- nullp(Gene_Vector2, "mm9", "geneSymbol", plot.fit = F)

GO_Results <- goseq(pwf2, "mm9", "geneSymbol", test.cats = "GO:BP")
#Turn off warning messages in RM

#Plot the top 10 GO term hits.
GO_Results %>% 
  top_n(10, wt = -over_represented_pvalue) %>% #Young et al., 2010.
  mutate(hitsPerc = numDEInCat*100/numInCat) %>% 
  ggplot(aes(x = hitsPerc, 
             y = term, 
             color = over_represented_pvalue, 
             size = numDEInCat)) +
  geom_point() +
  scale_x_continuous(breaks = c(0,10,20,30)) +
  expand_limits(x = 0) +
  labs(x ="Hits (%)", 
       y = "GO Terms", 
       color = "p-value", 
       size = "Count", 
       title = "GO Enrichment Analysis")+
  theme(aspect.ratio = 1.7/1, 
        axis.text.y = element_text(size = 10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) #Adapted argument from my first assignment


#We can see that the top 10 GO terms seen in this image showcase a high magnitude immune response. I believe we can interpret this as the differentially expressed genes seen in our data set, showcase an effect seen on the immune response. This makes sense as we are analyzing the initiation and progression of tumor formation in the mice mammary glands. I will further explain this in the results and discussion section.



#### 6- RESULTS AND DISCUSSION ----




#The question I would like evaluate for this script is, what are the biological processes of the differentially expressed genes in this data set? Are the DE genes related to the biological function of BRCA1, i.e., are they part of the individual networks for DNA repair and damage or are their functions on a broader spectrum like immune response or apoptosis? I plan to carry out gene expression analysis using the data set obtained by Sun et al. by utilizing the packages limma and goseq by Bioconductor.

#Objective met - heatmap shows up-and down-regulated genes. Research a few and explain the reason of their expression. GO enrichment analysis shows the processes significantly affected - mostly response to external stimulus and immune response - immune and defense response initiated.

#Voom was carried out on expression comparison of tumor to luminal - significance of these responses makes sense as the body of the mice is getting ready to fight, knows something is wrong.






#Compare GO to fig 3E
#Functional profiling of enrichment analysis in the last figure displays the biological process association seen in the differentially expressed genes. According to the official Gene Ontology website, p-value in enrichment analysis calculates the probability that the GO term is seen in the total list of genes inputted. The smaller the p-value, the closer the GO term is associated with the total genes. In our case, the darker the blue, the more significant is the term association in comparison to the whole list of GO terms.


#Reflection paragraph


#### 7- ACKNOWLEDGEMENTS ----








#### 8- REFERENCES ----

#https://master.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
#RNASeq Bioconductor vignette

#https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html
#Design matrices vignette

#https://bioconductor.org/packages/devel/bioc/vignettes/goseq/inst/doc/goseq.pdf
#goseq vigenette

#https://bioconductor.org/packages/release/data/annotation/manuals/TxDb.Mmusculus.UCSC.mm39.refGene/man/TxDb.Mmusculus.UCSC.mm39.refGene.pdf
#mm39 genome vignette

#https://bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf
#To use TxDB object

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148569
#Data

#Sun H, Zeng J, Miao Z, Lei KC, Huang C, Hu L, Su SM, Chan UI, Miao K, Zhang X, Zhang A, Guo S, Chen S, Meng Y, Deng M, Hao W, Lei H, Lin Y, Yang Z, Tang D, Wong KH, Zhang XD, Xu X, Deng CX. Dissecting the heterogeneity and tumorigenesis of BRCA1 deficient mammary tumors via single cell RNA sequencing. Theranostics 2021; 11(20):9967-9987.
#Paper of data

#Hayashizaki Y. (2003). RIKEN mouse genome encyclopedia. Mechanisms of ageing and development, 124(1), 93–102.

#Evans, C., Hardin, J., & Stoebel, D. M. (2018). Selecting between-sample RNA-Seq normalization methods from the perspective of their assumptions. Briefings in bioinformatics, 19(5), 776–792.

#Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome biology, 15(2), R29.

#Giavarina D. (2015). Understanding Bland Altman analysis. Biochemia medica, 25(2), 141–151.

#Yoshida K, Miki Y. Role of BRCA1 and BRCA2 as regulators of DNA repair, transcription, and cell cycle in response to DNA damage. Cancer Sci. 2004 Nov;95(11):866-71.

#Jhanwar-Uniyal M. BRCA1 in cancer, cell cycle and genomic stability. Front Biosci. 2003 Sep 1;8:s1107-17.

#Young, Matthew D, Matthew J Wakefield, Gordon K Smyth, and Alicia Oshlack. 2010. “Gene Ontology Analysis for Rna-Seq: Accounting for Selection Bias.” Genome Biology 11: R14.








### QUESTIONS TO ASK
