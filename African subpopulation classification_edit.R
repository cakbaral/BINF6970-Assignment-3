library("VariantAnnotation")
library("tidyverse")
library("ggplot2")
library("ggfortify")
#Load data in the file
APOBEC3G <- readVcf("gene vcf files/APOBEC3G.vcf.gz")
LCT <- readVcf("gene vcf files/LCT.vcf.gz")
Dock3 <- readVcf("gene vcf files/DOCK3.vcf.gz")

#Data Exploration
APOBEC3G
#dim(vcf)
samples(header(APOBEC3G))
seqlevels(rowRanges(APOBEC3G))
rowRanges(APOBEC3G)
#Connvert to data frame
Df_APOBEC3G <-  as.data.frame(geno(APOBEC3G)$GT)
DF_LCT <- as.data.frame(geno(LCT)$GT)
DF_Dock3 <- as.data.frame(geno(Dock3)$GT)

#Function to convert genotype information to numeric allele dosage
#Homozygous major allele 0, Heterozygous 1, Homozygous minor allele 2.
convert_gt <- function(x) {
  x <- as.character(x)
  x[x %in% c("0|0", "0/0")] <- "0"
  x[x %in% c("0|1", "1|0", "0/1", "1/0")] <- "1"
  x[x %in% c("1|1", "1/1")] <- "2"
  x[x %in% c("./.", ".|.")] <- NA
  as.numeric(x)
}
#Applying function 
Df_APOBEC3G[] <- as.data.frame(lapply(Df_APOBEC3G, convert_gt))
DF_LCT[] <- as.data.frame(lapply(DF_LCT, convert_gt))
DF_Dock3[] <- as.data.frame(lapply(DF_Dock3, convert_gt))

#Adding metadata 
metadata <- read.table(
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",
  header = TRUE,
  stringsAsFactors = FALSE
)

head(metadata)

#Transpose SNP data to have sample IDs as rownames
Df_APOBEC3G <- as.data.frame(t(Df_APOBEC3G))
DF_LCT <- as.data.frame(t(DF_LCT))
DF_Dock3 <- as.data.frame(t(DF_Dock3))

#Combine all data frames
Genotype_df <-  cbind(metadata,Df_APOBEC3G,DF_LCT,DF_Dock3)

#Working with the African superpopulation (AFR)
#We intend to use SNP data from these Ancestry informative genes to classify the AFR population into its sub-populations
AFR_DF <- Genotype_df %>% filter(super_pop == "AFR")


#Data exploration
ggplot(AFR_DF, aes(x = pop)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Subpopulations within African dataset", x = "Subpopulation", y = "Count")

#Remove SNPs with NA values
Cleaned_AFR_DF <- AFR_DF[, colSums(is.na(AFR_DF)) == 0, drop = FALSE]

#PCA
# Remove columns(SNPs) with zero variance
X <- Cleaned_AFR_DF[, -(1:4)]

X_filtered <- X[, apply(X, 2, var) != 0]

# Run PCA
PCA <- prcomp(X_filtered, center = TRUE, scale. = TRUE)
#PCA <- prcomp(Cleaned_AFR_DF[ ,-(1:4)], center = T, scale. = T)
var_explained <- (PCA$sdev^2) / sum(PCA$sdev^2)
percent_var <- var_explained * 100

percent_var
#Plot
autoplot( PCA , data = Cleaned_AFR_DF, colour = "pop", main = "PCA: PC1 vs PC2" )
0

#Test output for downstream analysis
X_filt_head <- X_filtered[1:10, ]
Y_head <- Cleaned_AFR_DF[1:10, ]

