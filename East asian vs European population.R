library("VariantAnnotation")
library("tidyverse")
library("ggplot2")
library("ggfortify")
#Load data in the file
SLC24A5<- readVcf("gene vcf files/SLC24A5.vcf.gz")
EDAR <- readVcf("gene vcf files/EDAR.vcf.gz")
#Data Exploration
SLC24A5
samples(header(SLC24A5))
seqlevels(rowRanges(SLC24A5))
rowRanges(SLC24A5)
#Connvert to data frame
DF_SLC24A5 <-  as.data.frame(geno(SLC24A5)$GT)
DF_EDAR <- as.data.frame(geno(EDAR)$GT)

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
DF_SLC24A5[] <- as.data.frame(lapply(DF_SLC24A5, convert_gt))
DF_EDAR[] <- as.data.frame(lapply(DF_EDAR, convert_gt))

#Adding metadata 
metadata <- read.table(
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",
  header = TRUE,
  stringsAsFactors = FALSE
)

head(metadata)

#Transpose SNP data to have sample IDs as rownames
DF_SLC24A5 <- as.data.frame(t(DF_SLC24A5))
DF_EDAR <- as.data.frame(t(DF_EDAR))


#Combine all data frames
Genotype_df <-  cbind(metadata,DF_SLC24A5,DF_EDAR)

#Working with the European  an East asian superpopulation (AFR)
#We intend to use SNP data from these Ancestry informative genes to classify the population.
COMBINED_DF <- Genotype_df %>%
  filter(super_pop %in% c("EAS", "EUR"))

#Data exploration
ggplot(COMBINED_DF, aes(x = super_pop)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Populations within dataset", x = "Superpopulation", y = "Count")

#Remove SNPs with NA values
Cleaned_DF <- COMBINED_DF[, colSums(is.na(COMBINED_DF)) == 0, drop = FALSE]


# Remove columns(SNPs) with zero variance
X <- Cleaned_DF[, -(1:4)]

X_filtered <- X[, apply(X, 2, var) != 0]

#Filtering low minor allele frequency below 0.01
maf <- apply(X_filtered, 2, function(snp) {
  p <- sum(snp, na.rm = TRUE) / (2 * sum(!is.na(snp)))
  min(p, 1 - p)
})

X_filtered <- X_filtered[, maf > 0.01]


# Run PCA
PCA <- prcomp(X_filtered, center = TRUE)
var_explained <- (PCA$sdev^2) / sum(PCA$sdev^2)
percent_var <- var_explained * 100

percent_var
#Combine data with metadata for plot
metadata_sub <-  metadata %>%
  filter(super_pop %in% c("EAS", "EUR"))

X_filtered2 <- cbind(metadata_sub, X_filtered)

#Plot
autoplot( PCA , data = X_filtered2t , colour = "super_pop", main = "PCA: PC1 vs PC2" )


