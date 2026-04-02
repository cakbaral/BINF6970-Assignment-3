library(GGally)
library(cluster)
library(bpca)
library(pheatmap)
library(factoextra)
library(NbClust)
library(bpca)
########## Preliminary Analysis

#ggpairs(dat[,-c(1:2)])

# Biplot: use prcomp(data, scale=T) in biplot() to get a biplot
# Scree plot: plot(prcomp(data,scale=T)) plots the sdev^2 (var) explained by all the PC's

########## Hierarchical Clustering

# Clustering by var: use cor(data) via agnes (); library(cluster)
# Clustering by brand: use dist(data) via agnes()
# Agglomerative: use agnes()
# Divisive: use diana()

##########  Non-Hierarchical Clustering
# k-means: kmeans(data,3) # 3-means 
# k-medoids: pam(dat.pu,3) # 


##=============================
## Prepare the data
##=============================

dat <- X_filtered

ggpairs(dat)

raw.brand <- as.matrix(dat) ## extract just the variables
rownames(raw.brand) <- 1:dim(dat)[sample]

## data with manufacturer row labels
raw.manu <- raw.brand; rownames(raw.manu) = as.character(dat[,2])

scale.brand <- scale(raw.brand) # standardized data with 1:43 as brand labels
scale.manu <- scale(raw.manu) # standardized data with G,K,Q as manufacturer labels 


## Visualize the "dis-similarity" matrices:

## For columns (variables)
fviz_dist(as.dist(cor(raw.brand)), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

## For rows (brands)
fviz_dist(dist(raw.brand), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

### Although not performed here, you should repeat the following analysis 
### without standardized data as well!


##=============================
######### 2D Biplot
##=============================

par(mfrow=c(1,2))
biplot(prcomp(scale.brand,scale=F),main= "Biplot by Brand")


#Scree Plot
barplot(prcomp(scale.brand,scale=F)[[1]]^2,
        main= "Scree Plot", names.arg=1:8, xlab="PCs", ylab="variance")
par(mfrow=c(1,1))

##=============================
######### 3D Biplot
##=============================

plot(bpca(as.data.frame(scale.brand),
          method='hj',
          d=1:3),
     rgl.use=TRUE,
     var.col='brown',
     var.factor=.3,
     var.cex=1.2,
     obj.names=FALSE,
     obj.cex=.8,
     obj.col=c('red', 'green3', 'orange')[dat[,2]],
     simple.axes=FALSE,
     box=TRUE)

##==================================
######### Agglomerative Clustering
##==================================

## Cluster the variables (various linkage types produce similar results)

?dist  ## compute distance matrix - you can try "manhattan distance" as well

plot(agnes(abs(cor(scale.manu)), method = "complete"), main="using agnes()",which.plots=2)

## Various Agglomerative plots

par(mfrow=c(2,2))
plot(agnes(dist(scale.brand), method = "single"), main="using agnes(brand-single)",which.plots=2)
plot(agnes(dist(scale.brand), method = "complete"), main="using agnes(brand-complete)",which.plots=2)
plot(agnes(dist(scale.brand), method = "average"), main="using agnes(brand-average)",which.plots=2)

### Divisive Clustering

plot(diana(dist(scale.brand)),which.plots=2, main="using diana(brand)")
par(mfrow=c(1,1))


### Identify the clusters and examine nutritional profiles

plot(hclust(dist(scale.brand), method = "complete"))
abline(h=6)
abline(h=5.5,col='blue')
abline(h=4.7,col='red')


cut_complete <- cutree(hclust(dist(scale.brand), method = "complete"), h=4.7)
cut_complete
table(cut_complete)

## get cluster memberships

clusters <- lapply(1:length(table(cut_complete)), 
                   function(i) names(cut_complete)[which(cut_complete==i)])

clusters

dat2 <- dat
dat2$cutcomp <- cut_complete

aggregate(.~cutcomp,dat2[,-c(1:2)],mean)


## Can you identify the healthy and unhealthy clusters?

### HEATMAP

heatmap(scale.brand, scale = "none")  ## base R
pheatmap(scale.brand,cutree_cols = 3) ## from library(pheatmap)