#################################################
# EXERCISE 2

#1) Input the Dataset 2

#2) Define a sequence object with elements in data columns 2:61
# and alphabet 1:6, using the following state names and labels
#    1 SNP "Single, childless",
#    2 SBP "Single, child b/separat.",
#    3 SAP "Single, child a/separat.",
#    4 UNP "Union, childless",
#    5 UBP "Union, child b/separat.",
#    6 UAP "Union, child a/separat."

#3) Compute the matrix of pairwise OM with constant cost distances between all 
# sequences and display the results for the first 5 sequences.

#4) Plot the first 2 sequences and check that the OM with constant cost distance 
# is the number of non matching positions between them.

#5) Check data that the LCS distance provides the same (non-normalized) 
# distances as OM with indel=1 and a constant substitution cost of 2 

#6) Define a substitution cost matrix reflecting what (according to your prior knowledge)
# are the distances between two states (i.e. customize state-dependent substitution costs)

#7) Compute the OM dissimilarity matrix using the previously derived substitution
# Set the indel cost as half the maximum substitution cost.

#8)  From the previously computed OM dissimilarity matrix, create a 
# hierarchical cluster tree object with Ward method. Display the hierarchical tree

#9) Calculate appropriate cluster cut-off criteria. 
# Assess what is an empirically optimal cluster solution

#10) Select the six-cluster solution from the Ward analysis, 
# check cluster consistency, and label the clusters by looking at the full sequence index plots 
# (or the relative frequency version) by cluster.

#11) Repeat steps 8-10 using a DHD dissimilarity matrix

#12) Compare the results between the OM and the DHD approaches




# Importing libraries

library(TraMineR)
library(ggplot2)
library(grDevices)
library(graphics)
library(foreign)
library(cluster)
library(Hmisc)
library(TraMineRextras)
library(WeightedCluster)
library(RColorBrewer)
library(colorspace)

# Importing data
data2 <- read.csv("SFS2018_Data2.csv", na.strings=c(".",".a",".b"))



#################################################
### Conversion across longitudinal data formats

#generate new data object where sequence element are in SPS format
data2SPS<-seqformat( data2,	
                     2:61, 	
                     from="STS", 
                     to="SPS")

# X-axis for exercise
xtlab <- seq(1,60, by=1)

#vector for the state labels
seqlab <-c("Single, childless", "Single, child b/separat.", 
           "Single, child a/separat.", "Union, childless", 
           "Union, child b/separat.", "Union, child a/separat.")

#vector of short state names (default would be alphabet labels)
sllist <- c("SNP","SBP","SAP","UNP", "UBP", "UAP")

#################################################
###  Generate sequence object
seqObj2 <- seqdef(data2, 
                  var=2:61,  
                  alphabet=c(1:6),
                  cpal=color1, 
                  states=sllist, 
                  labels=seqlab)


# 3) Compute the matrix of pairwise OM with constant cost distances between all 
# sequences and display the results for the first 5 sequences.

#OM with CONSTANT subcosts (OM with indel=1, subs=2)
Matrix.OM.Const <- seqdist(seqObj2,method="OM", indel=1, sm="CONSTANT")
#display matrix
print(Matrix.OM.Const[1:5,])



# 5) Check data that the LCS distance provides the same (non-normalized) 
# distances as OM with indel=1 and a constant substitution cost of 2 

#Longest common subsequence
Matrix.LCS <- seqdist(seqObj2,method="LCS")


# Comparing: 
Matrix.OM.Const - Matrix.LCS


# 6) Define a substitution cost matrix reflecting what (according to your prior knowledge)
# are the distances between two states (i.e. customize state-dependent substitution costs)


#OM with customized state-dependent subcosts 
submatrix <- matrix(c( 0,6,2,4,2,4,
                       6,0,6,4,6,4,
                       2,6,0,2,2,3,
                       4,4,2,0,4,2,
                       2,6,2,4,0,2,
                       4,4,3,2,2,0), nrow = 6, ncol = 6, byrow = TRUE)


# 7) Compute the OM dissimilarity matrix using the previously derived substitution
# Set the indel cost as half the maximum substitution cost.

Matrix.OM.State.dep <- seqdist(seqObj2,method="OM", indel=3, sm=submatrix)


# 8)  From the previously computed OM dissimilarity matrix, create a 
# hierarchical cluster tree object with Ward method. Display the hierarchical tree

ward.OM <- hclust(as.dist(Matrix.OM.State.dep), method = "ward.D2")

plot(ward.OM)


# 9) Calculate appropriate cluster cut-off criteria. 
# Assess what is an empirically optimal cluster solution

### Generate an object with 1-10 cluster solutions for each prior anal
wardrange.OM <-as.clustrange(ward.OM, diss=Matrix.OM.State.dep, ncluster=10)
summary(wardrange.OM, max.rank=3)
plot(wardrange.OM, stat=c("ASW", "HGSD", "PBC"), norm="zscore")

# 10) Select the six-cluster solution from the Ward analysis, 
# Check cluster consistency, and label the clusters by looking at the full sequence index plots 
# (or the relative frequency version) by cluster.

### store cluster solutions with best empirical fits
wardrange.OM.6 <- cutree(ward.OM , k=6)

# OM 6-cluster solution
silh.OM.6 <- silhouette(wardrange.OM.6, dmatrix = Matrix.OM.State.dep)
summary(silh.OM.6)

# USE WesAnderson palette!!!
plot(silh.OM.6, main= "Silhouette - OM 6 cluster", col=c("#E2E2E2", "#D3D3D3" ,"#B8B8B8", "#969696", "#707070", "#E2E2E2"), border=NA)


# OM 6-cluster solution
seqIplot(seqObj, group = wardrange.OM.6,use.layout=FALSE, border=NA,xtlab=xtlab, space=0,yaxis=FALSE)


