#################################################
# EXERCISE 1

#1) Input the Dataset 2

#2) Define a sequence object with elements in data columns 2:61
# and alphabet 1:6, using the following state names and labels
#    1 SNP "Single, childless",
#    2 SBP "Single, child b/separat.",
#    3 SAP "Single, child a/separat.",
#    4 UNP "Union, childless",
#    5 UBP "Union, child b/separat.",
#    6 UAP "Union, child a/separat."

#3) Display (print) the first 10 sequences in extended and compact form

#4) plot a full representation of sequences, and order them from the first state

#5) plot the 5 most frequent sequences. Comment the plot

#6) Create a state distribution plot for each birthcohort (BIRTHCOH)
# What are the cross-cohort differences in the distribution of states overtime? 

#7) What are the most frequent states one and five years after break-up?
# Use a modal state plot for illustration.

#8) Assess the cross-sectional state diversity plotting a measure of entropy
# At what time after separation is the cross-sectional diversity of the states at its highest?

#9) Display side by side in a same plot area the mean times spent 
#in each of the states and the sequence of modal states.

#10) Compute the (overall) transition rate matrix 
# What is the largest transition rate between two different states?

#11) compute the sequence length, the number of transitions, 
#the number of subsequences and the longitudinal entropy

#12) Using summary(), look at the min, max, mean, median and quartiles 
#of the distribution of each of the computed longitudinal characteristics.


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
                     2:26, 	
                     from="STS", 
                     to="SPS")


