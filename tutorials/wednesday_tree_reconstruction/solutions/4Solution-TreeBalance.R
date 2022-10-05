#-------------------------------------------------------------------------------
################################################################################
#                             ViroInf Workshop
################################################################################
#-------------------------------------------------------------------------------
# Tutorial on various tree reconstruction methods using R
# by Sophie Kersting
# Wednesday, October 5, 2022
################################################################################
#-------------------------------------------------------------------------------
# 4) Brief introduction to tree balance (indices)
#-------------------------------------------------------------------------------
################################################################################
# Load the required package:
library("treebalance")
################################################################################
# Edit the example path to set your working directory.
setwd("C:\\Users\\Sophie\\Desktop\\ViroInfTutorial")
# ----- load the R object
load(file="trees_ebov_cds_rooted.RData")
load(file="trees_ebov_ig_rooted.RData")
# ----- or load the Newick format
trees_ebov_cds_rooted <- read.tree(file="trees_ebov_cds_rooted.txt")
trees_ebov_ig_rooted <- read.tree(file="trees_ebov_ig_rooted.txt")
#_______________________________________________________________________________
# Task: Analyze how well the Yule process describes the evolutionary history
# of the samples in the ebov dataset for the coding and non-coding regions.
#+++++++++++++++++++++++++++++++++++
# Here are a few guiding steps:
# 1.) Generate N=2,000 Yule trees with the same amount of leaves.
Yule_trees <- lapply(1:2000, function(x){
  ape::rtree(n = length(trees_ebov_cds_rooted[[1]]$tip.label))})
# 2.) Compute the Colless index values of the Yule trees and the reconstructed
#     trees.
Yule_values <- sapply(1:2000, function(x){
  collessI(Yule_trees[[x]])})
cds_values <- sapply(1:length(trees_ebov_cds_rooted), function(x){
  collessI(trees_ebov_cds_rooted[[x]])})
ig_values <- sapply(1:length(trees_ebov_ig_rooted), function(x){
  collessI(trees_ebov_ig_rooted[[x]])})
# 3.) Compute the quantiles and visualize everything with a bar plot. 
#     Adjust the plot to your liking. Add vertical lines (with text) for the 
#     different reconstructed trees.
Yule_quantiles <- quantile(Yule_values,  probs = c(0.025,0.975))
hist(Yule_values, xlim = c(min(Yule_values,cds_values,ig_values)-2,
                           max(Yule_values,cds_values,ig_values)+2),
     main = "Hypothesis test for ebov dataset based on the Colless index",
     xlab = "Colless index")
abline(v = Yule_quantiles, lty=1)
abline(v = cds_values, lty=2, lwd=3, col= "blue")
abline(v = ig_values, lty=3, lwd=3, col= "red")
legend("topright", legend=c("quantiles", "cds trees", "ig trees"),
       col=c("black", "blue", "red"), lty = c(1,2,3), lwd=c(1,3,3), cex=1.2)
# 4.) Interpret the result.
#_______________________________________________________________________________
# Optional task: Try other indices, e.g. the Sackin or total cophenetic index
# (collessI or totCophI). See help pages for index definition.
#_______________________________________________________________________________
Yule_values <- sapply(1:2000, function(x){
  sackinI(Yule_trees[[x]])})
cds_values <- sapply(1:length(trees_ebov_cds_rooted), function(x){
  sackinI(trees_ebov_cds_rooted[[x]])})
ig_values <- sapply(1:length(trees_ebov_ig_rooted), function(x){
  sackinI(trees_ebov_ig_rooted[[x]])})
Yule_quantiles <- quantile(Yule_values,  probs = c(0.025,0.975))
hist(Yule_values, xlim = c(min(Yule_values,cds_values,ig_values)-2,
                           max(Yule_values,cds_values,ig_values)+2),
     main = "Hypothesis test for ebov dataset based on the Sackin index",
     xlab = "Sackin index")
abline(v = Yule_quantiles, lty=1)
abline(v = cds_values, lty=2, lwd=3, col= "blue")
abline(v = ig_values, lty=3, lwd=3, col= "red")
legend("topright", legend=c("quantiles", "cds trees", "ig trees"),
       col=c("black", "blue", "red"), lty = c(1,2,3), lwd=c(1,3,3), cex=1.2)

