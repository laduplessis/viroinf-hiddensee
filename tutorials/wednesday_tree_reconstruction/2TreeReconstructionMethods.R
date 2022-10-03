#-------------------------------------------------------------------------------
################################################################################
#                             ViroInf Workshop
################################################################################
#-------------------------------------------------------------------------------
# Tutorial on various tree reconstruction methods using R
# by Sophie Kersting
# Wednesday, 5th October, 2022
################################################################################
#-------------------------------------------------------------------------------
# 2) Tree reconstruction methods
#-------------------------------------------------------------------------------
################################################################################
# Load the required packages:
library("phangorn")
################################################################################
# a) The dataset: 15 Ebola virus (EBOV) genomes -> Adapt the paths!
#   CDS : Coding regions of the 15 EBOV genomes (14517 bp)
align_ebov_cds <- read.phyDat(paste0("C:\\Users\\Sophie\\Documents\\GitHub\\",
                                 "viroinf-hiddensee\\datasets\\ebov\\",
                                 "EBOV_reference_set_15_cds.nexus"), 
                              format = "nexus")
#   IG : Non-coding regions of the 15 EBOV genomes (4449 bp)
align_ebov_ig <- read.phyDat(paste0("C:\\Users\\Sophie\\Documents\\GitHub\\",
                                 "viroinf-hiddensee\\datasets\\ebov\\",
                                 "EBOV_reference_set_15_ig.nexus"), 
                             format = "nexus")
# ----- What do the datasets look like?
align_ebov_cds
align_ebov_ig 
# ----- Convert to DNAbin to have more models (JC69, K2P, ...) available.
dna_ebov_cds <- as.DNAbin(align_ebov_cds)
as.character(dna_ebov_cds[1:4,1:10]) # t instead of u
#_______________________________________________________________________________
# Task: Explore this data structure. Then:
# 1.) Convert the non-coding regions to DNAbin and print the first and last
# five characters of the first 6 sequences.
# 2.) Do both datasets observe the same samples? Use the labels() function on
# the alignments (or DNAbins) and compare them.
#_______________________________________________________________________________

################################################################################
# b) Neighbour Joining - reconstruct trees based on distance-matrices
# ----- create lists to save all future reconstructed trees
trees_ebov_cds <- list()
attr(trees_ebov_cds, "class") <- "multiPhylo"
trees_ebov_ig <- list()
attr(trees_ebov_ig, "class") <- "multiPhylo"
# ----- create the distance matrix (choose evolutionary model)
dm_ebov_cds_JC  <- dist.dna(dna_ebov_cds, model = "JC69")
trees_ebov_cds$treeNJ_JC <- NJ(dm_ebov_cds_JC)
ape::plot.phylo(trees_ebov_cds$treeNJ_JC, type = "unrooted", 
                cex=0.8, lab4ut = "axial")
#_______________________________________________________________________________
# Task: Experiment with NJ. For example, try:
# 1.) Create the NJ tree for the non-coding regions.
trees_ebov_ig$treeNJ_JC <- NULL # replace NULL with your code
# 2.) Create for both cds and ig the NJ trees with a different model "K80"
#     Kimura's 2-parameter model
trees_ebov_cds$treeNJ_K80 <- NULL # replace NULL with your code
trees_ebov_ig$treeNJ_K80 <- NULL # replace NULL with your code
# 3.) Plot the 4 trees and find similarities and differences:
par(mfrow=c(2,2),mar=c(1.1,2,2,2.8)) # plot in a 2x2 grid
# insert all 4 plot commands in here
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1)) # reset to normal
#_______________________________________________________________________________

################################################################################
# c) Robinson-Foulds distance - How different are these reconstructed trees?
# ----- enter two single trees
RF.dist(trees_ebov_cds$treeNJ_JC, tree2=trees_ebov_ig$treeNJ_JC, 
        check.labels=TRUE, rooted=FALSE)
# ----- or a list of trees
RF.dist(trees_ebov_cds, # cds trees
        check.labels=TRUE, rooted=FALSE)
RF.dist(append(trees_ebov_cds,trees_ebov_ig), # full list
        check.labels=TRUE, rooted=FALSE)
#_______________________________________________________________________________
# Task: For now, the trees are still unrooted. You can check this with
is.rooted.phylo(trees_ebov_cds$treeNJ_JC)
# However, we actually know where the root of these trees should be positioned:
# The two samples from the 1970s (Yambuku-Mayinga and Bonduni) are the outgroup.
# 1.) Edit the command below to correctly root the unrooted trees and plot
#     the rooted version.
trees_ebov_cds_rooted <- list()
attr(trees_ebov_cds_rooted, "class") <- "multiPhylo"
trees_ebov_cds_rooted$treeNJ_JC_r <- root.phylo(trees_ebov_cds$treeNJ_JC, 
                                outgroup = c("leaf label 1",
                                             "leaf label 2"), 
                                resolve.root = T)
trees_ebov_cds$treeNJ_JC$tip.label
# 2.) Similarly root the other three trees.
# 3.) Calculate the RF distances between the rooted trees.
#_______________________________________________________________________________


################################################################################
# d) Maximum parsimony - reconstruct trees based on the alignment
# ----- Compute the maximum parsimony score for a single given tree:
parsimony(trees_ebov_cds$treeNJ_JC, data=align_ebov_cds, method = "fitch")
# ----- Now, find the "optimal" MP tree (use NJ_JC tree as a starting point)
trees_ebov_cds$treeMP <- optim.parsimony(tree=trees_ebov_cds$treeNJ_JC, 
                                         data=align_ebov_cds, method = "fitch", 
                                         trace = 10, rearrangements = "SPR")
parsimony(trees_ebov_cds$treeMP, data=align_ebov_cds, method = "fitch")
par(mfrow=c(1,2),mar=c(1.1,2,2,2.8))
ape::plot.phylo(trees_ebov_cds$treeNJ_JC, main="CDS NJ tree (JC69)", 
                type = "unrooted", cex=0.6, lab4ut = "axial")
ape::plot.phylo(trees_ebov_cds$treeMP, main="CDS MP tree", 
                type = "unrooted", cex=0.6, lab4ut = "axial")
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
RF.dist(trees_ebov_cds$treeNJ_JC, tree2=trees_ebov_cds$treeMP, 
        check.labels=TRUE, rooted=FALSE) # only edge lengths different
# ----- An improved algorithm to find the "optimal" MP tree is parsimony ratchet
#       (Nixon 1999). We will not go into details here
trees_ebov_cds$treeMPratched <- pratchet(align_ebov_cds, 
                                         start=trees_ebov_cds$treeNJ_JC,
                                         trace = 0, minit=100)
parsimony(trees_ebov_cds$treeMPratched, align_ebov_cds) # no improvement
#_______________________________________________________________________________
# Task: Explore maximum parsimony. For example, try:
# 1.) Compute the maximum parsimony scores for the other 3 distance-based trees.
#     Can you explain the differences and similarities?
# 2.) Compute the the MP tree for the non-coding sequences and compare the two
#     MP trees for cds and ig and the respective NJ trees (JC).
#_______________________________________________________________________________


################################################################################
# e) Maximum likelihood - reconstruct trees based on the alignment
# ----- We can compute the likelihood for a given tree:
fitJC <- pml(trees_ebov_cds$treeNJ_JC, data=align_ebov_cds,
           model = "JC")
fitJC
# ----- Then we can optimize the tree for the given model
fitJC_opt <- optim.pml(fitJC, rearrangement="NNI")
trees_ebov_cds$treeML_JC  <- fitJC_opt$tree
# ----- We can also use other models: Like Kimura's 2-parameter model
fitK80 <- pml(trees_ebov_cds$treeNJ_JC, data=align_ebov_cds,
           model = "K80")
fitK80_opt <- optim.pml(fitK80, rearrangement="NNI")
trees_ebov_cds$treeML_K80  <- fitK80_opt$tree
# ----- or GTR+G(4)+I
fitGTRGI <- pml(trees_ebov_cds$treeNJ_JC, data=align_ebov_cds, k=4,
            model = "GTR")
fitGTRGI_opt <- optim.pml(fitGTRGI, model="GTR", optInv=TRUE, 
                          optGamma=TRUE, k=4)
trees_ebov_cds$treeML_GTRGI  <- fitGTRGI_opt$tree
# ----- Use the following command to get an overview over different models 
# -> How to find the best model (we will not go into detail here)
mt <- modelTest(align_ebov_cds, model=c("JC", "F81", "K80","HKY", "SYM", "GTR"), 
                control = pml.control(trace = 0))
mt
#_______________________________________________________________________________
# Task: Explore maximum likelihood. For example, try:
# 1.) Compare the likelihoods of these three ML trees.
# 2.) Construct all three ML trees for the non-coding sequences.
# 3.) Compare all ML trees using the RF distances.
# 4.) Root two trees with the biggest distance as before and plot them to 
#     observe the differences.
#_______________________________________________________________________________


################################################################################
# f) Save the rooted trees for later use in Exercise 4.
# ---- root the ig ML trees
trees_ebov_ig_rooted$treeML_K80 <- root.phylo(trees_ebov_ig$treeML_K80,  
                                    outgroup = c("KC242791_Bonduni_DRC_1977_06",
                                    "KR063671_Yambuku_Mayinga_DRC_1976_10_01"), 
                                    resolve.root = T)
trees_ebov_ig_rooted$treeML_GTRGI <- root.phylo(trees_ebov_ig$treeML_GTRGI,  
                                    outgroup = c("KC242791_Bonduni_DRC_1977_06",
                                    "KR063671_Yambuku_Mayinga_DRC_1976_10_01"), 
                                    resolve.root = T)
setwd("C:\\Users\\Sophie\\Desktop\\ViroInfTutorial") # <- edit path
# ----- save the R object (to load use load(file="trees_ebov_cds_rooted.RData"))
save(trees_ebov_cds_rooted, file="trees_ebov_cds_rooted.RData")
save(trees_ebov_ig_rooted, file="trees_ebov_ig_rooted.RData")
# ----- or save in Newick format
write.tree(phy=trees_ebov_cds_rooted, file="trees_ebov_cds_rooted.txt",
           tree.names = TRUE)
write.tree(phy=trees_ebov_ig_rooted, file="trees_ebov_ig_rooted.txt",
           tree.names = TRUE)
#_______________________________________________________________________________

