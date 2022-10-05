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
# 3) Evidence for a criminal case
#-------------------------------------------------------------------------------
################################################################################
# Load the required package:
library("phangorn")
#_______________________________________________________________________________
# Task: Your main task is to find the evidence for the criminal case.
# Here are some potential steps you might want to 
# 1.) Load the HIV-1 pol and env datasets ("hiv_env.fasta" and "hiv_pol.fasta").
# Alignment HIV env
align_hiv_env <- read.phyDat(paste0("C:\\Users\\Sophie\\Documents\\GitHub\\",
                                    "viroinf-hiddensee\\datasets\\hiv\\",
                                    "hiv_env.fasta"), format = "fasta")
align_hiv_env
dna_hiv_env <- as.DNAbin(align_hiv_env)

# Alignment HIV pol
align_hiv_pol <- read.phyDat(paste0("C:\\Users\\Sophie\\Documents\\GitHub\\",
                                    "viroinf-hiddensee\\datasets\\hiv\\",
                                    "hiv_pol.fasta"), format = "fasta")
align_hiv_pol
dna_hiv_pol <- as.DNAbin(align_hiv_pol)
# 3.) For pol choose a first tree reconstruction method, e.g. NJ 
#     i.e. distance-based) and reconstruct trees.
dm_hiv_pol_JC  <- dist.dna(dna_hiv_pol, model = "JC69")
treeNJ_hiv_pol<- NJ(dm_hiv_pol_JC)
# 4.) Have a look at the labels for pol. 
#     All labels with 'V' are from the victim and with the ones with 'P' are 
#     from the patient ('LA' is data from environment)
labels(align_hiv_pol)
#     Assume that you know that the outgroup is 'AY156771.1|LA02.RT ' and root
#     the tree accordingly.
treeNJ_hiv_pol_r <- root.phylo(treeNJ_hiv_pol,  
                              outgroup = c("AY156771.1|LA02.RT "), 
                              resolve.root = T)
# 5.) Plot the rooted tree and look for the evidence. Color the patient tips in
#     blue and the victim tips in red, you can use your myColsStrings function
#     from 1IntroductionToR.R.
plot.phylo(treeNJ_hiv_pol_r, cex=0.6)
plot.phylo(treeNJ_hiv_pol_r,  cex=0.6,
           tip.color = myColsStrings(treeNJ_hiv_pol_r$tip.label, c("V","P"),
                                     c("red","blue")),
           main="NJ tree for HIV pol with colored tips")
legend("bottomleft", legend=c("Patient", "Victim", "Environment"),
       col=c("blue", "red", "black"), pch=16, cex=0.8)
# 6.) Check one of the other methods (MP or ML) similarly. Do they yield the
# same result?
# ----- Maximum Parsimony:
parsimony(tree=treeNJ_hiv_pol, data=align_hiv_pol)
treeMP_hiv_pol <- optim.parsimony(tree=treeNJ_hiv_pol, 
                                  data=align_hiv_pol, method = "fitch", 
                                  trace = 10, rearrangements = "SPR")
parsimony(treeMP_hiv_pol, data=align_hiv_pol)
treeMP_hiv_pol_r <- root.phylo(treeMP_hiv_pol,  
                               outgroup = c("AY156771.1|LA02.RT "), 
                               resolve.root = T)

plot.phylo(treeMP_hiv_pol_r, cex=0.6,
           tip.color = myColsStrings(treeMP_hiv_pol_r$tip.label, c("V","P"),
                                     c("red","blue")), 
           main="MP tree for HIV pol with colored tips")
legend("bottomleft", legend=c("Patient", "Victim", "Environment"),
       col=c("blue", "red", "black"), pch=16, cex=0.8)
# ----- Maximum Likelihood:
fitGTRGI <- pml(treeNJ_hiv_pol, data=align_hiv_pol, k=4,
                model = "GTR")
fitGTRGI_opt <- optim.pml(fitGTRGI, model="GTR", optInv=TRUE, 
                          optGamma=TRUE, k=4)
treeMLGTRGI_hiv_pol  <- fitGTRGI_opt$tree

treeMLGTRGI_hiv_pol_r <- root.phylo(treeMLGTRGI_hiv_pol,  
                               outgroup = c("AY156771.1|LA02.RT "), 
                               resolve.root = T)

plot.phylo(treeMLGTRGI_hiv_pol_r, cex=0.6,
           tip.color = myColsStrings(treeMLGTRGI_hiv_pol_r$tip.label, c("V","P"),
                                                            c("red","blue")), 
           main="ML (GTR+G(4)+I) tree for HIV pol with colored tips")
legend("bottomleft", legend=c("Patient", "Victim", "Environment"),
       col=c("blue", "red", "black"), pch=16, cex=0.8)
# 7.) Explore the difference between the phylogenies resulting from
# the pol and env data sets for NJ. Assume 
# that you know that the outgroup is 'AY156764.1|LA26.EN'.
dm_hiv_env_JC  <- dist.dna(dna_hiv_env, model = "JC69")
treeNJ_hiv_env<- NJ(dm_hiv_env_JC)
labels(align_hiv_env)
# root the tree
treeNJ_hiv_env_r <- root.phylo(treeNJ_hiv_env,  
                               outgroup = c("AY156764.1|LA26.EN"), 
                               resolve.root = T)
# plot the tree
plot.phylo(treeNJ_hiv_env_r, cex=0.6,
           tip.color = myColsStrings(treeNJ_hiv_env_r$tip.label, c("V","P"),
                                                       c("red","blue")), 
           main="NJ tree for HIV env with colored tips")
legend("bottomleft", legend=c("Patient", "Victim", "Environment"),
       col=c("blue", "red", "black"), pch=16, cex=0.8)

#_______________________________________________________________________________

