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
# 1.) Load the HIV-1 pol and env datasets.
# Alignment HIV env
align_hiv_env <- NULL # replace NULL with your code

# Alignment HIV pol
align_hiv_pol <- NULL # replace NULL with your code

# 3.) For pol choose a first tree reconstruction method, e.g. NJ 
#     i.e. distance-based) and reconstruct trees.
treeNJ_hiv_pol <- NULL # replace NULL with your code
# 4.) Have a look at the labels for pol. 
#     All labels with 'V' are from the victim and with the ones with 'P' are 
#     from the patient ('LA' is data from environment)
labels(align_hiv_pol)
#     Assume that you know that the outgroup is 'AY156771.1|LA02.RT ' and root
#     the tree accordingly.
treeNJ_hiv_pol_r <- NULL # replace NULL with your code
# 5.) Plot the rooted tree and look for the evidence. Color the patient tips in
#     blue and the victim tips in red, you can use your myColsStrings function
#     from 1IntroductionToR.R.

# 6.) Check one of the other methods (MP or ML) similarly. Do they yield the
#     same result?

# 7.) Explore the difference between the phylogenies resulting from
# the pol and env data sets for NJ. Assume 
# that you know that the outgroup is 'AY156764.1|LA26.EN'.

#_______________________________________________________________________________

