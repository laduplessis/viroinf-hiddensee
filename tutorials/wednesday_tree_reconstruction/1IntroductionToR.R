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
# 1) Introduction to R
#-------------------------------------------------------------------------------
################################################################################

# Create your own directory for the tutorial and copy all R scripts!
# Now, you can set your working directory there (edit the example path).
setwd("C:\\Users\\Sophie\\Desktop\\ViroInfTutorial")
# Useful shortcuts to remember:
# Execute next command:   Ctrl+Enter (Windows/Linux)  Cmd+Return (Mac)
# Create "<-" symbol:     Alt+- (Windows/Linux)       Option+- (Mac)
# Show help for marked function: F1
# Show code of marked function: F2
# use the command 'dev.off()' to reset the plot settings

################################################################################
# a) Variables and first functions
# ----- all variables are defined with the "<-" symbol
a <- 2      
b <- "hello"
c3_Po <- 5.7*sqrt(12)
a+c3_Po
# ----- we can create vectors with c(element1, element2, ...)
my_vec <- c("hello", "HI", "high", "virus") 
rep_vec <- rep(c("Group A", "Group B"), 10)
num_vec <- 21:36
c(my_vec, rep_vec)
num_vec[3]
my_vec[1:2]
my_vec[c(4,3,1,2)]
length(my_vec)
# ----- these are some functions to deal with strings
grep(pattern = "i", my_vec, ignore.case = TRUE)
grepl(pattern = "i", my_vec)
paste("hello","hi", sep = " ")
paste0("hello","hi")
paste(my_vec, collapse = "', '")

#_______________________________________________________________________________
# Task: Experiment with the variables. For example, try...
# 1.) creating a vector containing strings and numbers,
test <- c(2, 3, -4, "test")
# 2.) filtering all of the following labels of dna sequences to get the ones 
# from the philippines,
LABS <- c('D4Brazi82', 'D4ElSal83', 'D4ElSal94', 'D4Indon76', 'D4Indon77', 
            'D4Mexico84', 'D4NewCal81', 'D4Philip64', 'D4PhilIp56', 
            'D4Philip84', 'D4PRico86', 'D4SLanka78', 'D4Tahiti79', 
            'D4Tahiti85', 'D4Thai63', 'D4Thai78', 'D4Thai84')
# 3.) accessing the last two elements of the LABS vector,
# 4.) and creating a variable containing the path to the dengue data set.
#_______________________________________________________________________________

################################################################################
# b) Trees and lists (format "phylo")
# the ape package can give us a random tree (generated with Yule process)
mytree <- ape::rtree(n=5, rooted = TRUE)
mytree2 <-  ape::rtree(n=10, rooted = TRUE)
ape::plot.phylo(mytree, type = "cladogram")
# access the parts of the tree:
mytree$edge
mytree$tip.label
mytree$Nnode
mytree$edge.length
# plot the tree
mytree$tip.label <- c("l1","l2","l3","l4","l5")
mytree$node.label <- c("i6","i7","i8","i9")
ape::plot.phylo(mytree, show.node.label = TRUE)
# lists of trees
moretrees <- list(firsttree = ape::rtree(n=3, rooted = T), 
                  secondtree = ape::rtree(n=4, rooted = T),
                  thirdtree = ape::rtree(n=16, rooted = T))
attr(moretrees, "class") <- "multiPhylo"
moretrees$firsttree
moretrees[[2]]
#_______________________________________________________________________________
# Task: Experiment with trees and lists. For example, try...
# 1.) experimenting with different numbers of leaves.
# 2.) setting rooted=FALSE and observing the difference 
# (plot with type="unrooted", "fan" or "radial).
# 3.) plotting the next tree (MYTREE) as a cladogram with half of the labels 
# colored in red, the other half in blue, larger text size and bold letters
# (-> consult the help page for more information on such options).
set.seed(42); MYTREE <- ape::rtree(n=10, rooted = T)
ape::plot.phylo(MYTREE)
#_______________________________________________________________________________


################################################################################
# c) Functions, if-else
# ------ Functions and if-else
if(3>4){ 
  print("hello")
}else{
  print("bye")
}
# ------ Functions (and if-else)
myMaxFunction <- function(x,y){
  if(x>y){
    # swap x and y
    temp <- y
    y <- x
    x <- temp
  } 
  return(y)
}
myMaxFunction(3,5); myMaxFunction(3,-5); myMaxFunction(1,2)

myGreatFunction <- function(how_often, strings, add_star){
  if(add_star){
    print(paste(strings, collapse = " * "))
  } else {
    print(rep(strings[1], how_often))
  }
}
myGreatFunction(3, c("hello", "hi", "virus"), TRUE)
myGreatFunction(strings = c("hello", "hi", "virus"), add_star = FALSE, 
                how_often = 3)
# ----- For-loops
for(i in c(-2,1,-4,5)){
    print(i)
}
for(i in 1:5){
  if(myMaxFunction(3,i)!=max(3,i)){
    print(paste("Oh no! Functions do not agree for i =",i))
  }
}
#_______________________________________________________________________________
# Task: Experiment with all of this. For example, try:
# 1.) Modify the For-Loop such that the error message is triggered.
# 2.) Create a function that receives a vector of strings, a vector of patterns 
# and a vector of color names as input.
# The function should create a vector of color names as long as the vector of
# strings. If the i-th string does not contain any pattern it gets color black,
# if it contains the first pattern it gets the first color, etc.
# Example: c("T","C","A","T") with pattern c("A","C") and colors 
# c("red", "blue") the function should yield c("black","blue","red","black")
myColsStrings <- function(strings, patterns, cols){
  #insert your code here
}
# 3.) Test it with the labels we already used in part a):
LABS <- c('D4Brazi82', 'D4ElSal83', 'D4ElSal94', 'D4Indon76', 'D4Indon77', 
          'D4Mexico84', 'D4NewCal81', 'D4Philip64', 'D4PhilIp56', 
          'D4Philip84', 'D4PRico86', 'D4SLanka78', 'D4Tahiti79', 
          'D4Tahiti85', 'D4Thai63', 'D4Thai78', 'D4Thai84')
myColsStrings(strings = LABS, patterns = c("Indon","Tahiti"), 
             cols = c("orange","green"))
set.seed(42); MYTREE <- ape::rtree(n=17, rooted = T)
MYTREE$tip.label <- LABS
ape::plot.phylo(MYTREE, type = "cladogram", 
                tip.color = myColsStrings(strings = MYTREE$tip.label, 
                                          patterns = c("Indon","Tahiti"), 
                                          cols = c("orange","green")), 
                font = 2, cex=1)
#_______________________________________________________________________________


################################################################################
# d) For-loops, the sapply function and a little bit of plotting with R
# ------ For-loops and the sapply function
# let's test our first function from part c) for 5 values:
myresults1 <- NULL
myvalues <- c(-4,-2,0,2,4)
for(myindex in 1:5){
  myresults1[myindex] <- myMaxFunction(3,myvalues[myindex])
}
myresults1
# let's do something similar with the sapply function
myresults2 <- sapply(c(-4,-2,0,2,4), FUN = function(x){myMaxFunction(3,x)})
myresults2
myresults1==myresults2
# ----- simple plot
plot(x=myvalues, y=myresults2, main="Example plot", xlab="Test values", 
     ylab="Maximum of 3 and test value")
lines(x=myvalues, y=myresults2)
#_______________________________________________________________________________
# Task: Experiment with all of this. For example, try:
# 1.) Use lapply (has the same usage as sapple but saves all results in a list)
# to create a list of 100 Yule trees with n=15 leaves using the ape:rtree 
# function. 
Yule_trees <- NULL # replace NULL with your code
# You can access the i-th tree by Yule_trees[[i]]
# 2.) How do the label sets of the trees look like? Use sapply to extract the 
# position (index) of label "t2" of each Yule tree.  (Hint: grep function)
Yule_trees[[1]]$tip.label
position_t2 <- NULL # replace NULL with your code
# 3.) Plot the results using a bar plot/histogram. Modify the histogram so that 
# we have 15 bars, add a title and a label for the x-axis, and then use
# abline to add a horizontal line indicating where the uniform distribution is.
hist(position_t2)
abline(h=100/15)
# 4.) Is the position of "t2" random? Do the same test for 1000 trees.
#_______________________________________________________________________________


