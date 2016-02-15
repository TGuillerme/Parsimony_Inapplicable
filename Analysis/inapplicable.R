# Playing around with the inapplicable package

# Get the package
library(devtools)
install("../inapplicable")
library(inapplicable)

# Walkthrough the package example


# Loading inbuilt data
data(SigSut); my.data <- SigSut.data
# List of 34 taxa with 5 tokens
length(my.data) ; unique(unlist(lapply(my.data, unique)))

## A contrast matrix translates the tokens used in your dataset to the character states to 
##    which they correspond: for example decoding 'A' to {01}.
##    For more details, see the 'phangorn-specials' vignette in the phangorn package, accesible 
##    by typing '?phangorn' in the R prompt and navigating to index > package vignettes.

# The contrast matrix is a a token translator code (e.g. 'A' = {01})
# The matrix is a presence absence matrix with the rows being the states in the matrix (0, 1, A, ?, etc.) and the columns being the real characters to be interpreted (i.e. ? = 0, 1, -)
contrast.matrix <- matrix(data=c(
  1,0,0, # 0 = {0}
  0,1,0, # 1 = {1}
  0,0,1, # - = {-}
  1,1,0, # A = {01}
  1,1,0, # + = {01}
  1,1,1  # ? = {01-}
), ncol=3, byrow=TRUE)

# Dimnames must be the list of tokens (in the matrix) and the list of character states
dimnames(contrast.matrix) <- list(
  c(0, 1, '-', 'A', '+', '?'), # A list of the tokens corresponding to each rows in the contrast matrix
  c(0, 1, '-') # A list of the character-states corresponding to the columns in the contrast matrix
)

# Transforming my.data from list to matrix
my.data <- matrix(data = unlist(my.data) , nrow= length(my.data), byrow=TRUE, dimnames=list(names(my.data)))

# Applying the contrast matrix to the matrix (through phangorn::phyDat)
my.phyDat <- phyDat(my.data, type='USER', contrast=contrast.matrix)
	#TG: There seems to be a bug here... Matrix doesn't seem properly transformed/converted in phyDat.


# Setting up the phyDat object to be passed to the following inapplicable functions
my.prepdata <- prepare.data(my.phyDat)

# Setting an outgroup
my.outgroup <- c('taxon1', 'taxon2')

# Generate a random zero branch length starting tree
tree <- rtree(length(my.phyDat), rooted=TRUE, tip.label=names(my.phyDat), br=NULL)
	#TG: "or use neighbour joining to generate a starting tree?"
	#TG: "or load a bifurcating tree?"
	#TG: Is it better to start with a random tree? I'll argue yes but then run multiple chains to start from different starting point?


# Calculate the tree's parsimony score
parsimony.inapp(tree, my.prepdata) # score depends on the starting tree right?

# Search for a better tree
better.tree <- tree.search(tree, my.prepdata, outgroup=my.outgroup, method='SPR')
	#TG: outgroup needs to be specified as an option (it doesn't seem to be passed to tree.search anyway...)

# Try a sectorial search (Goloboff, 1999)
better.tree <- sectorial.search(better.tree, my.prepdata)
	#TG: outgroup doesn't work!
	#TG: error line 237, ("method=") no method argument needed.

# Try the parsimony ratchet (Nixon, 1999)
better.tree <- pratchet.inapp(better.tree, my.prepdata, my.outgroup, maxhits=50, k=20)

# The default parameters may not be enough to find the most parsimonious tree; type 
#    ?pratchet.inapp or ?sectorial.search to view all search parameters.


## View the results
plot(better.tree)

## Once you have reached the most parsimonious tree, retain multiple trees to determine the consensus: 
best.trees <- tree.search(better.tree, my.prepdata,
      method='SPR', forest.size=100, maxiter=10000, maxhits=100)

## Calculate and display the consensus tree
plot(my.consensus <- consensus(best.trees))

## End(Not run)
