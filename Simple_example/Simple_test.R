matrix <- read.csv("martins_matrix.csv", header = F, row.names=1)
library(inapplicable)

# Transforming a morphological matrix into a phyDat object (automatically!)
get.contrast.matrix <- function(matrix) {
    
    # Extracting the states
    states <- sort(unique(as.vector(matrix)))
    
    # Check if there is a "?" token
    if(any(states == "?")) {
        # remove the "?" state
        states_num <- states[-which(states == "?")]
        # Create a simple square matrix with 0s...
        contrast_matrix <- matrix(data = rep(0, length(states_num)*length(states_num)), ncol = length(states_num), dimnames = list(as.character(states_num), as.character(states_num)))
        # Set the diagonal to 0 
        diag(contrast_matrix) <- 1
        # Add the joker character as a row full of 1s
        joker_matrix <- matrix(data = rep(1, length(states_num)), ncol = length(states_num), dimnames = list("?", as.character(states_num)))
        contrast_matrix <- rbind(contrast_matrix, joker_matrix)
    } else {
        # Create a simple square matrix with 0s...
        contrast_matrix <- matrix(data = rep(0, length(states)*length(states)), ncol = length(states), dimnames = list(as.character(states), as.character(states)))
        # Set the diagonal to 0 
        diag(contrast_matrix) <- 1
    }

    return(contrast_matrix)
}




my.data <- apply(matrix, MARGIN = 2, function(x) return(as.character(x)))
rownames(my.data) <- rownames(matrix)

contrast.matrix <- get.contrast.matrix(my.data)

## To see the annotated contrast matrix, type
contrast.matrix

## Apply the contrast matrix, using the phyDat format...
my.phyDat <- phyDat(my.data, type='USER', contrast=contrast.matrix)

## ... and optimize the data for analysis 
my.prepdata <- prepare.data(my.phyDat)

## Specify the names of the outgroup taxa
my.outgroup <- c('taxon_1')

tree <- root(nj(dist.hamming(my.phyDat)), my.outgroup, resolve.root=TRUE); tree$edge.length <- NULL;
tree <- rtree(12, tip.label = rownames(my.data))

## Calculate the tree's parsimony score
parsimony.inapp(tree, my.prepdata)

# Search for a better tree
better.tree <- tree.search(tree, my.prepdata, outgroup=my.outgroup, method='SPR')
    #TG: outgroup needs to be specified as an option (it doesn't seem to be passed to tree.search anyway...)

## Once you have reached the most parsimonious tree, retain multiple trees to determine the consensus: 
best.trees <- tree.search(better.tree, my.prepdata,
      method='SPR', forest.size=100, maxiter=10000, maxhits=100)

## View the results
plot(better.tree)

## Once you have reached the most parsimonious tree, retain multiple trees to determine the consensus: 
best.trees <- tree.search(better.tree, my.prepdata, my.outgroup,
      method='SPR', forest.size=100, maxiter=10000, maxhits=100)

## Calculate and display the consensus tree
plot(my.consensus <- consensus(best.trees))