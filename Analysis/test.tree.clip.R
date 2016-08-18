
# Number of expected trees from SPR
number.of.SPR.trees <- function(n) {
    return(2*(n-3)*(2*(n-3)-1))
}

number.of.unrooted.trees <- function(n){
    return( factorial(2*n - 5) / (2^(n-3) * factorial(n-3)))
}

max.number.of.TBR <- function(n) {
    return((2*n-3)*(n-3)^2)
}


# Create a random tree with n taxa (unrooted + with all edges lengths = 1)
make.test.tree <- function(n) {
    return(rtree(n = n, rooted = FALSE, br = rep(1, (2*n - 3))))
}

# Creating a "one tip tree".
create.tip.tree <- function(tip.label) {
    new.tip <- list(edge = matrix(c(2,1),1,2), tip.label = tip.label, edge.length = 1, Nnode = 1)
    class(new.tip) <- "phylo"
    return(new.tip)
}

# Rebranching a clip on edges of tree
rebranch.on.edges <- function(tree, clip, edges) {
    #Creating an empty tree buffer
    buffer <- list()

    #Loop through the available edges
    for(edge in edges){
        buffer[[edge]] <- bind.tree(tree, clip, where = tree$edge[edge,2], position = 0.5)
        buffer[[edge]]$edge.length <- NULL
    }    

    return(buffer)
}


# Get a random tree
set.seed(12345) # WARNING DOES NOT WORK FOR EVERY RANDOM SEED YET! (problem with the way to clip on the edge in R trees)
tree <- make.test.tree(5)
buffer <- list()

#First clip (1 tip)
current_clip <- create.tip.tree(tree$tip.label[1])
current_tree <- unroot(drop.tip(tree, tip = tree$tip.label[1], rooted = FALSE))
buffer_tmp <- rebranch.on.edges(tree = current_tree, clip = current_clip, edges = 1:nrow(current_tree$edge))
#remove the first tree
buffer_tmp[[1]] <- NULL
buffer <- c(buffer, buffer_tmp)

#Initiate the clade clip
node_clip <- tree$edge[which(tree$edge[,2] == 1),1]

#Do the clade clip
while(length(node_clip) != 0) {
    #Clade clips
    current_clip <- extract.clade(tree, node = node_clip)
    if(Ntip(current_clip) < (Ntip(tree)-2)) {
        current_tree <- unroot(drop.tip(tree, tip  = current_clip$tip.label))
        buffer_tmp <- rebranch.on.edges(tree = current_tree, clip = current_clip, edges = 1:nrow(current_tree$edge))
        #remove the first tree
        buffer_tmp[[1]] <- NULL
        buffer <- c(buffer, buffer_tmp)
    }
    node_clip <- tree$edge[which(tree$edge[,2] == node_clip),1]
}

#Now proceed to the other clips but not allowing to branch on the neighbours

for(tip in 2:Ntip(tree)){
    current_clip <- create.tip.tree(tree$tip.label[tip])
    current_tree <- unroot(drop.tip(tree, tip = tree$tip.label[tip], rooted = FALSE))
    buffer_tmp <- rebranch.on.edges(tree = current_tree, clip = current_clip, edges = 4:nrow(current_tree$edge)) # No branching on the first three edges
    #remove the first tree
    buffer_tmp[[1]] <- NULL
    buffer <- c(buffer, buffer_tmp)

    #Initiate the clade clip
    node_clip <- tree$edge[which(tree$edge[,2] == tip),1]

    #Do the clade clip
    while(length(node_clip) != 0) {
        #Clade clips
        current_clip <- extract.clade(tree, node = node_clip)
        if(Ntip(current_clip) < (Ntip(tree)-2)) {
            current_tree <- unroot(drop.tip(tree, tip  = current_clip$tip.label))
            buffer_tmp <- rebranch.on.edges(tree = current_tree, clip = current_clip, edges = 3:nrow(current_tree$edge)) # No branching on the first three edges
            #remove the first tree
            buffer_tmp[[1]] <- NULL
            buffer <- c(buffer, buffer_tmp)
        }
        node_clip <- tree$edge[which(tree$edge[,2] == node_clip),1]
    }
}

#Transform the buffer in multiphylo
trees <- buffer[!sapply(buffer, is.null)]
class(trees) <- "multiPhylo"
