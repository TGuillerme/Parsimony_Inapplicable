library(phangorn)

#Get the edge table of a tree and sort it
get.edges.table <- function(tree, sort = TRUE) {
    edge_table <- tree$edge
    if(sort != FALSE){
        edge_table <- edge_table[order(edge_table[,2]), ]
    }
    return(edge_table)
}

#Transform the edge table into a hash value
get.hashed.edge <- function(edge_table) {
    hash_vector <- as.vector(apply(edge_table, 1, c))
    return(as.numeric(paste(hash_vector, collapse = "")))
}

#Get hashes from list
get.hash.from.tree <- function(tree, ...) {
    return(get.hashed.edge(get.edges.table(tree, ...)))
}

# Get all trees for n taxa
length_hash <- NULL
length_hash_unique <- NULL
length_tree <- NULL
results <- NULL
for (n in 4:10) {
    all_tree_list <- allTrees(n)

    # Get all the hashes
    hash_vector <- unlist(lapply(all_tree_list, get.hash.from.tree))

    # Store the results
    length_hash[n-3] <- length(hash_vector)
    length_hash_unique[n-3] <- length(unique(hash_vector))
    length_tree[n-3] <- length(all_tree_list)

    # Check if they're all unique
    results[n-3] <- length(unique(hash_vector)) == length(all_tree_list)
}