#' @title Generates a birth death tree.
#'
#' @description Generates a birth death tree with a set number of taxa and a random birth and death rate.
#'
#' @param n an integer giving the number of tips in the tree.
#' 
#' @author Thomas Guillerme
#' @export

rtree.bd <- function(n) {
    #Random parameters selector
    rand.birth.death <- function() {
        lambda <- runif(1)
        mu <- runif(1, 0, lambda)
        return(cbind(lambda, mu))
    }

    tree <- diversitree::tree.bd(rand.birth.death(), max.taxa = n)

    # Make sure tree is not null
    while(is.null(tree)) {
        tree <- diversitree::tree.bd(rand.birth.death(), max.taxa = n)
    }

    return(tree)
}


#' @title Generates a contrast matrix.
#'
#' @description Creates a contrast matrix using the observed character states in an input matrix.
#'
#' @param matrix a discrete character matrix.
#' 
#' @author Thomas Guillerme
#' @export

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

