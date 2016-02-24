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
