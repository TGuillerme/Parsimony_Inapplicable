#' @title Check morphological matrix consistency.
#'
#' @description Performs a quick and dirty checking of the phylogenetic signal in a morphological matrix using parsimony.
#'
#' @param matrix A discrete morphological matrix.
#' @param parsimony Either the parsimony algorithm to be passed to \code{\link[phangorn]{optim.parsimony}} or a parsimony function that can take a \code{\link[phangorn]{phyDat}} object as an input.
#' @param first.tree A list of functions to generate the first parsimonious tree (default = \code{c(\link[phangorn]{dist.hamming}, \link[phangorn]{NJ})}; see details).
#' @param orig.tree Optional, the input tree to measure the distance between the parsimony and the original tree.
#' @param distance Optional, if orig.tree is provide, the function to use for measuring distance between the trees (default = \code{link[phangorn]{RF.dist}}).
#' @param ... Any additional arguments.
#' @param contrast.matrix An optional contrast matrix. By default, the function recognises any character states token as different apart from \code{?} that is treated as a joker character.
#'
#' @return
#' Returns the parsimony score (using \code{\link[phangorn]{parsimony}}), the consistency and retention indices (using \code{\link[phangorn]{CI}} and \\code{\link[phangorn]{RI}}) from quick most parsimonious tree obtained from the matrix.
#' Can also return the topological distance from the original tree if provided.
#' 
#' #' @details
#' \itemize{
#' \item The \code{first.tree} argument must be a list of functions to be use in a cascade to transform the matrix (as a \code{\link[phangorn]{phyDat}} object) into a tree using the functions iteratively.
#' For example the default \code{c(\link[phangorn]{dist.hamming}, \link[phangorn]{NJ})} will apply the following to the matrix: \code{\link[phangorn]{NJ}(\link[phangorn]{dist.hamming}(matrix))}
#' }
#' 
#' @examples
#' ## Generating a random tree
#' random.tree <- rcoal(10)
#' ## Generating a random matrix
#' random.matrix <- make.matrix(random.tree, characters = 50, model = "ER", rates = c(rgamma, 1, 1))
#'
#' ## Checking the matrix scores
#' check.matrix(random.matrix, orig.tree = random.tree)
#' 
#' @author Thomas Guillerme
#' @export

check.matrix <- function(matrix, orig.tree, parsimony = "fitch", first.tree = c(phangorn::dist.hamming, phangorn::NJ), distance = phangorn::RF.dist, ..., contrast.matrix) {
    #SANITIZNG

    #matrix
    check.class(matrix, "matrix")

    #parsimony
    if(class(parsimony) != "function") {
        #model is not a sure function
        implemented_parsimony <- c("fitch", "sankoff")
        if(all(is.na(match(parsimony, implemented_parsimony)))) stop("The parsimony argument must be either a user's function or one of the following: ", paste(implemented_parsimony, collapse=", "), sep="")
        #setting the parsimony algorithm
        use.optim.parsimony <- TRUE
        parsimony.algorithm <- phangorn::optim.parsimony
        method <- parsimony
    } else {
        stop("User functions not implemented yet for model argument.")
        use.optim.parsimony <- FALSE
        parsimony.algorithm <- parsimony
    }

    #first.tree
    if(class(first.tree) != "function") {
        if(any(unlist(lapply(first.tree, class)) != "function")) {
            stop("first.tree argument must be a list of functions to calculate the first tree.")
        }
    }

    #orig.tree
    if(!missing(orig.tree)) {
        check.class(orig.tree, "phylo")
        #must be same size as the matrix
        if(any(sort(row.names(matrix)) != sort(orig.tree$tip.label))) {
            stop("Provided orig.tree has not the same number of taxa as the matrix.")
        }
        #distance
        check.class(distance, "function")
    }

    #contrast.matrix
    if(!missing(contrast.matrix)) {
        check.class(contrast.matrix, "matrix")
    }


    #CHECKING THE MATRIX

    #Creating the contrast matrix
    if(missing(contrast.matrix)) {
        constrast.matrix <- get.contrast.matrix(matrix)
    }

    #Creating the phyDat object
    matrix_phyDat <- phangorn::phyDat(data = matrix, type = "USER", contrast = constrast.matrix)

    #Calcualte the first tree
    if(class(first.tree) == "function") {
        #Calculate the first tree from the phyDat
        first_tree <- first.tree(matrix_phyDat)
    } else {
        #Calculate the operation on the phyDat
        first_tree <- first.tree[[1]](matrix_phyDat)
        for(operations in 2:length(first.tree)) {
            #Transform the first_tree using the other operations
            first_tree <- first.tree[[operations]](first_tree)
        }
    }

    #Get the quick and dirty most parsimonious tree
    if(use.optim.parsimony == TRUE) {
        cat("Most parsimonious tree search:\n")
        MP_tree <- parsimony.algorithm(tree = first_tree, data = matrix_phyDat, method = method, ...)
        #MP_tree <- parsimony.algorithm(tree = first_tree, data = matrix_phyDat, method = method) ; warning("DEBUG")
    } else {
        MP_tree <- parsimony.algorithm(tree = first_tree, data = matrix_phyDat)
    }

    #Get the parsimony score
    parsimony_score <- phangorn::parsimony(MP_tree, matrix_phyDat)

    #Get the CI
    consistency_index <- phangorn::CI(MP_tree, matrix_phyDat)

    #Get the retention Index
    retention_index <- phangorn::RI(MP_tree, matrix_phyDat)    

    if(!missing(orig.tree)) {
        #Get the distance between the trees
        tree_distance <- phangorn::RF.dist(MP_tree, unroot(orig.tree))

        #Get the data out vectors (with distance)
        data_out <- c(parsimony_score, consistency_index, retention_index, tree_distance)
        out_names <- c("Maximum parsimony", "Consistency index", "Retention index", "Robinson-Foulds distance")

    } else {

        #Get the data out vectors (without distance)
        data_out <- c(parsimony_score, consistency_index, retention_index)
        out_names <- c("Maximum parsimony", "Consistency index", "Retention index")
    }

    return(matrix(data = data_out, ncol = 1, dimnames = list(out_names, c(""))))
}
