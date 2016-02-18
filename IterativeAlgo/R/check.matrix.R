#' @title Check morphological matrix consistency.
#'
#' @description Performs a quick and dirty checking of the phylogenetic signal in a morphological matrix using parismony.
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

check.matrix <- function(matrix, parismony = "fitch", first.tree = c(dist.hamming, NJ), orig.tree, distance = RF.dist, ...) {
    #SANITIZNG

    #matrix
    check.class(matrix, "matrix")

    #parsimony
    if(class(parsimony) != "function") {
        #model is not a sure function
        implemented_parsimony <- c("fitch", "sankoff")
        if(all(is.na(match(parsimony, implemented_parsimony)))) stop("The parsimony argument must be either a user's function or one of the following: ", paste(implemented_parsimony, collapse=", "), sep="")
        #setting the parsimony algorithm
        use.optim.parisomy <- TRUE
        parsimony.algorithm <- phangorn::optim.parsimony
        method <- parsimony
    } else {
        stop("User functions not implemented yet for model argument.")
        use.optim.parisomy <- FALSE
        parsimony.algorithm <- parsimony
    }

    #first.tree
    if(class(first.tree) != "function") {
        if(any(unlist(lapply(first.tree, class))) != "function") {
            stop("first.tree argument must be a list of functions to calculate the first tree.")
        }
    }

    #orig.tree
    check.class(orig.tree, "phylo")
    #must be same size as the matrix
    if(nrow(matrix) != Ntip(orig.tree)) {
        stop("Provided orig.tree has not the same number of taxa as the matrix.")
    }

    #distance
    check.class(distance, "function")


    #CHECKING THE MATRIX

    #Creating the contrast matrix
    constrast.matrix <- get.contrast.matrix(matrix)

    #Creating the phyDat object
    matrix_phyDat <- phangorn::phyDat(data, type = "USER", contrast = constrast.matrix)

}




# contrast.matrix <- matrix(data=c(
#   1,0,0, # 0 = {0}
#   0,1,0, # 1 = {1}
#   0,0,1, # - = {-}
#   1,1,0, # A = {01}
#   1,1,0, # + = {01}
#   1,1,1  # ? = {01-}
# ), ncol=3, byrow=TRUE)

# # Dimnames must be the list of tokens (in the matrix) and the list of character states
# dimnames(contrast.matrix) <- list(
#   c(0, 1, '-', 'A', '+', '?'), # A list of the tokens corresponding to each rows in the contrast matrix
#   c(0, 1, '-') # A list of the character-states corresponding to the columns in the contrast matrix
# )

# # Transforming my.data from list to matrix
# my.data <- matrix(data = unlist(my.data) , nrow= length(my.data), byrow=TRUE, dimnames=list(names(my.data)))

# # Applying the contrast matrix to the matrix (through phangorn::phyDat)
# my.phyDat <- phyDat(my.data, type='USER', contrast=contrast.matrix)

# # phangorn::parsimony

# dm = phangorn::dist.hamming(Laurasiatherian)
# tree = phangorn::NJ(dm)

# phangorn::optim.parsimony(tree, data, method="fitch", cost=NULL, trace=1, rearrangements="SPR", ...)
# phangorn::CI(tree, data, cost = NULL, sitewise=FALSE)
# phangorn::RI(tree, data, cost = NULL, sitewise=FALSE)