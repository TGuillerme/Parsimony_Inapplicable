#' @title Apply inapplicable characters to a matrix.
#'
#' @description Apply inapplicable characters to discrete morphological matrix.
#'
#' @param matrix A discrete morphological matrix.
#' @param inapplicables Optional, a vector of characters inapplicability source (either \code{"character"} or \code{"clade"}; see details). The length of this vector must be at maximum half the total number of characters.
#' @param invariant Whether to allow any invariant sites among the characters with inapplicable data.
#' @param ... Any additional arguments.
#' 
#' @details
#' \itemize{
#' \item The \code{inapplicables} argument intakes a vector of character inapplicability source rendering a number of characters inapplicable using the following sources:
#'      \itemize{
#'          \item \code{"character"} draws inapplicable characters directly from the character matrix, ignoring the phylogeny (i.e. for a random character X, an other random character Y will have inappicable characters fro each character states 0 for character X).
#'          \item \code{"clade"} draws inapplicable characters from the phylogeny: it will randomly apply inapplicable characters states for some characters by randomly selecting clades from the provided tree.
#'      }
#' For example \code{inapplicables = c(rep("character", 2), rep("clade", 2))} will generate 4 characters with inapplicable data, two using previous characters and two other using random clades.
#' 
#' }
#' 
#' @examples
#'
#' @author Thomas Guillerme
#' @export

apply.inapplicables <- function(...) {

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