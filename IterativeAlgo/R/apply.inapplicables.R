#' @title Apply inapplicable characters to a matrix.
#'
#' @description Apply inapplicable characters to discrete morphological matrix.
#'
#' @param matrix A discrete morphological matrix.
#' @param inapplicables Optional, a vector of characters inapplicability source (either \code{"character"} or \code{"clade"}; see details). The length of this vector must be at maximum half the total number of characters.
#' @param tree If any inapplicable source is \code{"clade"}, a tree from where to select the clades.
#' @param invariant Whether to allow invariant sites among the characters with inapplicable data. If \code{invariant = FALSE} the algorithm will try to remove such characters (if possible).
#' @param verbose Whether to be verbose or not.
#' @param ... Any additional arguments.
#' 
#' @details
#' \itemize{
#' \item The \code{inapplicables} argument intakes a vector of character inapplicability source rendering a number of characters inapplicable using the following sources:
#'      \itemize{
#'          \item \code{"character"} draws inapplicable characters directly from the character matrix, ignoring the phylogeny (i.e. for a random character X, an other random character Y will have inappicable characters fro each character states 0 for character X).
#'          \item \code{"clade"} draws inapplicable characters from the phylogeny: it will randomly apply inapplicable characters states for some characters by randomly selecting clades from the provided tree. The algorithm randomly assigns an inapplicable token for this character for all taxa in this clade or all taxa outside this clade.
#'      }
#' For example \code{inapplicables = c(rep("character", 2), rep("clade", 2))} will generate 4 characters with inapplicable data, two using previous characters and two other using random clades.
#' 
#' }
#' 
#' @examples
#' set.seed(4)
#' ## A random tree with 15 tips
#' tree <- rcoal(15)
#' ## setting up the parameters
#' my_rates = c(rgamma, 1, 1) # A gamma rate distribution with of shape alpha = 0.5
#' my_substitutions = c(runif, 2, 2) # A fixed substitution rate of 2 (T/T ratio in HKY)
#'
#' ## A Mk matrix (10*50)
#' matrixMk <- make.matrix(tree, characters = 100, model = "ER", states = c(0.85, 0.15), rates = my_rates)
#' 
#' ## Setting the number and source of inapplicable characters
#' my_inapplicables <- c(rep("character", 5), rep("clade", 5))
#' 
#' ## Apply some inapplicable characters to the matrix
#' matrix <- apply.inapplicables(matrixMk, my_inapplicables, tree)
#' @author Thomas Guillerme
#' @export

apply.inapplicables <- function(matrix, inapplicables, tree, invariant = FALSE, verbose = FALSE, ...) {

    #SANITIZING
    #matrix
    check.class(matrix, "matrix")

    #inapplicables
    check.class(inapplicables, "character")
    inap.source_options <- c("character", "clade")
    if(all(is.na(match(unique(inapplicables), inap.source_options)))) {
        stop("inapplicables argument must be a vector containing at least one of the following: ", paste(inap.source_options, collapse=", "), sep="")
    }        
    if(length(inapplicables) > ncol(matrix)/2) {
        stop("Only half the number of characters can be inapplicables")
    }

    #tree
    if(any(inapplicables == "clade") && missing(tree)) {
        stop("Tree argument is missing for applying inapplicable characters on random clades.")
    } else {
        #tree must be same size as the matrix
        if(any(sort(row.names(matrix)) != sort(tree$tip.label))) {
            stop("Provided tree has not the same number of taxa as the matrix.")
        }
    }

    #invariant
    check.class(invariant, "logical")

    #verbose
    check.class(verbose, "logical")

    #APPLY INAPPLICABLE CHARACTERS
    #Setting the output matrix
    matrix_out <- matrix

    #From the characters first (if any)
    if(any(inapplicables == "character")) {
        #Get the number of inapplicable characters
        inapplicables_characters <- length(which(inapplicables == "character"))
        #Create the characters to make inapplicable
        target_characters <- seq(from = 2, to = inapplicables_characters*2, by = 2)
        #Create the characters to match inapplicability from
        pattern_characters <- seq(from = 1, to = inapplicables_characters*2, by = 2)

        #Get the target characters with inapplicables
        matrix_inapplicable <- mapply(mapply.inap.character, as.list(target_characters), as.list(pattern_characters), MoreArgs=list(matrix, invariant))

        #Invariant warning
        if(invariant == FALSE) {
            invariants <- length(which(apply(matrix_inapplicable, 2, function(X) length(unique(X))) <= 2))
            if(invariants != 0) {
                message(paste(invariants, "characters are now invariant due inapplicable data."))
            }
        }

        #Include these characters in the matrix
        matrix_out[,target_characters] <- matrix_inapplicable
    } else {
        inapplicables_characters <- 0
    }

    #From the clades (if any)
    if(any(inapplicables == "clade")) {
        #Get the number of inapplicable characters
        inapplicables_clades <- length(which(inapplicables == "clade"))
        #Create the characters to make inapplicable (after the inapplicable characters if any)
        target_characters <- seq(from = inapplicables_characters*2+1, to = inapplicables_characters*2+inapplicables_clades)

        #Get the target characters with inapplicables
        matrix_inapplicable <- matrix(unlist(lapply(as.list(target_characters), lapply.inap.clade, matrix, tree, invariant)), ncol = inapplicables_clades, byrow = FALSE, dimnames = list(rownames(matrix)))

        #Include these characters in the matrix
        matrix_out[,target_characters] <- matrix_inapplicable
    }

    return(matrix_out)
}