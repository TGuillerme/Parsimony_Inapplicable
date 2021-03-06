#' @title Generates a morphological matrix.
#'
#' @description Generates a morphological matrix using \code{\link[ape]{rTraitDisc}} function with some optional inapplicable characters and allows to save it in various formats out of \code{R} environment.
#'
#' @param tree A phylogenetic tree to use for generating the characters.
#' @param characters The number of morphological characters.
#' @param model Either an implemented (\code{"ER"} or \code{"HKY"}; see details) or user defined model for generating one discrete morphological characters with at least the following arguments: \code{tree}, \code{states}, \code{rates}, \code{substitution}.
#' @param states A string of probabilities for the number of states for each characters (\code{default = 1}; i.e. 100\% binary state characters; see details).
#' @param rates A function an it's parameters for the rates distribution (see details).
#' @param substitution A function an it's parameters for the substitutions distribution (see details; \code{default = c(runif, 2, 2)}).
##' @param ... Any optional arguments to be passed to the model argument.
#' @param invariant Whether to allow any invariant sites.
##' @param inapplicable Optional, a vector of characters inapplicability source (either \code{"character"} or \code{"clade"}; see details). The length of this vector must be at maximum half the total number of characters.
##' @param output Optional, an output file name for writing the matrix out of the \code{R} environement in \code{nexus} format.
#' @param verbose Whether to be verbose or not.
#'
#' @details
#' \itemize{
#' 
#' \item The \code{model} arguments must be either a user's defined function for generating the discrete morphological characters (that intakes the states, rates and substitution arguments) or one of the two following:
#'      \itemize{
#'          \item \code{"ER"} uses the \code{ape::rTraitDisc} function with the \code{"ER"} model argument (= Mk model).
#'          \item \code{"HKY"} uses the \code{phyclust::gen.seq.HKY} function with \code{kappa} sampled from the \code{substitution} argument, \code{pi = runif(4)} (divided by \code{sum(runif(4))}), \code{rate.scale} sampled from the \code{rates} distribution and \code{L} being the number of \code{characters} and transforms the purines (A, G) into 0 and the pyrimidines (C, T) into 1.
#'      }
#'
#' \item The \code{states} argument attributes a number of states to each character by using the given probability vector for each number of states starting from 2.
#' For example \code{states = c(0.7, 0.2, 0.1)} will generate 70% of characters with 2 states, 20% of characters with 3 states and 10% of characters with 4 states. 
#' 
#' \item The \code{rates} and \code{substitution} arguments attributes a distribution function and it's optional parameters to a model. For example \code{rates = c(runif, 1, 10)} attributes a uniform distribution between 1 and 10 for the rates distribution.
#' 
# ' \item The \code{inapplicables} argument intakes a vector of character inapplicability source rendering a number of characters inapplicable using the following sources:
# '      \itemize{
# '          \item \code{"character"} draws inapplicable characters directly from the character matrix, ignoring the phylogeny (i.e. for a random character X, an other random character Y will have inappicable characters fro each character states 0 for character X).
# '          \item \code{"clade"} draws inapplicable characters from the phylogeny: it will randomly apply inapplicable characters states for some characters by randomly selecting clades from the provided tree.
# '      }
# ' For example \code{inapplicables = c(rep("character", 2), rep("clade", 2))} will generate 4 characters with inapplicable data, two using previous characters and two other using random clades.
#' 
#' }
#' 
#' @examples
#' ## A random tree with 10 tips
#' tree <- rcoal(10)
#' ## setting up the parameters
#' my_rates = c(rgamma, 1, 1) # A gamma rate distribution with of shape alpha = 0.5
#' my_substitutions = c(runif, 2, 2) # A fixed substitution rate of 2 (T/T ratio in HKY)
#'
#' ## Mk matrix (10*50)
#' matrixMk <- make.matrix(tree, characters = 50, model = "ER", rates = my_rates)
#' ## HKY binary (10*50)
#' matrixHKY <- make.matrix(tree, characters = 50, model="HKY", rates = my_rates, substitution = my_substitutions)
#' @author Thomas Guillerme
#' @export


make.matrix <- function(tree, characters, states = 1, model = "ER", rates, substitution = c(runif, 2, 2), invariant = FALSE, verbose = FALSE) { #...

    #SANITIZNG
    #tree
    check.class(tree, "phylo")

    #characters
    check.class(characters, "numeric")
    check.length(characters, 1, " must be a single numeric value.")

    #states
    check.class(states, "numeric")
    if(sum(states) != 1) {
        stop("States argument must sum up to 1.")
    }
    #Set to 1 if model is HKY
    if(length(states) > 1 && model == "HKY") {
        states <- 1
        warning("The HKY model only allows the use of two states characters (binary).")
    }

    #model
    if(class(model) != "function") {
        #model is not a sure function
        implemented_models <- c("ER", "HKY")
        if(all(is.na(match(model, implemented_models)))) stop("The model must be either a user's function or one of the following: ", paste(implemented_models, collapse=", "), sep="")
        #Setting up the model
        if(class(model) != "function" && model == "ER") {
            model <- rTraitDisc.mk
            #Warning on the substitutions:
            substitution <- c(runif, 1, 1)
            #message("Substitution parameter is ignored for the ER model.")
        }
        if(class(model) != "function" && model == "HKY") model <- gen.seq.HKY.binary

    } else {
        stop("User functions not implemented yet for model argument.")
        #Add checker for arguments to be passed to users function
    }

    #rate
    if(length(rates) == 1) {
        # Is only a function (default function arguments)
        check.class(rates, "function")
    } else {
        # Is a function with options
        check.class(rates[[1]], "function")
        # Check if arguments work
        test <- NULL ; try(test <- sample.distribution(1, rates), silent = TRUE)
        if(length(test) != 1) {
            stop("Error in rates argument format. Should be c(function, arg1, arg2, ...).")
        }
    }

    #substitution
    if(length(substitution) == 1) {
        # Is only a function (default function arguments)
        check.class(substitution, "function")
    } else {
        # Is a function with options
        check.class(substitution[[1]], "function")
        # Check if arguments work
        test <- NULL ; try(test <- sample.distribution(1, substitution), silent = TRUE)
        if(length(test) != 1) {
            stop("Error in substitution argument format. Should be c(function, arg1, arg2, ...).")
        }
    }

    #invariant
    check.class(invariant, "logical")

    #verbose
    check.class(verbose, "logical")


    #GENERATING THE CHARACTERS
    #Isolating the arguments
    #arguments <- as.list(substitute(list(tree = tree, states = states, rates = rates, substitution = substitution, ...)))[-1L]

    #Creating the matrix
    if(verbose == TRUE) cat(paste("Generating a matrix of ", characters, " characters for ", Ntip(tree), " taxa:...", sep=""))
    #matrix <- replicate(characters, do.call(model, arguments))
        matrix <- replicate(characters, model(tree = tree, states = states, rates = rates, substitution = substitution))
    if(verbose == TRUE) cat("Done.\n")


    if(invariant == FALSE) {
        if(any(apply(matrix, 2, is.invariant))) {
            if(verbose == TRUE) cat("Re-simulating ", length(which(apply(matrix, 2, is.invariant)) == TRUE), " invariant characters:", sep="") 
            #Repeat the invariant characters sampling
            while(any(apply(matrix, 2, is.invariant))) {
                #matrix[, which(apply(matrix, 2, is.invariant) == TRUE)] <- replicate(length(which(apply(matrix, 2, is.invariant) == TRUE)), do.call(model, arguments))
                matrix[, which(apply(matrix, 2, is.invariant) == TRUE)] <- replicate(length(which(apply(matrix, 2, is.invariant) == TRUE)), model(tree = tree, states = states, rates = rates, substitution = substitution))
                if(verbose == TRUE) cat(".")
            }
            if(verbose == TRUE) cat("Done.\n")
        }
    }

    #Adding the row names
    rownames(matrix) <- tree$tip.label
    return(matrix)
    

    
}