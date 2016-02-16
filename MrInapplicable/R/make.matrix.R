#' @title Generates a morphological matrix.
#'
#' @description Generates a morphological matrix using \code{\link[ape]{rTraitDisc}} function with some optional inaplicable characters and allows to save it in various formats out of \code{R} environment.
#'
#' @param tree A phylogenetic tree to use for generating the characters.
#' @param characters The number of morphological characters.
#' @param model Either an implemented (\code{"ER"} or \code{"HKY"}; see details) or user defined model for generating discrete morphological characters.
#' @param states A string of probabilities for the number of states for each characters (\code{default = 1}; i.e. 100% binary state characters; see details).
#' @param rates A function an it's parameters for the rates distribution (see details).
#' @param substitution A function an it's parameters for the substitutions distribution (see details).
#' @param ... Any optional arguments to be passed to the model argument.
#' @param CI A threshold consistency index value.
#' @param inapplicable Optional, a vector of characters inapplicability source (either \code{"character"} or \code{"clade"}; see details). The length of this vector must be at maximum half the total number of characters.
#' @param output Optional, an output file name for writing the matrix out of the \code{R} environement in \code{nexus} format.
#'
#' @details
#' \itemize {
#' 
#' \item The \code{model} arguments must be either a user's defined function for generating the discrete morphological characters (that intakes the states, rates and substitution arguments) or one of the two following:
#'      \itemize{
#'          \item \code{"ER"} uses the \code{ape::rTraitDisc} function with the \code{"ER"} model argument (= Mk model).
#'          \item \code{"HKY"} uses the \code{phyclust::gen.seq.HKY} function with \code{kappa = 2}, \code{pi = runif(4)} (divided by \code{sum(runif(4))}), \code{rate.scale = 1} and \code{L} being the number of \code{characters} and transforms the purines (A, G) into 0 and the pyrimidines (C, T) into 1.
#'      }
#'
#' \item The \code{states} argument attributes a number of states to each character by using the given probability vector for each number of states starting from 2.
#' For example \code{states = c(0.7, 0.2, 0.1)} will generate 70% of characters with 2 states, 20% of characters with 3 states and 10% of characters with 4 states. 
#' 
#' \item The \code{rates} and \code{substitution} arguments attributes a distribution function and it's optional parameters to a model. For example \code{rates = c(runif, 1, 10)} attributes a uniform distribution between 1 and 10 for the rates distribution.
#' 
#' \item The \code{inapplicable} argument intakes a vector of character inapplicability source rendering a number of characters inapplicable using the following sources:
#'      \itemize{
#'          \item \code{"character"} draws inapplicable characters directly from the character matrix, ignoring the phylogeny (i.e. for a random character X, an other random character Y will have inappicable characters fro each character states 0 for character X).
#'          \item \code{"clade"} draws inapplicable characters from the phylogeny: it will randomly apply inapplicable characters states for some characters by randomly selecting clades from the provided tree.
#'      }
#' For example \code{inapplicable = c(rep("character", 2), rep("clade", 2))} will generate 4 characters with inapplicable data, two using previous characters and two other using random clades.
#' }
#' 
#' @author Thomas Guillerme
#' @export


make.matrix <- function(tree, characters, states = 1, model = "ER", rates, substitution, ..., CI, inapplicables, output) {

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
        if(model == "ER") model <- rTraitDisc.mk
        if(model == "HKY") model <- gen.seq.HKY.binary
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

    #inapplicable
    if(!missing(inapplicable)) {
        check.class(inapplicable, "character")
        inap.source_options <- c("character", "clade")
        if(all(is.na(match(unique(inapplicable), inap.source_options)))) {
            stop("inapplicable argument must be a vector containing at least one of the following: ", paste(inap.source_options, collapse=", "), sep="")
        }        
        if(inapplicable > characters/2) {
            stop("Only half the number of characters can be inapplicable.")
        }
    }

    #output
    if(!missing(output)) {
        check.class(output, "character")
        check.length(output, 1, " must be a single string of characters.")
    }


    rTraitDisc.mk <- function(characters, tree, states, rates, model, ...) {
        replicate(characters, rTraitDisc(tree, k = k.sampler(states), rate = sample.distribution(1, rates), model = "ER", states = seq(from=0, to=(length(states))), ...))
    }

    #GENERATING THE CHARACTERS
    matrix <- replicate(characters, rTraitDisc(tree, k = k.sampler(states), rate = rates(1, ...), model = model, states = seq(from=0, to=(length(states)))))
    #matrix <- replicate(characters, rTraitDisc(tree, k = k.sampler(states), rate = rates(1,    ), model = model, states = seq(from=0, to=(length(states))))) ; warning("DEBUG MODE")

    #Repeat the invariant characters sampling
    while(any(apply(matrix, 2, is.invariant))) {
        matrix[, which(apply(matrix, 2, is.invariant) == TRUE)] <- replicate(length(which(apply(matrix, 2, is.invariant) == TRUE)), rTraitDisc(tree, k = k.sampler(states), rate = rates(1, ...), model = model, states = seq(from=0, to=(length(states)))))
        #matrix[, which(apply(matrix, 2, is.invariant) == TRUE)] <- replicate(length(which(apply(matrix, 2, is.invariant) == TRUE)), rTraitDisc(tree, k = k.sampler(states), rate = rates(1,    ), model = model, states = seq(from=0, to=(length(states))))) ; warning("DEBUG MODE")
    }

    #Add inapplicable characters (from character)
    if(!missing(inapplicable) & inap.source == "character") {
        #Select all the pairs of characters
        characters1 <- seq(from=1, to=characters, by=2)
        characters2 <- seq(from=2, to=characters, by=2)

        #Take a percentage of the shortest vector and use them for "inapplicability"
        if(length(characters2) < length(characters1)) {
            inapplicators <- characters2[1:inapplicable]
        } else {
            inapplicators <- characters1[1:inapplicable]
        }

        #making characters inapplicable
        for(inap in 1:inapplicable) {
            matrix[which(matrix[,inapplicators[inap]] == "0"), inapplicators[inap]+1] <- "-"
        }
    }

    #Add inapplicable characters (from clades)
    if(!missing(inapplicable) & inap.source == "clade") {
        #Select some clades
        clades <- replicate(inapplicable, select.clade(tree), simplify = FALSE)
        #Make them inapplicable in the matrix
        for(inap in 1:inapplicable) {
            matrix[clades[[inap]], inap] <- "-"
        }
    }


    #OUTPUTING THE FILE
    if(!missing(output)) {
        if(output == "nexus") {
            #Modifying write.nexus.data entry
            write.nexus.data.tmp <- write.nexus.data
            body(write.nexus.data.tmp)[[2]] <- quote(format <- match.arg(toupper(format), c("DNA", "PROTEIN", "STANDARD")))
            #Saving the file as a nexus
            write.nexus.data.tmp(matrix, file=paste(output.name, ".nex", sep=""), format = "standard")
            #Verbose
            cat(paste("Matrix saved as ", output.name, ".nex in current directory.\n", sep=""))
        }
    }
    
    return(matrix)

}