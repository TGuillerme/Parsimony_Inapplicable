#' @title Generates a morphological matrix.
#'
#' @description Generates a morphological matrix using \code{\link[ape]{rTraitDisc}} function with some optional inaplicable characters and allows to save it in various formats out of \code{R} environment.
#'
#' @param tree A phylogenetic tree to use for generating the characters.
#' @param characters The number of morphological characters.
#' @param states A string of probabilities for the number of states for each characters (see details).
#' @param inapplicable Optional, a number of inapplicable characters to generate.
#' @param inap.source Optional, where to draw the inapplicable characters from: \code{"character"} or \code{"clade"|.
#' @param rates A distribution from which to draw the characters rates.
#' @param ... Any optional arguments to be passed to rates.
#' @param model The characters evolution model.
#' @param output An optional output format among the following \code{"nexus"}, \code{"tnt"}, \code{"phylip"}.
#' @param output.name If output is not missing, the name of the file to write out of \code{R} environement.
#'
#' @details
#' For the \code{states} argument, the function attributes a number of states to each character by using the given probability vector for each number of states starting from 2.
#' For example \code{states = c(0.7, 0.2, 0.1)} will generate 70% of characters with 2 states, 20% of characters with 3 states and 10% of characters with 4 states.
#' For the \code{inap.source} argument, \code{"character"} option draws inapplicable characters directly from the character matrix, ignoring the phylogeny (i.e. for a random character X, an other random character Y will have inappicable characters fro each character states 0 for character X).
#' The \code{"clade"} option draws inapplicable characters from the phylogeny: it will randomly apply inapplicable characters states for some characters by randomly selecting clades from the provided tree.
#'
#' 
#' @author Thomas Guillerme
#' @export


make.matrix <- function(tree, characters, states = 1, inapplicable, inap.source = "character", rates = runif, ..., model = "ER", output, output.name) {
    #SANITIZNG
    #tree
    check.class(tree, "phylo")

    #characters
    check.class(characters, "numeric")
    check.length(characters, 1, " must be a single numeric value.")

    #states
    check.class(states, "numeric")
    if(sum(states) != 1) {
        stop("states argument must sum up to 1.")
    }

    #inapplicable
    if(!missing(inapplicable)) {
        check.class(inapplicable, "numeric")
        check.length(inapplicable, 1, " must be a single numeric value.")
        if(inapplicable > characters/2) {
            stop("Only half the number of characters can be inapplicable.")
        }
    }

    #inap.source
    if(!missing(inapplicable)) {
        inap.source_options <- c("character", "clade")
        check.class(inap.source, "character")
        check.length(inap.source, 1, paste(" must be one of the following: ", paste(inap.source_options, collapse=", "), sep=""))
        if(all(is.na(match(inap.source, inap.source_options)))) {
            stop("inap.source argument must be one of the following: ", paste(inap.source_options, collapse=", "), sep="")
        }
    }

    #rates
    check.class(rates, "function")
    # allow multiple distributions?

    #output
    if(!missing(output)) {
        check.class(output, "character", paste(" must be one of the following: ", paste(output_options, collapse=", "), sep=""))
        check.length(output, 1, paste(" must be one of the following: ", paste(output_options, collapse=", "), sep=""))
        output_options <- c("nexus") #add more options
        if(all(is.na(match(output, output_options)))) {
            stop("The output format must be one of the following: ", paste(output_options, collapse=", "), sep="")
        }
    }

    #output.name
    if(!missing(output)) {
        if(missing(output.name)) {
            stop("Output file name is missing!")
        } else {
            check.class(output.name, "character")
            check.length(output.name, 1, " must be a single string of characters.")
        }
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