library(devtools)
install("../Iterative_algo")
library(Iterative_algo)
library(Claddis)
install_github("cran/phylosim")
#Analysis plan


# Step 1: get a matrix
#matrix <- make.matrix(rtree(10), 50, states=c(0.85, 0.15), inapplicable=10, inap.source="character")
matrix <- ReadMorphNexus("empirical.nex")

# Remove characters with inapplicable tokens

# Separate inapplicable from applicable characters

detect.inapplicable <- function(character) {
    return(any(character == "-"))
}

#Inapplciable characters
matrix_inapplicable <- matrix[, apply(matrix, 2, detect.inapplicable)]
#All applicable characters
matrix_applicable <- matrix[, !apply(matrix, 2, detect.inapplicable)]

# Find patterns in the inapplicable characters




# Extracting inapplicable characters
get.inapplicable.pattern <- function(character) {
    return(character == "-")
}
# Matching a character to an other one
match.character <- function(character, match) {
    return(all(character == match))
}
# Matching characters to a repeated pattern
match.pattern <- function(duplicated_pattern, patterns, matrix_inapplicable) {
    # Getting the duplicated characters
    duplicated <- which(apply(patterns, 2, match.character, match=patterns[,duplicated_pattern]))

    # Getting the duplicated matrix
    return(list("matrix" = matrix_inapplicable[, duplicated], "characters" = duplicated))
}




# Finding inapplicable patterns
extract.inapplicable.pattern <- function(matrix_inapplicable) {

    # Extracting the patterns
    patterns <- apply(matrix_inapplicable, 2, get.inapplicable.pattern)

    if(any(duplicated(patterns, MARGIN = 2))) {
        # Extracting the duplicated patterns (if any)
        matrix_patterns <- lapply(as.list( which(duplicated(patterns, MARGIN = 2))), match.pattern, patterns, matrix_inapplicable)

        # Getting the matrices
        matrices <- lapply(matrix_patterns, `[[`, 1)
        characters <- lapply(matrix_patterns, `[[`, 2)
    }

    #
}

#mb file.nex