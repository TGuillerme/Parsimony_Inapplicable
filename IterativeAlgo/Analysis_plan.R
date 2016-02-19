library(devtools)
install("../IterativeAlgo")
document()
test()
library(IterativeAlgo)
library(phangorn)
#Analysis plan


# Step 1: get a matrix
#matrix <- ReadMorphNexus("empirical.nex")

# Generating a realistic birth-death tree
# Fancy way
# birth <- function(t) 1/(1 + exp(0.2*t - 1)) #TG: define the birth process distribution (here logistic)
# death <- 0.07 #TG: define the death process distribution
# tree <- rbdtree(b, mu, Tmax = 2) # Prefer diversitree::tree.bd





parameters_names <- 
size_names <- c("10*50", "10*100", )



# Random (non-biological) way:
tree <- rcoal(10)
## setting up the parameters
my_rates = c(rgamma, 0.5, 1) # A gamma rate distribution with of shape alpha = 5
my_substitutions = c(runif, 2, 2) # A fixed substitution rate of 2 (T/T ratio in HKY)

## A pure Mk matrix (10*50)
matrixMk <- make.matrix(tree, characters = 50, model = "ER", rates = my_rates, substitution = my_substitutions, verbose = TRUE)
## HKY binary (10*50)
matrixHKY <- make.matrix(tree, characters = 50, model = "HKY", rates = my_rates, substitution = my_substitutions, , verbose = TRUE)


#Generating a bunch of coalescent trees

trees <- replicate(100, rcoal(15), simplify = FALSE)
class(trees) <- "multiPhylo"
my_rates = c(rgamma, 0.5, 1)
my_substitutions = c(runif, 2, 2)

matrices_Mk <- lapply(trees, make.matrix, characters = 100, model = "ER", rates = my_rates, substitution = my_substitutions, verbose = TRUE)
matrices_HKY <- lapply(trees, make.matrix, characters = 100, model = "HKY", rates = my_rates, substitution = my_substitutions, verbose = TRUE)

check_Mk <- mapply(check.matrix, matrices_Mk, trees)
check_HKY <- mapply(check.matrix, matrices_HKY, trees)





ylim_pars <- c(min(c(check_Mk[1,], check_HKY[1,])), max(c(check_Mk[1,], check_HKY[1,])))
ylim_CI <- c(min(c(check_Mk[2,], check_HKY[2,])), max(c(check_Mk[2,], check_HKY[2,])))
ylim_RI <- c(min(c(check_Mk[3,], check_HKY[3,])), max(c(check_Mk[3,], check_HKY[3,])))
ylim_RF <- c(min(c(check_Mk[4,], check_HKY[4,])), max(c(check_Mk[4,], check_HKY[4,])))
xlim <- c(0, 100)

pdf(paste("~/Projects/Parsimony_Inapplicable/IterativeAlgo/Testings_parameters/", size_names,"-",parameters_names, ".pdf", sep=""), width = 11.69, height = 8.27)

par(mfrow = c(4,1), bty = "n")

plot(1:100, sort(check_Mk[1,]), ylim = ylim_pars , xlab = "", ylab = "Parsimony Score", type = "l",
    main = paste(size_names, "-", parameters_names), xlim = xlim)
legend(min(xlim), max(ylim_pars), col = c("black", "red"), legend = c("Mk", "HKY.binary"), lty = c(1,1), cex = 0.8)
lines(1:100, sort(check_HKY[1,]), type = "l", col = "red")

plot(1:100, sort(check_Mk[2,]), ylim = ylim_CI , xlab = "", ylab = "Consistency Index", type = "l", xlim = xlim)
lines(1:100, sort(check_HKY[2,]), type = "l", col = "red")

plot(1:100, sort(check_Mk[3,]), ylim = ylim_RI , xlab = "", ylab = "Retention Index", type = "l", xlim = xlim)
lines(1:100, sort(check_HKY[3,]), type = "l", col = "red")

plot(1:100, sort(check_Mk[4,]), ylim = ylim_RF , xlab = "Iteration (ordered)", ylab = "RF distance", type = "l", xlim = xlim)
lines(1:100, sort(check_HKY[4,]), type = "l", col = "red")

dev.off()











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