library(IterativeAlgo)
library(phangorn)

#################
# TESTING THE TREE EFFECT
#################

#################
# A good parameters combination seems to be:
my_rates = c(rgamma, 1, 0.1)
my_substitutions = c(runif, 2, 2)

# With ER being better than HKY
my_model = "ER"

# And states proportions of
my_states = c(0.85, 0.15)
#################

# Creating the list for storing the results
test_list_tmp <- list("Parsimony" = NULL, "CI" = NULL, "RI" = NULL, "RF" = NULL)
parameters_test <- list("rtree" = test_list_tmp, "rcoa" = test_list_tmp, "rtree.bd" = test_list_tmp)


#Doing 100 random replicates for each tree
taxa = 15 ; nchar = 150
for (rep in 1:1000) {

    #Gnerate the trees
    tree_rtree <- rtree(taxa)
    tree_rcoal <- rcoal(taxa)
    tree_rtree.bd <- rtree.bd(taxa)

    #Generate the matrices
    matrix_rtree <- make.matrix(tree_rtree, characters = nchar, model = my_model, states = my_states, rates = my_rates, substitution = my_substitutions, verbose = FALSE)
    matrix_rcoal <- make.matrix(tree_rcoal, characters = nchar, model = my_model, states = my_states, rates = my_rates, substitution = my_substitutions, verbose = FALSE)
    matrix_rtree.bd <- make.matrix(tree_rtree.bd, characters = nchar, model = my_model, states = my_states, rates = my_rates, substitution = my_substitutions, verbose = FALSE)

    #Checking the matrices
    tmp_rtree <- check.matrix(matrix_rtree, tree_rtree)
    tmp_rcoal <- check.matrix(matrix_rcoal, tree_rcoal)
    tmp_rtree.bd <- check.matrix(matrix_rtree.bd, tree_rtree.bd)

    #Storing the results
    #Parsimony
    parameters_test[[1]][[1]] <- c(parameters_test[[1]][[1]], tmp_rtree[1,1])
    parameters_test[[2]][[1]] <- c(parameters_test[[2]][[1]], tmp_rcoal[1,1])
    parameters_test[[3]][[1]] <- c(parameters_test[[3]][[1]], tmp_rtree.bd[1,1])
    #CI
    parameters_test[[1]][[2]] <- c(parameters_test[[1]][[2]], tmp_rtree[2,1])
    parameters_test[[2]][[2]] <- c(parameters_test[[2]][[2]], tmp_rcoal[2,1])
    parameters_test[[3]][[2]] <- c(parameters_test[[3]][[2]], tmp_rtree.bd[2,1])
    #RI
    parameters_test[[1]][[3]] <- c(parameters_test[[1]][[3]], tmp_rtree[3,1])
    parameters_test[[2]][[3]] <- c(parameters_test[[2]][[3]], tmp_rcoal[3,1])
    parameters_test[[3]][[3]] <- c(parameters_test[[3]][[3]], tmp_rtree.bd[3,1])
    #RF
    parameters_test[[1]][[4]] <- c(parameters_test[[1]][[4]], tmp_rtree[4,1])
    parameters_test[[2]][[4]] <- c(parameters_test[[2]][[4]], tmp_rcoal[4,1])
    parameters_test[[3]][[4]] <- c(parameters_test[[3]][[4]], tmp_rtree.bd[4,1])
}


op <- par(mfrow = c(2,2), bty = "n")

xvar <- summary(c(density(parameters_test[[1]][[1]])$x, density(parameters_test[[2]][[1]])$x, density(parameters_test[[3]][[1]])$x))
yvar <- summary(c(density(parameters_test[[1]][[1]])$y, density(parameters_test[[2]][[1]])$y, density(parameters_test[[3]][[1]])$y))
plot(density(parameters_test[[1]][[1]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "Parsimony score", main = "")
lines(density(parameters_test[[2]][[1]]), col = "red")
lines(density(parameters_test[[3]][[1]]), col = "blue")
legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Yule", "Coalescent", "Birth-Death"), lty = c(1,1,1))

xvar <- summary(c(density(parameters_test[[1]][[2]])$x, density(parameters_test[[2]][[2]])$x, density(parameters_test[[3]][[2]])$x))
yvar <- summary(c(density(parameters_test[[1]][[2]])$y, density(parameters_test[[2]][[2]])$y, density(parameters_test[[3]][[2]])$y))
plot(density(parameters_test[[1]][[2]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "Consistency index", main = "")
lines(density(parameters_test[[2]][[2]]), col = "red")
lines(density(parameters_test[[3]][[2]]), col = "blue")
#legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Yule", "Coalescent", "Birth-Death"), lty = c(1,1,1))

xvar <- summary(c(density(parameters_test[[1]][[3]])$x, density(parameters_test[[2]][[3]])$x, density(parameters_test[[3]][[3]])$x))
yvar <- summary(c(density(parameters_test[[1]][[3]])$y, density(parameters_test[[2]][[3]])$y, density(parameters_test[[3]][[3]])$y))
plot(density(parameters_test[[1]][[3]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "Retention index", main = "")
lines(density(parameters_test[[2]][[3]]), col = "red")
lines(density(parameters_test[[3]][[3]]), col = "blue")
#legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Yule", "Coalescent", "Birth-Death"), lty = c(1,1,1))

xvar <- summary(c(density(parameters_test[[1]][[4]])$x, density(parameters_test[[2]][[4]])$x, density(parameters_test[[3]][[4]])$x))
yvar <- summary(c(density(parameters_test[[1]][[4]])$y, density(parameters_test[[2]][[4]])$y, density(parameters_test[[3]][[4]])$y))
plot(density(parameters_test[[1]][[4]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "RF distance", main = "")
lines(density(parameters_test[[2]][[4]]), col = "red")
lines(density(parameters_test[[3]][[4]]), col = "blue")
#legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Yule", "Coalescent", "Birth-Death"), lty = c(1,1,1))

par(op)


#################
# TESTING THE MODEL EFFECT
#################

#################
# A good parameters combination seems to be:
my_rates = c(rgamma, 1, 0.1)
my_substitutions = c(runif, 2, 2)

# And states proportions of
my_states = 1
#################

# Creating the list for storing the results
test_list_tmp <- list("Parsimony" = NULL, "CI" = NULL, "RI" = NULL, "RF" = NULL)
parameters_test <- list("Mk" = test_list_tmp, "HKY" = test_list_tmp, "MKvar" = test_list_tmp)


#Doing 100 random replicates for each tree
taxa = 15 ; nchar = 150
for (rep in 1:1000) {

    #Gnerate the trees
    tree_rtree <- rtree(taxa)

    #Generate the matrices
    matrix_Mk <- make.matrix(tree_rtree, characters = nchar, model = "ER", states = 1, rates = my_rates, substitution = my_substitutions, verbose = FALSE)
    matrix_HKY <- make.matrix(tree_rcoal, characters = nchar, model = "HKY", states = 1, rates = my_rates, substitution = my_substitutions, verbose = FALSE)
    matrix_Mkvar <- make.matrix(tree_rtree, characters = nchar, model = "ER", states = c(0.85, 0.15), rates = my_rates, substitution = my_substitutions, verbose = FALSE)

    #Checking the matrices
    tmp_Mk <- check.matrix(matrix_Mk, tree_rtree)
    tmp_HKY <- check.matrix(matrix_HKY, tree_rtree)
    tmp_Mkvar <- check.matrix(matrix_Mkvar, tree_rtree)

    #Storing the results
    #Parsimony
    parameters_test[[1]][[1]] <- c(parameters_test[[1]][[1]], tmp_Mk[1,1])
    parameters_test[[2]][[1]] <- c(parameters_test[[2]][[1]], tmp_HKY[1,1])
    parameters_test[[3]][[1]] <- c(parameters_test[[3]][[1]], tmp_Mkvar[1,1])
    #CI
    parameters_test[[1]][[2]] <- c(parameters_test[[1]][[2]], tmp_Mk[2,1])
    parameters_test[[2]][[2]] <- c(parameters_test[[2]][[2]], tmp_HKY[2,1])
    parameters_test[[3]][[2]] <- c(parameters_test[[3]][[2]], tmp_Mkvar[2,1])
    #RI
    parameters_test[[1]][[3]] <- c(parameters_test[[1]][[3]], tmp_Mk[3,1])
    parameters_test[[2]][[3]] <- c(parameters_test[[2]][[3]], tmp_HKY[3,1])
    parameters_test[[3]][[3]] <- c(parameters_test[[3]][[3]], tmp_Mkvar[3,1])
    #RF
    parameters_test[[1]][[4]] <- c(parameters_test[[1]][[4]], tmp_Mk[4,1])
    parameters_test[[2]][[4]] <- c(parameters_test[[2]][[4]], tmp_HKY[4,1])
    parameters_test[[3]][[4]] <- c(parameters_test[[3]][[4]], tmp_Mkvar[4,1])
}


op <- par(mfrow = c(2,2), bty = "n")

xvar <- summary(c(density(parameters_test[[1]][[1]])$x, density(parameters_test[[2]][[1]])$x, density(parameters_test[[3]][[1]])$x))
yvar <- summary(c(density(parameters_test[[1]][[1]])$y, density(parameters_test[[2]][[1]])$y, density(parameters_test[[3]][[1]])$y))
plot(density(parameters_test[[1]][[1]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "Parsimony score", main = "")
lines(density(parameters_test[[2]][[1]]), col = "red")
lines(density(parameters_test[[3]][[1]]), col = "blue")
legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Mk", "HKY", "Mk-var"), lty = c(1,1,1))

xvar <- summary(c(density(parameters_test[[1]][[2]])$x, density(parameters_test[[2]][[2]])$x, density(parameters_test[[3]][[2]])$x))
yvar <- summary(c(density(parameters_test[[1]][[2]])$y, density(parameters_test[[2]][[2]])$y, density(parameters_test[[3]][[2]])$y))
plot(density(parameters_test[[1]][[2]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "Consistency index", main = "")
lines(density(parameters_test[[2]][[2]]), col = "red")
lines(density(parameters_test[[3]][[2]]), col = "blue")
#legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Yule", "Coalescent", "Birth-Death"), lty = c(1,1,1))

xvar <- summary(c(density(parameters_test[[1]][[3]])$x, density(parameters_test[[2]][[3]])$x, density(parameters_test[[3]][[3]])$x))
yvar <- summary(c(density(parameters_test[[1]][[3]])$y, density(parameters_test[[2]][[3]])$y, density(parameters_test[[3]][[3]])$y))
plot(density(parameters_test[[1]][[3]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "Retention index", main = "")
lines(density(parameters_test[[2]][[3]]), col = "red")
lines(density(parameters_test[[3]][[3]]), col = "blue")
#legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Yule", "Coalescent", "Birth-Death"), lty = c(1,1,1))

xvar <- summary(c(density(parameters_test[[1]][[4]])$x, density(parameters_test[[2]][[4]])$x, density(parameters_test[[3]][[4]])$x))
yvar <- summary(c(density(parameters_test[[1]][[4]])$y, density(parameters_test[[2]][[4]])$y, density(parameters_test[[3]][[4]])$y))
plot(density(parameters_test[[1]][[4]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "RF distance", main = "")
lines(density(parameters_test[[2]][[4]]), col = "red")
lines(density(parameters_test[[3]][[4]]), col = "blue")
#legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Yule", "Coalescent", "Birth-Death"), lty = c(1,1,1))

par(op)


#################
# TESTING THE RATES DISTRIBUTION SHAPE
#################

#################
# A good parameters combination seems to be:
my_substitutions = c(runif, 2, 2)

# With ER being better than HKY
my_model = "ER"

# And states proportions of
my_states = c(0.85, 0.15)
#################

# Creating the list for storing the results
test_list_tmp <- list("Parsimony" = NULL, "CI" = NULL, "RI" = NULL, "RF" = NULL)
parameters_test <- list("gamma" = test_list_tmp, "norm" = test_list_tmp, "unif" = test_list_tmp)


#Doing 100 random replicates for each tree
taxa = 15 ; nchar = 150
for (rep in 1:1000) {

    #Gnerate the trees
    tree_rtree <- rtree(taxa)

    #Generate the matrices
    matrix_gamma <- make.matrix(tree_rtree, characters = nchar, model = my_model, states = my_states, rates = c(rgamma, 1, 1), substitution = my_substitutions, verbose = FALSE)
    matrix_norm <- make.matrix(tree_rtree, characters = nchar, model = my_model, states = my_states, rates = c(rlnorm, 0, 1), substitution = my_substitutions, verbose = FALSE)
    matrix_unif <- make.matrix(tree_rtree, characters = nchar, model = my_model, states = my_states, rates = c(runif, 0.01, 10), substitution = my_substitutions, verbose = FALSE)

    #Checking the matrices
    tmp_gamma <- check.matrix(matrix_gamma, tree_rtree)
    tmp_norm <- check.matrix(matrix_norm, tree_rtree)
    tmp_unif <- check.matrix(matrix_unif, tree_rtree)

    #Storing the results
    #Parsimony
    parameters_test[[1]][[1]] <- c(parameters_test[[1]][[1]], tmp_gamma[1,1])
    parameters_test[[2]][[1]] <- c(parameters_test[[2]][[1]], tmp_norm[1,1])
    parameters_test[[3]][[1]] <- c(parameters_test[[3]][[1]], tmp_unif[1,1])
    #CI
    parameters_test[[1]][[2]] <- c(parameters_test[[1]][[2]], tmp_gamma[2,1])
    parameters_test[[2]][[2]] <- c(parameters_test[[2]][[2]], tmp_norm[2,1])
    parameters_test[[3]][[2]] <- c(parameters_test[[3]][[2]], tmp_unif[2,1])
    #RI
    parameters_test[[1]][[3]] <- c(parameters_test[[1]][[3]], tmp_gamma[3,1])
    parameters_test[[2]][[3]] <- c(parameters_test[[2]][[3]], tmp_norm[3,1])
    parameters_test[[3]][[3]] <- c(parameters_test[[3]][[3]], tmp_unif[3,1])
    #RF
    parameters_test[[1]][[4]] <- c(parameters_test[[1]][[4]], tmp_gamma[4,1])
    parameters_test[[2]][[4]] <- c(parameters_test[[2]][[4]], tmp_norm[4,1])
    parameters_test[[3]][[4]] <- c(parameters_test[[3]][[4]], tmp_unif[4,1])
}


op <- par(mfrow = c(2,2), bty = "n")

xvar <- summary(c(density(parameters_test[[1]][[1]])$x, density(parameters_test[[2]][[1]])$x, density(parameters_test[[3]][[1]])$x))
yvar <- summary(c(density(parameters_test[[1]][[1]])$y, density(parameters_test[[2]][[1]])$y, density(parameters_test[[3]][[1]])$y))
plot(density(parameters_test[[1]][[1]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "Parsimony score", main = "")
lines(density(parameters_test[[2]][[1]]), col = "red")
lines(density(parameters_test[[3]][[1]]), col = "blue")
legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Gamma", "Log-normal", "Uniform"), lty = c(1,1,1))

xvar <- summary(c(density(parameters_test[[1]][[2]])$x, density(parameters_test[[2]][[2]])$x, density(parameters_test[[3]][[2]])$x))
yvar <- summary(c(density(parameters_test[[1]][[2]])$y, density(parameters_test[[2]][[2]])$y, density(parameters_test[[3]][[2]])$y))
plot(density(parameters_test[[1]][[2]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "Consistency index", main = "")
lines(density(parameters_test[[2]][[2]]), col = "red")
lines(density(parameters_test[[3]][[2]]), col = "blue")
#legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Yule", "Coalescent", "Birth-Death"), lty = c(1,1,1))

xvar <- summary(c(density(parameters_test[[1]][[3]])$x, density(parameters_test[[2]][[3]])$x, density(parameters_test[[3]][[3]])$x))
yvar <- summary(c(density(parameters_test[[1]][[3]])$y, density(parameters_test[[2]][[3]])$y, density(parameters_test[[3]][[3]])$y))
plot(density(parameters_test[[1]][[3]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "Retention index", main = "")
lines(density(parameters_test[[2]][[3]]), col = "red")
lines(density(parameters_test[[3]][[3]]), col = "blue")
#legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Yule", "Coalescent", "Birth-Death"), lty = c(1,1,1))

xvar <- summary(c(density(parameters_test[[1]][[4]])$x, density(parameters_test[[2]][[4]])$x, density(parameters_test[[3]][[4]])$x))
yvar <- summary(c(density(parameters_test[[1]][[4]])$y, density(parameters_test[[2]][[4]])$y, density(parameters_test[[3]][[4]])$y))
plot(density(parameters_test[[1]][[4]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "RF distance", main = "")
lines(density(parameters_test[[2]][[4]]), col = "red")
lines(density(parameters_test[[3]][[4]]), col = "blue")
#legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Yule", "Coalescent", "Birth-Death"), lty = c(1,1,1))

par(op)


#################
# TESTING THE RATES GAMMA SHAPE
#################

#################
# A good parameters combination seems to be:
my_substitutions = c(runif, 2, 2)

# With ER being better than HKY
my_model = "ER"

# And states proportions of
my_states = c(0.85, 0.15)
#################



# Now checking the importance of the input tree (between rcoal, rtree and tree.bd)

# Creating the list for storing the results
test_list_tmp <- list("Parsimony" = NULL, "CI" = NULL, "RI" = NULL, "RF" = NULL)
parameters_test <- list("(0.5,1)" = test_list_tmp, "(1,1)" = test_list_tmp, "(0.5,0.1)" = test_list_tmp, "(1,0.1)" = test_list_tmp)


#Doing 100 random replicates for each tree
taxa = 15 ; nchar = 150
for (rep in 1:1000) {

    #Gnerate the trees
    tree_rtree <- rtree(taxa)

    #Generate the matrices
    matrix_0.5_1 <- make.matrix(tree_rtree, characters = nchar, model = my_model, states = my_states, rates = c(rgamma, 0.5, 1), substitution = my_substitutions, verbose = FALSE)
    matrix_1_1 <- make.matrix(tree_rtree, characters = nchar, model = my_model, states = my_states, rates = c(rgamma, 1, 1), substitution = my_substitutions, verbose = FALSE)
    matrix_0.5_0.1 <- make.matrix(tree_rtree, characters = nchar, model = my_model, states = my_states, rates = c(rgamma, 0.5, 0.1), substitution = my_substitutions, verbose = FALSE)
    matrix_1_0.1 <- make.matrix(tree_rtree, characters = nchar, model = my_model, states = my_states, rates = c(rgamma, 1, 0.1), substitution = my_substitutions, verbose = FALSE)

    #Checking the matrices
    tmp_0.5_1 <- check.matrix(matrix_0.5_1, tree_rtree)
    tmp_1_1 <- check.matrix(matrix_1_1, tree_rtree)
    tmp_0.5_0.1 <- check.matrix(matrix_0.5_0.1, tree_rtree)
    tmp_1_0.1 <- check.matrix(matrix_1_0.1, tree_rtree)

    #Storing the results
    #Parsimony
    parameters_test[[1]][[1]] <- c(parameters_test[[1]][[1]] , tmp_0.5_1[1,1])
    parameters_test[[2]][[1]] <- c(parameters_test[[2]][[1]] , tmp_1_1[1,1])
    parameters_test[[3]][[1]] <- c(parameters_test[[3]][[1]] , tmp_0.5_0.1[1,1])
    parameters_test[[4]][[1]] <- c(parameters_test[[4]][[1]] , tmp_1_0.1[1,1])
    #CI
    parameters_test[[1]][[2]] <- c(parameters_test[[1]][[2]] , tmp_0.5_1[2,1])
    parameters_test[[2]][[2]] <- c(parameters_test[[2]][[2]] , tmp_1_1[2,1])
    parameters_test[[3]][[2]] <- c(parameters_test[[3]][[2]] , tmp_0.5_0.1[2,1])
    parameters_test[[4]][[2]] <- c(parameters_test[[4]][[2]] , tmp_1_0.1[2,1])
    #RI
    parameters_test[[1]][[3]] <- c(parameters_test[[1]][[3]] , tmp_0.5_1[3,1])
    parameters_test[[2]][[3]] <- c(parameters_test[[2]][[3]] , tmp_1_1[3,1])
    parameters_test[[3]][[3]] <- c(parameters_test[[3]][[3]] , tmp_0.5_0.1[3,1])
    parameters_test[[4]][[3]] <- c(parameters_test[[4]][[3]] , tmp_1_0.1[3,1])
    #RF
    parameters_test[[1]][[4]] <- c(parameters_test[[1]][[4]] , tmp_0.5_1[4,1])
    parameters_test[[2]][[4]] <- c(parameters_test[[2]][[4]] , tmp_1_1[4,1])
    parameters_test[[3]][[4]] <- c(parameters_test[[3]][[4]] , tmp_0.5_0.1[4,1])
    parameters_test[[4]][[4]] <- c(parameters_test[[4]][[4]] , tmp_1_0.1[4,1])
}


op <- par(mfrow = c(2,2), bty = "n")

xvar <- summary(c(density(parameters_test[[1]][[1]])$x, density(parameters_test[[2]][[1]])$x, density(parameters_test[[3]][[1]])$x, density(parameters_test[[4]][[1]])$x))
yvar <- summary(c(density(parameters_test[[1]][[1]])$y, density(parameters_test[[2]][[1]])$y, density(parameters_test[[3]][[1]])$y, density(parameters_test[[4]][[1]])$y))
plot(density(parameters_test[[1]][[1]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "Parsimony score", main = "")
lines(density(parameters_test[[2]][[1]]), col = "red")
lines(density(parameters_test[[3]][[1]]), col = "blue")
lines(density(parameters_test[[4]][[1]]), col = "green")
legend(xvar[1], yvar[6], col = palette()[1:6], legend = c("(0.5,1)", "(1,1)", "(0.5,0.1)", "(1,0.1)"), lty = c(rep(1, 4)))

xvar <- summary(c(density(parameters_test[[1]][[2]])$x, density(parameters_test[[2]][[2]])$x, density(parameters_test[[3]][[2]])$x, density(parameters_test[[4]][[2]])$x))
yvar <- summary(c(density(parameters_test[[1]][[2]])$y, density(parameters_test[[2]][[2]])$y, density(parameters_test[[3]][[2]])$y, density(parameters_test[[4]][[2]])$y))
plot(density(parameters_test[[1]][[2]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "Consistency index", main = "")
lines(density(parameters_test[[2]][[2]]), col = "red")
lines(density(parameters_test[[3]][[2]]), col = "blue")
lines(density(parameters_test[[4]][[2]]), col = "green")
#legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Yule", "Coalescent", "Birth-Death"), lty = c(1,1,1))

xvar <- summary(c(density(parameters_test[[1]][[3]])$x, density(parameters_test[[2]][[3]])$x, density(parameters_test[[3]][[3]])$x, density(parameters_test[[4]][[3]])$x))
yvar <- summary(c(density(parameters_test[[1]][[3]])$y, density(parameters_test[[2]][[3]])$y, density(parameters_test[[3]][[3]])$y, density(parameters_test[[4]][[3]])$y))
plot(density(parameters_test[[1]][[3]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "Retention index", main = "")
lines(density(parameters_test[[2]][[3]]), col = "red")
lines(density(parameters_test[[3]][[3]]), col = "blue")
lines(density(parameters_test[[4]][[3]]), col = "green")
#legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Yule", "Coalescent", "Birth-Death"), lty = c(1,1,1))

xvar <- summary(c(density(parameters_test[[1]][[4]])$x, density(parameters_test[[2]][[4]])$x, density(parameters_test[[3]][[4]])$x, density(parameters_test[[4]][[4]])$x))
yvar <- summary(c(density(parameters_test[[1]][[4]])$y, density(parameters_test[[2]][[4]])$y, density(parameters_test[[3]][[4]])$y, density(parameters_test[[4]][[4]])$y))
plot(density(parameters_test[[1]][[4]]), xlim = xvar[c(1,6)], ylim = yvar[c(1,6)], xlab = "RF distance", main = "")
lines(density(parameters_test[[2]][[4]]), col = "red")
lines(density(parameters_test[[3]][[4]]), col = "blue")
lines(density(parameters_test[[4]][[4]]), col = "green")
#legend(xvar[1], yvar[6], col = c("black", "red", "blue"), legend = c("Yule", "Coalescent", "Birth-Death"), lty = c(1,1,1))

par(op)