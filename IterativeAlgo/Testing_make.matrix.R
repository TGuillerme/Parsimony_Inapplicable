library(devtools)
install("../IterativeAlgo")
library(IterativeAlgo)
library(phangorn)
#Analysis plan

taxa_param <- c(10, 20 ,30)
size_param <- c(50, 100, 150)
rates_param <- list(c(rgamma, 0.1, 1), c(rgamma, 1, 1), c(rgamma, 10, 1), c(runif, 0.01, 10))
rates_param_names <- c("G0.1","G1","G10","U0.01-10")
substitution_param <- list(c(runif, 1, 1), c(runif, 2, 2))
substitution_param_names <- c("TT1", "TT2")
iterations <- 100

for(taxa in 1:length(taxa_param)) {
    for(size in 1:length(size_param)) {
        for(rate in 1:length(rates_param)) {
            for(substitution in 1:length(substitution_param)) {

                test_name <- paste(taxa_param[taxa], size_param[size], rates_param_names[rate], substitution_param_names[substitution], sep="_")

                trees <- replicate(iterations, rcoal(taxa_param[taxa]), simplify = FALSE)
                my_rates = rates_param[[rate]]
                my_substitutions = substitution_param[[substitution]]

                matrices_Mk <- lapply(trees, make.matrix, characters = size_param[size], model = "ER", rates = my_rates, substitution = my_substitutions, verbose = TRUE)
                matrices_HKY <- lapply(trees, make.matrix, characters = size_param[size], model = "HKY", rates = my_rates, substitution = my_substitutions, verbose = TRUE)

                check_Mk <- mapply(check.matrix, matrices_Mk, trees)
                check_HKY <- mapply(check.matrix, matrices_HKY, trees)

                ylim_pars <- c(min(c(check_Mk[1,], check_HKY[1,])), max(c(check_Mk[1,], check_HKY[1,])))
                ylim_CI <- c(min(c(check_Mk[2,], check_HKY[2,])), max(c(check_Mk[2,], check_HKY[2,])))
                ylim_RI <- c(min(c(check_Mk[3,], check_HKY[3,])), max(c(check_Mk[3,], check_HKY[3,])))
                ylim_RF <- c(min(c(check_Mk[4,], check_HKY[4,])), max(c(check_Mk[4,], check_HKY[4,])))
                xlim <- c(0, iterations)

                pdf(paste("~/Projects/Parsimony_Inapplicable/IterativeAlgo/Testings_parameters/", test_name, ".pdf", sep=""), width = 11.69, height = 8.27)

                par(mfrow = c(4,1), bty = "n")

                plot(1:iterations, sort(check_Mk[1,]), ylim = ylim_pars , xlab = "", ylab = "Parsimony Score", type = "l",
                    main = test_name, xlim = xlim)
                legend(min(xlim), max(ylim_pars), col = c("black", "red"), legend = c("Mk", "HKY.binary"), lty = c(1,1), cex = 0.8)
                lines(1:iterations, sort(check_HKY[1,]), type = "l", col = "red")

                plot(1:iterations, sort(check_Mk[2,]), ylim = ylim_CI , xlab = "", ylab = "Consistency Index", type = "l", xlim = xlim)
                lines(1:iterations, sort(check_HKY[2,]), type = "l", col = "red")

                plot(1:iterations, sort(check_Mk[3,]), ylim = ylim_RI , xlab = "", ylab = "Retention Index", type = "l", xlim = xlim)
                lines(1:iterations, sort(check_HKY[3,]), type = "l", col = "red")

                plot(1:iterations, sort(check_Mk[4,]), ylim = ylim_RF , xlab = "Iteration (ordered)", ylab = "RF distance", type = "l", xlim = xlim)
                lines(1:iterations, sort(check_HKY[4,]), type = "l", col = "red")

                dev.off()

                savings <- list(list(check_Mk, check_HKY), list(matrices_Mk, matrices_HKY))
                save(savings, file = paste("~/Projects/Parsimony_Inapplicable/IterativeAlgo/Testings_parameters/", test_name, ".rda"))

            }
        }
    }
}





#################
# A good parameters combination seems to be:
my_rates = c(rgamma, 1, 0.1)
my_substitutions = c(runif, 2, 2)

# With ER being better than HKY
my_model = "ER"

# And states proportions of
my_states = c(0.85, 0.15)
#################


# Now checking the importance of the input tree (between rcoal, rtree and tree.bd)




tree <- rtree.bd(15)


tree.bd_pars <- NULL
tree.bd_CI <- NULL
tree.bd_RI <- NULL
tree.bd_RF <- NULL
for(i in 1:100) {
    # Proper birth death tree
    tree <- tree.bd(rand.birth.death(), max.taxa = 15) # Make sure is not null!
    while(is.null(tree)) {
        tree <- tree.bd(rand.birth.death(), max.taxa = 15) # Make sure is not null!
    }

    # A pure Mk matrix (10*50)
    matrixMk <- make.matrix(tree, characters = 150, model = my_model, states = my_states, rates = my_rates, substitution = my_substitutions, verbose = TRUE)
    tmp <- check.matrix(matrixMk, tree)
    tree.bd_pars[i] <- tmp[1,1]
    tree.bd_CI[i] <- tmp[2,1]
    tree.bd_RI[i] <- tmp[3,1]
    tree.bd_RF[i] <- tmp[4,1]
}


rcoal_pars <- NULL
rcoal_CI <- NULL
rcoal_RI <- NULL
rcoal_RF <- NULL
for(i in 1:100) {
    # Proper birth death tree
    tree <- rcoal(15)

    # setting up the parameters
    my_rates = c(rgamma, 1, 0.1) # A gamma rate distribution with of shape alpha = 5
    my_substitutions = c(runif, 2, 2) # A fixed substitution rate of 2 (T/T ratio in HKY)

    # A pure Mk matrix (10*50)
    matrixMk <- make.matrix(tree, characters = 150, model = my_model, states = my_states, rates = my_rates, substitution = my_substitutions, verbose = TRUE)
    tmp <- check.matrix(matrixMk, tree)
    rcoal_pars[i] <- tmp[1,1]
    rcoal_CI[i] <- tmp[2,1]
    rcoal_RI[i] <- tmp[3,1]
    rcoal_RF[i] <- tmp[4,1]
}
