## Loading the packages

if(!require(devtools)) install.packages("devtools")
install_github("TGuillerme/Parsimony_Inapplicable/IterativeAlgo")
library(IterativeAlgo)
if(!require(ape)) install.packages("ape")
if(!require(phangorn)) install.packages("phangorn")
if(!require(phyclust)) install.packages("phyclust")
if(!require(Claddis)){
    install_github("graemetlloyd/Claddis")
    library(Claddis)
}


## Generating the random tree
my_tree <- rcoal(500)

## Setting up the parameters
## Evolutionary model rates (gamma distribution with alpha = 0.5)
my_rates = c(rgamma, 0.5, 1)
## States: 85% binary and 15% 3 states
my_states = c(0.85, 0.15)
## Inapplicables: 10% inapplicable characters states based on characters and 10% based on clades
my_inapplicables <- c(rep("character", 100), rep("clade", 100))

## Generating the random matrix
morpho_mat <- make.matrix(my_tree, characters = 1000, model = "ER", rates = my_rates, states = my_states, invariant = FALSE, verbose = TRUE)

## Adding the Inapplicables
morpho_mat <- apply.inapplicables(morpho_mat, inapplicable = my_inapplicables, my_tree, invariant = TRUE, verbose = TRUE)

## Checking the matrix status
check.matrix(morpho_mat, my_tree)

## Adding the missing data (30%)
morpho_mat[sample(1:500000, 150000)] <- "?"

## ape::write.nexus.data() quick fix
write.morpho.nexus <- write.nexus.data
body(write.morpho.nexus)[[2]] <- quote(format <- match.arg(toupper(format), c("DNA", "PROTEIN", "STANDARD")))

## Write the matrix
write.morpho.nexus(morpho_mat, file = "500t_1000c.nex", interleaved = FALSE, format = "STANDARD")
