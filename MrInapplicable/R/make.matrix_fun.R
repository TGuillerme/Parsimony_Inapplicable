#sampling from a distribution
sample.distribution <- function(n,args) {
    fun <- args[[1]]
    args[[1]] <- n
    return(do.call(fun, args))
}

#scaling the results from a distribution sample (to be equal to 1)
proportional.distribution <- function(n, distribution, ...) {
    freq <- distribution(n, ...)
    return(freq/sum(freq))
}

#seqgen HKY binary
gen.seq.HKY.binary <- function(tree, characters) {
    matrix <- phyclust::gen.seq.HKY(tree, pi = proportional.distribution(4, runif), kappa = 2, L = characters, rate.scale = 1)
    #Removing the phylip header and separating the OTUs
    matrix <- apply(matrix(as.matrix(matrix)[-1,]), 1, strsplit, split = " ")
    #Extracting the data
    data <- matrix(unlist(apply(matrix(rapply(lapply(matrix, `[[`, 1), function(x) tail(x, 1))), 1, strsplit, split="")), ncol = characters, byrow = TRUE)

    #Transforming the base pairs
    #Purines
    data <- ifelse(data == "A" , "0", data)
    data <- ifelse(data == "G" , "0", data)
    #Pyrimidines
    data <- ifelse(data == "C" , "1", data)
    data <- ifelse(data == "T" , "1", data)

    #Adding row names
    rownames(data) <- lapply(lapply(matrix, `[[`, 1), `[[`, 1)

    return(data)
}



#sampling the number of characters states
k.sampler <- function(states) {
    if(length(states) == 1) {
        #Only binary characters
        return(2)
    } else {
        return(sample(2:(length(states)+1), 1, prob = states))
    }
}

#Invariant characters detector
is.invariant <- function(character) {
    if(length(unique(character)) == 1) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#Generating inapplicable characters
select.clade <- function(tree) {
    return(extract.clade(tree, node = sample(1:Nnode(tree), 1)+Ntip(tree))$tip.label)
}

# 