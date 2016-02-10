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