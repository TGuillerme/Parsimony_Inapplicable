#TEST make.matrix

context("make.matrix")

#Testing sample.distribution
test_that("sample.distribution works", {
    #errors
    expect_error(sample.distribution("a", c(runif,1,2)))
    expect_error(sample.distribution(1, "c(runif,1,2)"))
    expect_error(sample.distribution(1, c(aov,1,2)))

    #Returns the right number of values
    expect_equal(length(sample.distribution(1, c(runif))), 1)
    expect_equal(length(sample.distribution(1, c(runif, 1, 2))), 1)
    expect_equal(length(sample.distribution(1000, c(runif, 1, 2))), 1000)

    #Returns values in the range
    expect_equal(length(sample.distribution(1, c(runif))), 1)
    expect_less_than(max(sample.distribution(1000, c(runif, 1,2))), 2.0000000001)
    expect_more_than(min(sample.distribution(1000, c(runif, 1,2))), 0.9999999999)
})

#Testing proportional.distribution
test_that("proportional.distribution works", {
    #errors
    expect_error(proportional.distribution("a", runif))
    expect_error(proportional.distribution(4, "runif"))
    expect_error(proportional.distribution(4, runif, "a"))

    #sum(results) = 1
    expect_equal(sum(proportional.distribution(4, runif)), 1)
    expect_equal(sum(proportional.distribution(100, runif)), 1)
    expect_equal(sum(proportional.distribution(4, runif, 1000, 2000)), 1)
    expect_equal(sum(proportional.distribution(4, rnorm)), 1)
})

#Testing gen.seq.HKY.binary
test_that("gen.seq.HKY.binary works", {
    #errors
    expect_error(gen.seq.HKY.binary("a", c(runif, 2, 2), c(runif, 1, 1)))
    expect_error(gen.seq.HKY.binary(5, c(runif, 2, 2), c(runif, 1, 1)))
    expect_error(gen.seq.HKY.binary(rtree(5), runif, c(runif, 1, 1)))
    expect_error(gen.seq.HKY.binary(rtree(5), c(runif, 1, 1), runif))

    #results is a vector of length 5 (characters)
    expect_equal(length(gen.seq.HKY.binary(rtree(5), c(runif, 2, 2), c(runif, 1, 1))), 5)
    expect_is(gen.seq.HKY.binary(rtree(5), c(runif, 2, 2), c(runif, 1, 1)), "character")
    set.seed(1) ; expect_equal(unique(as.vector(gen.seq.HKY.binary(rtree(5), c(runif, 2, 2), c(runif, 1, 1)))), c("1", "0"))
})

#Testing k.sampler
test_that("k.sampler works", {
    #binary states (most of the cases)
    expect_equal(k.sampler("a"), 2)
    expect_equal(k.sampler(1), 2)
    expect_equal(k.sampler(0.5), 2)

    #multistates (up to 4 states)
    set.seed(1) ; expect_equal( sort(unique(replicate(100, k.sampler(c(0.34, 0.33, 0.33))))), c(2,3,4) )
    #Proportion respected
    set.seed(1) ; test <- replicate(10000, k.sampler(c(0.80, 0.15, 0.05)))
    expect_equal( sort(unique(test)), c(2,3,4) )
    expect_equal( length(which(test == 2))/10000, 0.7932 )
    expect_equal( length(which(test == 3))/10000, 0.1535 )
    expect_equal( length(which(test == 4))/10000, 0.0533 )
})


#Testing rTraitDisc.mk
test_that("rTraitDisc.mk works", {
    #

})

rTraitDisc.mk <- function(tree, substitution, rates, states, ...)

k.sampler <- function(states) {
    if(length(states) == 1) {
        #Only binary characters
        return(2)
    } else {
        return(sample(2:(length(states)+1), 1, prob = states))
    }
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