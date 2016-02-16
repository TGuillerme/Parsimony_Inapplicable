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
    expect_error(gen.seq.HKY.binary("a",5))
    expect_error(gen.seq.HKY.binary(5,5))

    #results is a 5 by 5 matrix
    expect_equal(dim(gen.seq.HKY.binary(rtree(5),5)), c(5,5))
    set.seed(1) ; expect_equal(unique(as.vector(gen.seq.HKY.binary(rtree(5),5))), c("1", "0"))
})



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