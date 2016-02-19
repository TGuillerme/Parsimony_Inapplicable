#TESTING cust.series

context("check.matrix")


#get.contrast.matrix
test_that("get.contrast.matrix works", {
    #Errors
    expect_error(
        get.contrast.matrix(mean)
        )
    # A simple 2 by 2 matrix (0 1)
    expect_equal(
        dim(get.contrast.matrix(matrix(data = c(0,1,0,1), ncol = 2 ))), c(2,2)
        )
    expect_true(
        all(get.contrast.matrix(matrix(data = c(0,1,0,1), ncol = 2 )) == matrix(data = c(1,0,0,1), ncol = 2 ))
        )

    # A 2 by 2 with ?
    expect_equal(
        dim(get.contrast.matrix(matrix(data = c("A","B","A","?"), ncol = 2 ))), c(3,2)
        )
    expect_true(
        all(get.contrast.matrix(matrix(data = c("A","B","A","?"), ncol = 2 )) == matrix(data = c(1,0,0,1,1,1), ncol = 2 , byrow=TRUE))
        )

    # A "complex" one with inapplicables
    expect_equal(
        dim(get.contrast.matrix(matrix(data = c("A","0","-","?", "!", "A"), ncol =3 ))), c(5,4)
        )
    expect_true(
        all(get.contrast.matrix(matrix(data = c("A","0","-","?", "!", "A"), ncol =3 )) == matrix(data = c(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1, 1,1,1,1), ncol = 4 , byrow=TRUE))
        )
})

test_that("check.matrix works", {

    set.seed(1)
    #Get a random tree
    random.tree <- rcoal(10)
    #Get a random matrix
    random.matrix <- make.matrix(random.tree, characters = 50, model = "HKY", rates = c(rgamma, 1, 1), substitution = c(runif, 2, 2))
    #Erros
    expect_error(
        check.matrix("a", parsimony = "fitch", first.tree = c(phangorn::dist.hamming, phangorn::NJ), orig.tree, distance = phangorn::RF.dist)
        )
    expect_error(
        check.matrix(matrix, parsimony = "a", first.tree = c(phangorn::dist.hamming, phangorn::NJ), orig.tree, distance = phangorn::RF.dist)
        )
    expect_error(
        check.matrix(matrix, parsimony = "fitch", first.tree = "a", orig.tree, distance = phangorn::RF.dist)
        )
    expect_error(
        check.matrix(matrix, parsimony = "fitch", first.tree = c(phangorn::dist.hamming, phangorn::NJ), "a", distance = phangorn::RF.dist)
        )
    expect_error(
        check.matrix(matrix, parsimony = "fitch", first.tree = c(phangorn::dist.hamming, phangorn::NJ), orig.tree, distance = "a")
        )

    #Output
    set.seed(1)
    test <- check.matrix(random.matrix, parsimony = "sankoff", orig.tree = random.tree) # Noisy!
    expect_equal(
        dim(test), c(4,1)
        )
    expect_equal(
       test[1,], 73
        )
    expect_equal(
       test[2,], 0
        )
    expect_equal(
       round(test[3,], digit = 4), round(0.5100671, digit = 4)
        )
    expect_equal(
       test[4,], 4
        )
})
