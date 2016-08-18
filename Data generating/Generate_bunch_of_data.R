matrices <- list()

fun.sample <- function(){sample(c(0,1,2,3,4,5,"-"), 16, replace = TRUE)}


for (i in 1:50) {
smp <- fun.sample()

while(length(grep("-", smp)) < 3) {
    smp <- fun.sample()
}

matrices[[i]] <- smp
}


trees <- rmtree(3, 16, br = 999)
trees[[4]] <- stree(16, "balanced")
trees[[5]] <- stree(16, "left")