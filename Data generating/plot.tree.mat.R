plot.tree.mat <- function(char1, A4 = TRUE) {

    if(A4 == TRUE) {quartz(width = 8.3, height = 11.7)}

    par(mfrow = (c(3,2)), bty = "n")

    tree <- read.tree(text="((((1,2),((3,4),5)),(6,7)),((((8,9),10),(11,(12,13))),(14,(15,16))));")
    plot(tree, type = "c", label.offset = 1, main = "tree 1")
    char <- c(substring(char1, seq(1,nchar(char1),1), seq(1,nchar(char1),1)))
    tiplabels(char[as.numeric(tree$tip.label)], frame = "none", adj = -1, col = "grey")

    tree <- read.tree(text="((((12,2),(7,16)),(5,4)),(((6,(9,10)),((((14,1),(13,11)),8),15)),3));")
    plot(tree, type = "c", label.offset = 1, main = "tree 2")
    char <- c(substring(char1, seq(1,nchar(char1),1), seq(1,nchar(char1),1)))
    tiplabels(char[as.numeric(tree$tip.label)], frame = "none", adj = -1, col = "grey")

    tree <- read.tree(text="(((8,4),10),(((5,(9,((6,14),(((11,3),7),2)))),((16,12),(13,15))),1));")
    plot(tree, type = "c", label.offset = 1, main = "tree 3")
    char <- c(substring(char1, seq(1,nchar(char1),1), seq(1,nchar(char1),1)))
    tiplabels(char[as.numeric(tree$tip.label)], frame = "none", adj = -1, col = "grey")

    tree <- read.tree(text="((((7,2),(10,1)),((9,13),(6,5))),(((15,11),(12,8)),((3,16),(14,4))));")
    plot(tree, type = "c", label.offset = 1, main = "tree 4")
    char <- c(substring(char1, seq(1,nchar(char1),1), seq(1,nchar(char1),1)))
    tiplabels(char[as.numeric(tree$tip.label)], frame = "none", adj = -1, col = "grey")

    tree <- read.tree(text="(7,(2,(10,(1,(9,(13,(6,(5,(15,(11,(12,(8,(3,(16,(14,4)))))))))))))));")
    plot(tree, type = "c", label.offset = 1, main = "tree 5")
    char <- c(substring(char1, seq(1,nchar(char1),1), seq(1,nchar(char1),1)))
    tiplabels(char[as.numeric(tree$tip.label)], frame = "none", adj = -1, col = "grey")
}