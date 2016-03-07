#' @title Export a matrix to tnt.
#'
#' @description Generates a tnt file that contains the matrix and the tnt instructions for running a tree.
#'
#' @param matrix A discrete morphological matrix.
#' @param chain.name A string of characters to be the tnt files chain name.
#' @param inapplicable \code{logical} whether to treat inapplicable as \code{-} (\code{TRUE}) or as \code{?} (\code{FALSE}).
#' @param run.n.read \code{logical} whether to run the matrix and read in the results (\code{default = FALSE}). For UNIX machines only: see details.
#'
#' @return
#' Generates a tnt file out of the R environment. If \code{run.n.read = TRUE} also reads in the resulting trees and parsimony score.
#' 
#' @details
#' The \code{run.n.read} arguments works by calling \code{system("tnt.command < chain.name.tnt")}. Make sure the \code{tnt.command} file is in the the \code{usr/bin/} or \code{usr/local/bin/} directory.
#' 
#' @examples
#' ## A random tree with 10 tips
#' tree <- rcoal(10)
#' ## setting up the parameters
#' my_rates = c(rgamma, 1, 1) # A gamma rate distribution with of shape alpha = 0.5
#' my_substitutions = c(runif, 2, 2) # A fixed substitution rate of 2 (T/T ratio in HKY)
#'
#' ## Mk matrix (10*50)
#' matrixMk <- make.matrix(tree, characters = 50, model = "ER", rates = my_rates)
#' 
#' ## Exporting to tnt
#' export.tnt(matrixMk, chain.name = "dummy")
#' 
#' ## Tidying
#' remove.files("dummy.tnt")
#' @author Thomas Guillerme
#' @export

#source("sanitizing.R")

export.tnt <- function(matrix, chain.name, inapplicable = TRUE, run.n.read = FALSE) { #...

    #SANITIZNG
    #tree
    check.class(matrix, "matrix")

    #characters
    check.class(chain.name, "character")
    check.length(chain.name, 1, " must be a character string.")

    #states
    check.class(inapplicable, "logical")

    #Creating the matrix to print out
    matrix_out <- matrix(data = apply(matrix, 1, function(X) paste(X, collapse = "")), ncol = 1, dimnames = list(c(rownames(matrix))))
    #Adding the matrix number of character and taxa
    colnames(matrix_out) <- paste(c(ncol(matrix), nrow(matrix)), collapse = " ")

    #Exporting the table
    #write.table(matrix_out, file = paste(chain.name, ".tmp", sep = ""), quote = FALSE)

    #Exporting the text
    sink(paste(chain.name, ".tnt", sep = ""))
    #Header
    cat("xread\n\n")
    #Matrix
    print(noquote(matrix_out))
    #Log file
    cat(paste("\nlog ", paste(chain.name, ".log", sep = ""), " ;\n", sep = ""))
    #Gaps
    if(inapplicable == TRUE) {
        cat("nstates GAPS ;\n")
    } else {
        cat("nstates NOGAPS ;\n")
    }
    #Tree search and scores
    cat("\nXmult ;\nscores ;\n")
    #Export trees
    cat("export - ", paste(chain.name, ".tre", sep = "")," ;\n", sep = "")
    #Strict consensus
    cat("\nnelsen ;\n")
    #Export consensus
    cat("export - ", paste(chain.name, ".con.tre", sep = "")," ;\n", sep = "")
    #End procedure
    cat("\nprocedure/;\n")

    #Quit
    if(run.n.read == TRUE) {
        cat("quit;\n")
    }

    #Stop sinking
    sink()

    #Runing and reading
    if(run.n.read == TRUE) {
        # Running the tnt file
        system(paste("tnt.command < ", paste(chain.name, ".tnt", sep = ""), sep = ""))
    }

    
}