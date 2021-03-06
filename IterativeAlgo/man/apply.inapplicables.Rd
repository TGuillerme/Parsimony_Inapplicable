% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apply.inapplicables.R
\name{apply.inapplicables}
\alias{apply.inapplicables}
\title{Apply inapplicable characters to a matrix.}
\usage{
apply.inapplicables(matrix, inapplicables, tree, invariant = FALSE,
  verbose = FALSE, ...)
}
\arguments{
\item{matrix}{A discrete morphological matrix.}

\item{inapplicables}{Optional, a vector of characters inapplicability source (either \code{"character"} or \code{"clade"}; see details). The length of this vector must be at maximum half the total number of characters.}

\item{tree}{If any inapplicable source is \code{"clade"}, a tree from where to select the clades.}

\item{invariant}{Whether to allow invariant sites among the characters with inapplicable data. If \code{invariant = FALSE} the algorithm will try to remove such characters (if possible).}

\item{verbose}{Whether to be verbose or not.}

\item{...}{Any additional arguments.}
}
\description{
Apply inapplicable characters to discrete morphological matrix.
}
\details{
\itemize{
\item The \code{inapplicables} argument intakes a vector of character inapplicability source rendering a number of characters inapplicable using the following sources:
     \itemize{
         \item \code{"character"} draws inapplicable characters directly from the character matrix, ignoring the phylogeny (i.e. for a random character X, an other random character Y will have inappicable characters fro each character states 0 for character X).
         \item \code{"clade"} draws inapplicable characters from the phylogeny: it will randomly apply inapplicable characters states for some characters by randomly selecting clades from the provided tree. The algorithm randomly assigns an inapplicable token for this character for all taxa in this clade or all taxa outside this clade.
     }
For example \code{inapplicables = c(rep("character", 2), rep("clade", 2))} will generate 4 characters with inapplicable data, two using previous characters and two other using random clades.

}
}
\examples{
set.seed(4)
## A random tree with 15 tips
tree <- rcoal(15)
## setting up the parameters
my_rates = c(rgamma, 1, 1) # A gamma rate distribution with of shape alpha = 0.5
my_substitutions = c(runif, 2, 2) # A fixed substitution rate of 2 (T/T ratio in HKY)

## A Mk matrix (10*50)
matrixMk <- make.matrix(tree, characters = 100, model = "ER", states = c(0.85, 0.15), rates = my_rates)

## Setting the number and source of inapplicable characters
my_inapplicables <- c(rep("character", 5), rep("clade", 5))

## Apply some inapplicable characters to the matrix
matrix <- apply.inapplicables(matrixMk, my_inapplicables, tree)
}
\author{
Thomas Guillerme
}

