## Time-stamp: <Mon Feb 22 11:33:53 2016 Ashton Trey Belew (abelew@gmail.com)>

## I just learned that R has not 2, but 3 OO systems.
## I knew of S3/S4 and am reasonably functional with them, but the newer 'RC'
## system looks to me to be more 'normal' vis a vis my previous experience
## in other languages; but S3 object are so wonderfully simple.

## I am sticking with the mtb example provided in the paper
## Thus I did:
## awk '{printf("%s %s\n", $1, $3)}' mtb_stuff.tsv | uniq > mtb_annotation.tsv
## This file I am leaving in inst/ as something to be used

## For the moment I am exporting a lot of the functions here so that I can poke at them with less typing.
## That is sort of lazy I admit, but I am lazy.

#' Create a TnSeq object given the camelcaseish nature of the original source
#' code, perhaps I will finally start following the google/etc code style guide.
#' The code from which I am starting provides a position sorted data frame with
#' gene names as a field, this is a bit peculiar and I will likely take some
#' liberties with it, especially given that other code I have written is able
#' to figure out the number of possible insertion sites denovo.
#'
#' This begs a question of Yoann, should I assume any other insertion sites
#' than TA?  All tnseq work I have read uses mariner, but that isn't necessarily
#' the only transposable element one works with?
#'
#' @param count_table either a data frame or filename with counts by position
#' @param annotation a data frame or filename with annotation information
#' @param genome a filename containing a fasta genome
#' @return a TnSeq object
#' @export
TnSeq <- function(count_table, annotation=NULL, genome=NULL, task="em_nocov", ...) {
    arglist <- list(...)
    me <- list(counts = ReadCounts(count_table, ...),
               annotation = ReadAnnotation(annotation, ...),
               genome = ReadGenome(genome, ...))
    class(me) <- append(class(me), "TnSeq")
    return(me)
}
