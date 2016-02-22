## Time-stamp: <Mon Feb 22 11:42:37 2016 Ashton Trey Belew (abelew@gmail.com)>

#' Given a set of hits/position and gene names
#' figure out how many possible insertion sites exist and return them in a very peculiar fashion
#'
#' @param me a tnseq object containing the data frame of interest
#'           I should probably have this return me
#' @return the na.omitted set of sites
InsertionSites <- function(me) {
    df <- me$dat
    mapping <- tapply(df[as.logical(df$gene_type) == FALSE, 1],
                      df[as.logical(df$gene_type) == FALSE, 1],
                      length)
    num_sites <- na.omit(as.numeric(mapping))
    return(num_sites)
}

#' Read an annotation file or data frame and make sure all is kosher
#'
#' @param input a data frame or filename of annotation information
#' @param ... extra arguments most likely passed to read.table
#' @return a data frame of annotation information
#' @export
ReadAnnotation <- function(input, ...) {
    arglist <- list(...) ## potentially useful stuff like 'sep' or 'header'
    if (is.null(arglist$sep)) { arglist$seq <- " " }
    if (is.null(arglist$header)) { arglist$header <- TRUE }
    df <- NULL
    if (is.null(input)) {
        return(NULL)
    } else if (class(input) == 'character' & grepl(pattern=".gff", x=input, perl=TRUE)) {
        df <- hpgltools::gff2df(input)
    } else if (class(input) == 'character') {  ## Then assume a filename
        df <- read.table(input, sep=arglist$sep, header=arglist$header)
    } else {
        df <- as.data.frame(input)  ## in case a matrix was passed or similar shenanigans
    }
    ## The TnSeq implementation I started with uses TRUE/FALSE to denote pseudogene/normal
    ## That confuses me without end, so I will state 'cds', 'pseudogene', 'ncrna', 'utr'.
    ## For the moment, change the TRUE/FALSE to pseudogene/cds and figure out how I want
    ## to deal with other annotations later
    if (is.null(df$gene_type)) {
        df$gene_type <- 'cds'
    } else if (is.logical(df$gene_type)) {
        pseudo_index <- isTRUE(df$gene_type)
        df$gene_type[pseudo_index] <- 'pseudogene'
        df$gene_type[!pseudo_index] <- 'cds'
    }
    return(df)
}

#' Read a count table and make sure that it has some stuff needed to work
#'
#' @param input a data frame of counts
#' @param ... extra arguments most likely passed to read.table
#' @return a data frame of counts by gene/position
#' @export
ReadCounts <- function(input, ...) {
    arglist <- list(...) ## potentially useful stuff like 'sep' or 'header'
    if (is.null(arglist$sep)) { arglist$sep <- " " }
    if (is.null(arglist$header)) { arglist$header <- TRUE }
    df <- NULL
    if (class(input) == 'character') {  ## Then assume a filename
        df <- read.table(input, sep=arglist$sep, header=arglist$header)
    } else {
        df <- as.data.frame(input)  ## in case a matrix was passed or similar shenanigans
    }
    return(df)
}

#' Gather a genome for analyses of possible insertion sites, gene lengths, whatever
#'
#' @param fasta a fast file to read
#' @param ... extra arguments likely to FaFile but maybe something else if I think on it
#' @return a fasta object I can use for PDict()
#' @export
ReadGenome <- function(fasta, ...) {
    arglist <- list(...) ## potentially useful stuff like 'sep' or 'header'
    if (is.null(fasta)) {
        return(NULL)
    }
    if (class(fasta) == 'character') {
        raw_seq <- try(Rsamtools::FaFile(fasta, ...))
        if (class(raw_seq)[1] == 'try-error') {
            warning(paste0("There was a problem reading: ", fasta))
            return(NULL)
        } else {
            return(raw_seq)
        }
    }
    warning("One should not get here.")
    return(NULL)
}
