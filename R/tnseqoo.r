## I just learned that R has not 2, but 3 OO systems.
## I knew of S3/S4 and am reasonably functional with them, but the newer 'RC'
## system looks to me to be more 'normal' vis a vis my previous experience
## in other languages; but S3 object are so wonderfully simple.

## I am sticking with the mtb example provided in the paper
## Thus I did:
## awk '{printf("%s %s\n", $1, $3)}' mtb_stuff.tsv | uniq > mtb_annotation.tsv
## This file I am leaving in inst/ as something to be used

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

PerformEM <- function(me, maxit=600, tolerance=1e-6) {
    finished_flag <- 0
    current_mean_counts <- me$initial_means
    current_regression <- me$initial_regression
    dat <- me$dat
    ## Iterate between expectation and maximization steps
    for(i in 1:maxit) {
	cur <- c(current_mean_counts, current_regression)
	current_mean_counts_G <- current_mean_counts[dat$ord]
	## X <- rep(1, maxit)
        StepEM <- function(count, mean_counts_G, regression) {
            step <- ifelse(count == 0,
                             inverse_logit(mean_counts_G + regression),
                             0)
            return(step)
        }
	a_step <- StepEM(dat$count,
                         current_mean_counts_G,
                         as.numeric(current_regression))
	new_mean_counts <- MaximizeReadCounts(me, a_step)
	new_regression <- MaximizeRegression(a_step)
	new_step <- c(new_mean_counts, new_regression)
        ## Stop iteration if the difference between the current and new estimates is less than a tolerance level
        if(all(abs(cur - new_step) < tolerance)) {
            finished_flag <- 1
            break
        }
        ## Otherwise continue iteration
        current_mean_counts <- new_mean_counts
        current_regression <- new_regression
    }
    if(!finished_flag) {
        warning("not converge\n")
    }
    me$final_mean_counts <- current_mean_counts
    me$final_regression <- current_regression
    return(me)
}

MaximizeReadCounts <- function(me, regression_step) {
    at <- tapply(me$dat$count * (1 - regression_step), me$dat$ord, sum)
    am <- tapply((1 - regression_step), me$dat$ord, sum)
    max_mean_counts <- at / am
    return(max_mean_counts)
}

MaximizeRegression <- function(me, regression_step) {
    oldw <- getOption("warn")
    options(warn = -1)
    fit <- glm(z_temp ~ 1, family=binomial)
    options(warn = oldw)
    maximized_regression <- fit$coef
    return(maximized_regression)
}

InsertionSites <- function(me) {
    df <- me$dat
    mapping <- tapply(df[df$gene_type==0, 1],
                      df[df$gene_type==0, 1],
                      length)
    length_n <- na.omit(as.numeric(mapping))
    return(length_n)
}

EstimateMeans <- function(me) {
    means <- me$final_mean_counts
    n <- me$normal_genes
    gene_types <- me$dat$gene_type
    next_means <- means[n + 1]
    mean_counts_n <- means[1:n]
    new_mean_counts <- rep(next_means, nrow(me$dat))
    in_sites <- InsertionSites(me)
    new_mean_counts[me$dat$gene_type == FALSE] <- rep(mean_counts_n, in_sites)
    return(new_mean_counts)
}

InverseLogit <- function(number, derivative=0) {
    ret <- NULL
    if (derivative == 0) {
        ret <- 1 / (1 + exp(-number))
    } else if (derivative == 1) {
        inverse_log <- number - 2 * log(1 + exp(number))
        ret <- exp(inverse_log)
    } else if (derivative == 2) {
        first_element <- number - 3 * log(1 + exp(number))
        second_element <- 2 * number - 3 * log(1 + exp(number))
        ret <- exp(first_element) - exp(second_element)
    } else {
        warning(paste0("Unknown derivative given: ", derivative))
        ret <- number
    }
    return(ret)
}

Likelihood <- function(num, counts) {
    ## p,mu: vector with the same length.
    ret <- num + (1 - num) * exp(-counts)
    return(ret)
}

DeriveMeanRegression <- function(me, nbet=1) {
    mean_counts <- EstimateMeans(me)
    regression <- me$final_regression
    dat <- me$dat
    p <- InverseLogit(regression)
    derivp1 <- exp(-mean_counts) / (Likelihood(p, mean_counts)) ^ 2
    derivp2 <- 1 / (p - 1)
    derivp <- ifelse(dat$count == 0, derivp1, derivp2) * InverseLogit(regression, derivative=2)
    derivp <- sum(as.numeric(derivp))
    derivp1_squared <- -derivp1 ^ 2  ## WTF!? why negative the thing if we are squaring it!?
    derivp2_squared <- -derivp2 ^ 2
    derivp2 <- sum(ifelse(dat$count == 0, derivp1_squared, derivp2_squared)) * (InverseLogit(regression, derivative=1)) ^ 2
    return(derivp + derivp2)
}

DeriveRegressionSquared <- function(me) {
    dat <- me$dat
    regression <- me$final_regression[1]
    mean_count_estimate <- me$final_mean_counts
    mean_counts <- EstimateMeans(me)
    p <- InverseLogit(regression)
    derivp1 <- (1 - exp(-mean_counts)) / Likelihood(p, mean_counts)
    derivp2 <- 1 / (p - 1)
    derivp <- ifelse(dat$count == 0, derivp1, derivp2) * InverseLogit(regression, derivative=2)
    derivp1_squared <- -derivp1 ^ 2
    derivp2_squared <- -derivp2 ^ 2
    derivp2 <- sum(ifelse(dat$count == 0, derivp1_squared, derivp2_squared)) * (InverseLogit(regression, derivative=1)) ^ 2
    return(derivp + derivp2)
}

DeriveMeanSquared <- function(me) {
    mean_counts <- EstimateMeans(me)
    regression <- me$final_regression
    dat <- me$dat
    p <- InverseLogit(regression)
    deriv_mean_squared_1 <- p * (1 - p) * exp(-mean_counts) / (Likelihood(p, mean_counts)) ^ 2
    deriv_mean_squared_2 <- -dat$count / mean_counts ^ 2
    means_squared <- ifelse(dat$count == 0, deriv_mean_squared_1, deriv_mean_squared_2)
    deriv2 <- rep(NA, me$normal_genes + 1)
    pos_pseudo <- which(dat$gene_type == TRUE)
    deriv2[1] <- sum(means_squared[pos_pseudo])
    deriv2[-1] <- tapply(means_squared[-pos_pseudo],
                         as.numeric(dat[dat$gene_type == FALSE, 1]),
                         sum)
    ret <- diag(as.numeric(deriv2))
    return(ret)
}

DeriveMeanSD <- function(me, nbet=1) {
    dat <- me$dat
    tem_mb <- DeriveMeanRegression(me)
    temp_bet2 <- DeriveRegressionSquared(me)
    temp_mu2 <- DeriveMeanSquared(me)
    N <- me$normal_genes
    hessian <- matrix(0, N+1+nbet, N+1+nbet)
    hessian[1:nbet, 1:nbet] <- temp_bet2
    hessian[(nbet+1):(N+1+nbet), (nbet+1):(N+1+nbet)] <- temp_mu2
    hessian[1:nbet, (nbet+1):(N+1+nbet)] <- tem_mb
    hessian[(nbet+1):(N+1+nbet), 1:nbet] <- t(tem_mb)
    fisher <- solve(-hessian)
    var <- diag(solve(-hessian))
    ## beta, diff
    ## var[1:nbet]
    var_diff <- rep(0, N)
    for(i in 1:N) {
        var_diff[i] <- fisher[i+1+nbet, i+1+nbet]
    }
    return(var_diff)
}

ScoreEM <- function(me, sig_cutoff=0.05, p_adjust="BH") {
    variances <- DeriveMeanSD(me)
    N <- me$normal_genes
    estimates <- as.numeric(mean_counts[1:N]) - 5.5
    standard_deviations <- sqrt(variances)
    pvalues <- 2.0 * stats::pnorm(-abs(estimates / standard_deviations))
    qvalues <- round(stats::p.adjust(pvalues, p_adjust), 4)
    res <- data.frame(
        est=as.numeric(estimates),
        sd=as.numeric(standard_deviations),
        pvalue=as.numeric(as.character(pvalues)),
        qvalue=as.numeric(as.character(qvalues))
    )
    res$essential <- ifelse(res$qvalue <= sig_cutoff & res$est <= 0,
                            1, 2)
    essential_idx <- res[, "qvalue"] > sig_cutoff
    res$essential[essential_idx] <- 0
    res$Gene <- as.character(unique(dat[dat$gene_type==0, 1]))
    return(res)
}

InitializeEM <- function(me, cutoff=5) {
    dat <- me$counts
    me$normal_genes <- length(unique(dat[dat$gene_type == FALSE, 1]))
    mapping <- tapply(dat[dat$gene_type == FALSE, 1],
                      dat[dat$gene_type == FALSE, 1],
                      length)
    me$normal_insertion_sites <- na.omit(as.numeric(mapping))
    me$pseudo_positions <- which(dat$gene_type == TRUE)
    me$num_pseudo_positions <- length(me$pseudo_positions)
    dat$ord <- me$normal_genes + 1
    dat$ord[dat$gene_type == FALSE] <- rep(1:me$normal_genes, me$normal_insertion_sites)
    me$initial_means <- rep(0, me$normal_genes + 1)
    me$nonzero_genes <- tapply(dat$count > 0, dat$ord, sum)
    me$sums_by_gene <- tapply(dat$count, dat$ord, sum)
    for (i in 1:(me$normal_genes + 1)) {
        score_mug <- function(mean_counts) {
            s <- me$sums_by_gene[i]
            ss <- me$nonzero_genes[i]
            ret <- s / mean_counts - (1 / (1 - exp(-mean_counts))) * ss
            return(ret)
        }
        if (me$sums_by_gene[i] <= cutoff) {
            me$initial_means[i] <- 0.01
        } else {
            me$initial_means[i] <- stats::uniroot(score_mug, interval=c(0.001, 10000000))$root
        }
    } ## End for loop
    initial_p <- (sum(dat$count==0) - exp(-me$initial_means) %*% c(me$normal_insertion_sites, me$num_pseudo_positions)) / length(dat$count)
    me$initial_regression <- log(initial_p / (1 - initial_p))
    me$dat <- dat
    return(me)
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

