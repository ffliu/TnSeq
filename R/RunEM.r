## Time-stamp: <Mon Feb 22 11:44:49 2016 Ashton Trey Belew (abelew@gmail.com)>

#' Given a tnseq object containing counts by positions and some annotation information
#' add to the tnseq object the number of normal genes, the number of insertion sites,
#' the number of pseudogene positions, and some initial estimates
#'
#' @param me the original tnseq object
#' @param cutoff an arbitrary 5 which I don't yet understand
#' @return a bigger tnseq object
#' @export
InitializeEM <- function(me, cutoff=5) {
    dat <- me$counts
    me$normal_genes <- length(unique(dat[as.logical(dat$gene_type) == FALSE, 1]))
    me$normal_insertion_sites <- InsertionSites(me)
    me$pseudo_positions <- which(as.logical(dat$gene_type) == TRUE)
    me$num_pseudo_positions <- length(me$pseudo_positions)
    dat$ord <- me$normal_genes + 1
    dat$ord[as.logical(dat$gene_type) == FALSE] <- rep(1:me$normal_genes, me$normal_insertion_sites)
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
    initial_p <- (sum(dat$count == 0) - exp(-me$initial_means) %*% c(me$normal_insertion_sites, me$num_pseudo_positions)) / length(dat$count)
    me$initial_regression <- log(initial_p / (1 - initial_p))
    me$dat <- dat
    return(me)
}

#' multiply all the read counts by the current regression state and sum them,
#' then sum the regression states, and take the ratio for each count position
#'
#' @param me a tnseq object
#' @param regression_step an individual iteration of the regression
#' @return maximized mean counts
#' @export
MaximizeReadCounts <- function(me, regression_step) {
    at <- tapply(me$dat$count * (1 - regression_step), me$dat$ord, sum)
    am <- tapply((1 - regression_step), me$dat$ord, sum)
    max_mean_counts <- at / am
    return(max_mean_counts)
}

#' Use a general linear model on the regression and return the coefficients
#' I have a fairly basic understanding of glm, I am not entirely certain how
#' this counts as a maximization, but that is probably due to my ignorance.
#'
#' @param regression_step a new regression, which is in turn a call to
#'        InverseLogit of the sum of mean_counts_G and the last regression
#' @return a maximized version of the regression step
#' @export
MaximizeRegression <- function(regression_step) {
    oldw <- getOption("warn")
    options(warn=-1)
    fit <- stats::glm(regression_step ~ 1, family=binomial)
    options(warn=oldw)
    maximized_regression <- fit$coef
    return(maximized_regression)
}

#' Perform the Expectation Maximization algorithm given a set of mean counts and tolerance
#'
#' @param me a tnseq object containing an initial set of mean counts by position
#' @param maxit the maximum number of iterations before giving up
#' @param tolerance the maximum difference between iterations before dropping out and returning an answer
#' @return me but with final_mean_counts and final_regression filled out
#' @export
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
                             InverseLogit(mean_counts_G + regression),
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

#' given a completed EM run, score the result, drop in some p/q values
#' and return a table of the results
#'
#' @param me the finalized EM run
#' @param sig_cutoff maximum adjusted p value
#' @param p_adjust the p adjustment to use
#' @return a dataframe of yes/no by gene with scores
#' @export
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
    res$Gene <- as.character(unique(
        dat[as.logical(dat$gene_type) == FALSE, 1]
        ))
    return(res)
}
