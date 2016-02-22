## Time-stamp: <Mon Feb 22 11:46:03 2016 Ashton Trey Belew (abelew@gmail.com)>

#' extract the standard deviation of a set of means using the evil nbets
#' yet another function which makes 0 sense to me
#'
#' @param me the tnseq object
#' @param nbet rawr1
#' @return the standard deviations
#' @export
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

#' Similar to DeriveRegressionSquared I do not understand what is going on here, either
#'
#' @param me the tnseq object with the final regression state and mean counts
#' @return a new set of means squared using fancy math
#' @export
DeriveMeanSquared <- function(me) {
    mean_counts <- EstimateMeans(me)
    regression <- me$final_regression
    dat <- me$dat
    p <- InverseLogit(regression)
    deriv_mean_squared_1 <- p * (1 - p) * exp(-mean_counts) / (Likelihood(p, mean_counts)) ^ 2
    deriv_mean_squared_2 <- -dat$count / mean_counts ^ 2
    means_squared <- ifelse(dat$count == 0, deriv_mean_squared_1, deriv_mean_squared_2)
    deriv2 <- rep(NA, me$normal_genes + 1)
    pos_pseudo <- which(as.logical(dat$gene_type) == TRUE)
    deriv2[1] <- sum(means_squared[pos_pseudo])
    deriv2[-1] <- tapply(means_squared[-pos_pseudo],
                         as.numeric(dat[as.logical(dat$gene_type) == FALSE, 1]),
                         sum)
    ## pull the diagonal matrix from a matrix of logits or counts/(mean_counts^2) hmm
    ret <- diag(as.numeric(deriv2))
    return(ret)
}

#' get a mean regression from a tnseq object
#' Once an EM has been performed, find the mean of the regressions?
#'
#' @param me the tnseq object including the final regression state
#' @param nbet I suspect this is related to the difference between a covariance model and not
#'        but as it stands it looks kind of dumb to me
#' @export
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

#' ok I don't understand what is the point of -derivp1 ^ 2
#' This is one of a shrinking majority of functions in this code for which I have no understanding
#'
#' @param me the tnseq object including the final regression state
#' @return the sum of derivp and derivp2 of course, gosh
#' @export
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

#' make a list of the normal genes and put into it the estimated (means/position)/gene
#'
#' @param me the tnseq object containing counts/position, the number of normal genes, and gene types
#' @return estimated means, probably should be changed to me
EstimateMeans <- function(me) {
    means <- me$final_mean_counts
    n <- me$normal_genes
    ## gene_types <- me$dat$gene_type
    next_means <- means[n + 1]
    mean_counts_n <- means[1:n]
    new_mean_counts <- rep(next_means, nrow(me$dat))
    in_sites <- InsertionSites(me)
    new_mean_counts[me$dat$gene_type == FALSE] <- rep(mean_counts_n, in_sites)
    return(new_mean_counts)
}

#' Given a number, calculate its inverse logit, the derivative of it, or second derivative
#'
#' @param number the number to modify
#' @param derivative 0 means just calculate the value, 1 or 2 as for derivatives
#' @return a new number
#' @export
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

#' a cheapo likelihood function
#'
#' @param num given a number
#' @param counts and some count data
#' @return where does that number fit with respect to the counts using an exponent?
#' @export
Likelihood <- function(num, counts) {
    ## p,mu: vector with the same length.
    ret <- num + (1 - num) * exp(-counts)
    return(ret)
}
