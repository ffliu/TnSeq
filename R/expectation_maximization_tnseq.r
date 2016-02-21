
#' calculate the number of insertion sites
#'
#' @param df a data frame of counts
#' @return a length data structure
#' @export
insertion_sites <- function(df) {
    mapping <- tapply(df[df$gene_type==0, 1],
                      df[df$gene_type==0, 1],
                      length)
    length_n <- na.omit(as.numeric(mapping))
    return(length_n)
}

#' initialize the list of means
#'
#' @param mean_estimates a data frame of counts
#' @param df a data frame of counts
#' @param in_sites the list of insertion sites
#' @return a length data structure
#' @export
initialize_means <- function(mean_estimates, df, in_sites=NULL) {
    mean_counts <- mean_estimates
    n <- length(unique(df[df$gene_type==0, 1]))
    if (is.null(in_sites)) {
        in_sites <- insertion_sites(df)
    }
    gene_types <- df$gene_type
    next_mean_counts <- mean_estimates[n + 1]
    mean_counts_n <- mean_counts[1:n]
    new_mean_counts <- rep(next_mean_counts, nrow(df))
    in_sites <- insertion_sites(df)
    new_mean_counts[df$gene_type==0] <- rep(mean_counts_n, in_sites)
    return(new_mean_counts)
}

#' Do an inverse logit
#'
#' @param number a number
#' @param derivative 0 means calculate the logit, 1 or 2 a 1st/2nd derivatives
#' @return another number
#' @export
inverse_logit <- function(number, derivative=0) {
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

#' A dirty likelihood function
#'
#' @param p a number
#' @param mean_counts a list of counts
#' @return another number
#' @export
lik0 <- function(p, mean_counts) {
    ## p,mu: vector with the same length.
    p + (1 - p) * exp(-mean_counts)
}

#' Expectation Step in the EM algorithm
#'
#' @param counts our count list
#' @param mean_counts_G dunno
#' @param regression the glm regression
#' @return z_temp which is 0 or an inverse_logit
#' @export
estep <- function(count, mean_counts_G, regression) {
    z_temp <- ifelse(count == 0,
                     inverse_logit(mean_counts_G + regression),
                     0)
    return(z_temp)
}

#' Maximization Step---mu in the EM algorithm
#'
#' @param dat the data frame
#' @param count the count list
#' @param z_temp something else
#' @return a maximum
#' @export
maximize_read_counts <- function(dat, count, z_temp) {
    at <- tapply(count * (1 - z_temp), dat$ord, sum)
    am <- tapply((1 - z_temp), dat$ord, sum)
    max_mean_counts <- at / am
    return(max_mean_counts)
}

#' Maximization Step---beta in the EM algorithm
#'
#' @param z_temp a number
#' @return a maximized regression
#' @export
maximize_regression_coef <- function(z_temp) {
    oldw <- getOption("warn")
    options(warn = -1)
    fit <- glm(z_temp ~ 1, family=binomial)
    options(warn = oldw)
    maximized_regression <- fit$coef
    return(maximized_regression)
}

#' the main function of EM Algorithm
#'
#' @param dat the data frame
#' @param initial_mean_counts this really doesn't belong here I think
#' @param initial_regression neither does this
#' @param maxit maximum number of iterations before it drops out
#' @param tolerance minimum difference between two iterations before
#'     it pops out and returns
#' @return stuff
#' @export
em.algo <- function(dat, initial_mean_counts, initial_regression,
                    maxit=600, tolerance=1e-6) {
    ## Initial parameter estimates
    flag <- 0
    current_mean_counts <- initial_mean_counts
    current_regression <- initial_regression
    ## Iterate between expectation and maximization steps
    for(i in 1:maxit) {
	cur <- c(current_mean_counts, current_regression)
	current_mean_counts_G <- current_mean_counts[dat$ord]
	## X <- rep(1, maxit)
	a_step <- estep(dat$count,
                        current_mean_counts_G,
                        as.numeric(current_regression))
	new_mean_counts <- maximize_read_counts(dat, dat$count, a_step)
	new_regression <- maximize_regression_coef(a_step)
	new_step <- c(new_mean_counts, new_regression)
        ## Stop iteration if the difference between the current and new estimates is less than a tolerance level
        if(all(abs(cur - new_step) < tolerance)) {
            flag <- 1
            break
        }
        ## Otherwise continue iteration
        current_mean_counts <- new_mean_counts
        current_regression <- new_regression
    }
    if(!flag) {
        warning("not converge\n")
    }
    list(mean_counts=current_mean_counts, regression=current_regression)
}

#######################################
## obtain standardard deviation for the difference between the mean for
## each of the normal genes with the common mean for pseudogenes
#######################################
## covariate is included, the dimension of coefficient regression is 1.
## nbet <- 1

#' the estimated second derivative of loglikelihood w.r.t beta
#'
#' @param dat the data frame
#' @param regression_estimate_coef coefficients of the regression
#' @param mean_count_estimate guess!
#' @return stuff
#' @export
deriv_beta2 <- function(dat, regression_estimate_coef, mean_count_estimate) {
    ## regression_estimate_coef: the maximum likelihood estimate of regression coefficient beta
    ## mean_count_estimate:   the maximum likelihood estimate of mean read count for each normal genes
    ## and common mean for pseudogenes.
    ## return the estimated second derivative of loglikelihood w.r.t beta
    mean_counts <- initialize_means(mean_count_estimate, dat)
    regression <- regression_estimate_coef
    ##mean_counts <- mean_count_estimate
    ##mean_counts_0 <- mean_counts[N + 1]
    ##mean_counts_N <- mean_counts[1:N]
    ## mu for each site
    ##mean_counts <- rep(mean_counts_0, nrow(dat))
    ##mean_counts[dat$gene_type==0] <- rep(mean_counts_N, lengthN)
    p <- inverse_logit(regression)
    derivp1 <- (1 - exp(-mean_counts)) / lik0(p, mean_counts)
    derivp2 <- 1 / (p - 1)
    derivp <- ifelse(dat$count==0, derivp1, derivp2) *
        inverse_logit(regression, derivative=2)
    derivp <- sum(as.numeric(derivp))
    derivp21 <- -derivp1 ^ 2
    derivp22 <- -derivp2 ^ 2
    ## derivp21 <- -(1 - exp(-mu)) ^ 2 / (lik0(p, mu)) ^ 2
    ## derivp22 <- -1 / (p - 1) ^ 2
    derivp2 <- sum(ifelse(dat$count==0, derivp21, derivp22)) * (inverse_logit(regression, derivative=1)) ^ 2
    return(derivp + derivp2)
}

#' the estimated second derivative of loglikelihood w.r.t mu
#'
#' @param dat the data frame
#' @param regression_estimate_coef the coefficients of the regression
#' @param mean_count_estimate guess!
#' @return stuff
#' @export
deriv_mu2 <- function(dat, regression_estimate_coef, mean_count_estimate) {
    ## regression_estimate_coef: the maximum likelihood estimate of regression coefficient beta
    ## mean_count_estimate:   the maximum likelihood estimate of mean read count for each normal genes
    ## and common mean for pseudogenes.
    ## return (N+1)*(N+1) matrix,
    ## the estimated second derivative of loglikelihood w.r.t mu.
    mean_counts <- initialize_means(mean_count_estimate, dat)
    regression <- regression_estimate_coef
    ##mean_counts <- mean_count_estimate
    ##mean_counts_0 <- mean_counts[N + 1]
    ##mean_counts_N <- mean_counts[1:N]
    ## mu for each site
    ##mean_counts <- rep(mean_counts_0, nrow(dat))
    ##mean_counts[dat$gene_type==0] <- rep(mean_counts_N, lengthN)
    p <- inverse_logit(regression)
    deriv_mus2_1 <- p * (1 - p) * exp(-mean_counts) / (lik0(p, mean_counts)) ^ 2
    deriv_mus2_2 <- -dat$count / mean_counts ^ 2
    deriv_mus2 <- ifelse(dat$count==0, deriv_mus2_1, deriv_mus2_2)
    ## I want to make this N into a slot in a tnseq object
    N <- length(unique(dat[dat$gene_type==0, 1]))
    deriv2 <- rep(NA, N + 1)
    pos_pseudo <- which(dat$gene_type==1)
    deriv2[1] <- sum(deriv_mus2[pos_pseudo])
    deriv2[-1] <- tapply(deriv_mus2[-pos_pseudo],
                         as.numeric(dat[dat$gene_type==0, 1]),
                         sum)
    return(diag(as.numeric(deriv2)))
}

#' the estimated second derivative of loglikelihood w.r.t beta and mu
#'
#' @param dat the data frame
#' @param regression_estimate_coef bingo!
#' @param mean_count_estimate the mean counts
#' @param nbet I suspect this has something to do with the covariance
#'     model or lack thereof
#' @return stuff
#' @export
deriv_betamu <- function(dat, regression_estimate_coef,
                         mean_count_estimate, nbet=1) {
    ## regression_estimate_coef: the maximum likelihood estimate of regression coefficient beta
    ## mean_count_estimate:   the maximum likelihood estimate of mean read count for each normal genes
    ## and common mean for pseudogenes.
    ## return nbet*N matrix
    ## the estimated second derivative of loglikelihood w.r.t beta and mu
    mean_counts <- initialize_means(mean_count_estimate, dat)
    regression <- regression_estimate_coef
    ##mean_counts <- mean_count_estimate
    ##mean_counts_0 <- mean_counts[N + 1]
    ##mean_counts_N <- mean_counts[1:N]
    ## mu for each site
    ##mean_counts <- rep(mean_counts_0, nrow(dat))
    ##mean_counts[dat$gene_type==0] <- rep(mean_counts_N, lengthN)
    p <- inverse_logit(regression)
    derivp1 <- exp(-mean_counts) / (lik0(p, mean_counts)) ^ 2
    derivp2 <- 0
    derivp <- ifelse(dat$count==0, derivp1, derivp2) *
        inverse_logit(regression, derivative=1)
    mat <- matrix(0, nbet, length(dat$count))
    derivp <- as.numeric(derivp)
    for(i in 1:nbet) {
        mat[i, ] <- derivp
    }
    res <- aggregate(t(mat), list(group=dat$ord), sum)
    return(t(res[, -1]))
}

#' derive the variance of the difference b.w. between the mean for each of
#' the normal genes with common mean for pseudogenes
#'
#' @param dat my data
#' @param regression_coef the regression coefficients
#' @param mean_counts the counts
#' @param nbet the nbetas or something
#' @return stuff
#' @export
deriv_SD_diff <- function(dat, regression_coef, mean_counts, nbet=1) {
    ## bet: the regression coefficients
    ## mu:   the mean read count for each normal genes
    ## and common mean for pseudogenes.
    tem_mb <- deriv_betamu(dat, regression_coef, mean_counts)
    temp_bet2 <- deriv_beta2(dat, regression_coef, mean_counts)
    temp_mu2 <- deriv_mu2(dat, regression_coef, mean_counts)
    N <- length(unique(dat[dat$gene_type==0, 1]))
    ## where does nbet come from?
    ## This is the hessian matrix for mu's, not for diff's
    hessian <- matrix(0, N+1+nbet, N+1+nbet)
    hessian[1:nbet, 1:nbet] <- temp_bet2
    hessian[(nbet+1):(N+1+nbet), (nbet+1):(N+1+nbet)] <- temp_mu2
    hessian[1:nbet, (nbet+1):(N+1+nbet)] <- tem_mb
    hessian[(nbet+1):(N+1+nbet), 1:nbet] <- t(tem_mb)
    fisher <- solve(-hessian)
    var <- diag(solve(-hessian))
    ## beta, diff
    ## var[1:nbet]
    var_diff <- rep(0,N)
    for(i in 1:N) {
        var_diff[i] <- fisher[i+1+nbet, i+1+nbet] + fisher[1+nbet, 1+nbet] -2*fisher[i+1+nbet, 1+nbet]
    }
return(var_diff)
}

#' derive the variance for the mean offor each of
#' the normal genes
#'
#' @param dat the data
#' @param regression_coef the regression coefficients
#' @param mean_counts the counts
#' @param nbet the nbets
#' @return stuff
#' @export
deriv_SD_mu <- function(dat, regression_coef, mean_counts, nbet=1) {
    ## Ok, there is a problem here
    ## In the original source file there is a global nbet <- 1
    ## I made it a parameter and set it for each function
    ## Unfortunately it appears that it was not getting used in every
    ## instance, as a result, the 'sd' numbers are +1 more
    ## for the iterations which use this code as a library
    ## as opposed to those which use it as a script
    ## having said that, it does not affect any results
    ## but does shift the p-values by ~ 1e-6
    ## This is another reason I need to make some of this a
    ## first-class object, so that I can change those types of
    ## parameters and see what happens to the data

    ## bet: the regression coefficients
    ## mu:   the mean read count for each normal genes
    ## and common mean for pseudogenes.
    tem_mb <- deriv_betamu(dat, regression_coef, mean_counts)
    temp_bet2 <- deriv_beta2(dat, regression_coef, mean_counts)
    temp_mu2 <- deriv_mu2(dat, regression_coef, mean_counts)
    ## This is the hessian matrix for mu's, not for diff's
    N <- length(unique(dat[dat$gene_type==0, 1]))
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

#' Given a dataframe with gene names, counts by position, and a
#' boolean describing whether they are pseudogenes, perform the
#' expectation maximization algorithm
#'
#' @param dat A data frame with (currently) 3 columns, gene, count,
#'        and gene_type (pseudogene is TRUE)
#' @param covariance when performing the linear fits, have a
#'        covariance?
#' @return a matrix (soon to be data frame) of the data
#'         probably will be made into a list with other data too
#' @export
perform_tnseq_em <- function(dat, covariance=FALSE) {
    ## number of normal genes.
    N <- length(unique(dat[dat$gene_type==0, 1]))
    ## number of possible insertion site for each of the normal genes
    lengthN <- insertion_sites(dat)
    pos_pseudo <- which(dat$gene_type==1)
    length_ps <- length(pos_pseudo)
    ## add a column "ord" to denote normal genes ID and use the sampe ID for all the pseudogenes.
    dat$ord <- N + 1
    dat$ord[dat$gene_type==0] <- rep(1:N, lengthN)
    ## Obtain initial parameter estimates for EM algorithm.
    initial_mean_counts <- rep(0, N + 1)
    ## y <- dat$count
    y_numn0_gene <- tapply(dat$count > 0, dat$ord, sum)
    y_sum_gene <- tapply(dat$count, dat$ord, sum)
    for (i in 1:(N + 1)) {
        score_mug <- function(mean_counts) {
            s <- y_sum_gene[i]
            ss <- y_numn0_gene[i]
            s / mean_counts - (1 / (1 - exp(-mean_counts))) * ss
        }
        if (y_sum_gene[i] <=5) {
            initial_mean_counts[i]==0.01
        } else {
            initial_mean_counts[i] <- uniroot(score_mug, interval=c(0.001, 10000000))$root
        }
    }
    p0 <- (sum(dat$count==0) - exp(-initial_mean_counts) %*% c(lengthN, length_ps)) / length(dat$count)
    initial_regression <- log(p0 / (1 - p0))
    ## Obtain MLE for parameters using EM Algorithm
    message("About to start the expectation maximization, this may take a while.")
    output <- em.algo(dat, initial_mean_counts, initial_regression)
    ## Obtain SD
    ## var_dif<-deriv_SD_diff(output$beta,output$mu)
    ## est_dif<- as.numeric(output$mu[1:N])-output$mu[N+1]
    ## Using median of pseudogenes (5.5), Obtain SD
    message("Finished the expectation maximization, creating results df.")
    res <- make_result_df(dat, output)
    ## var_dif <- deriv_SD_mu(dat, output$regression[1], output$mean_counts)
    ## est_dif <- as.numeric(output$mean_counts[1:N]) - 5.5
    ## res <- matrix(NA, length(var_dif), 6)
    ## res[, 1] <- est_dif
    ## res[, 2] <- sqrt(var_dif)
    ## p value
    ## res[, 3] <- 2 * pnorm(-abs(est_dif / sqrt(var_dif)))
    ## q value
    ## res[, 4] <- round(p.adjust(res[, 3], "BH"), 4)
    ## essential: 1, not significant: 0, the other side: 2
    ## res[, 5] <- ifelse(res[, 4] <= 0.05 & res[, 1] <= 0, 1, 2)
    ## res[res[, 4] > 0.05, 5] <- 0
    ## res[, 6] <- as.character(unique(dat[dat$gene_type==0, 1]))
    ##colnames(res) <- c("est","sd","pvalue","qvalue","essentialp05","Gene")
    return(res)
}

#' create a data frame of results from the output of the em algorithm
#'
#' @param dat the original data
#' @param em_output the result from em
#' @param sig_cutoff maximum qvalue to accept as 'essential'
#' @param p_adjust method to use to adjust the pvalues
#' @result a data frame of 6 columns
#' @export
make_result_df <- function(dat, em_output, sig_cutoff=0.05, p_adjust="BH") {
    regression <- em_output$regression[1]
    mean_counts <- em_output$mean_counts
    variances <- deriv_SD_mu(dat, regression, mean_counts)
    N <- length(unique(dat[dat$gene_type==0, 1]))
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
