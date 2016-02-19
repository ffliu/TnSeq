#'  Calculate the inverse logit
#'
#' yeah I will admit, I don't have a clue what this does nor why
#' But the code calls e^(-1 * z), adds 1 to the result, then flips it to 1/result
#'
#' @param z  An input number to be inversed
#' @return inv_logit
invlogit <- function(z) {
    sub_exp <- exp(-1 * z)
    inv_logit <- 1 / (1 + sub_exp)
    return(inv_logit)
}

#' Calculate the inverse logit derivative
#'
#' Same as above, I don't see why this is happening
#' But it takes a number and subtracts from it 2*log(1+exp(z)) for a result
#' Then takes that and does e^result
#'
#' @param z An input number to be inverse logitted
#' @return an inverse derivative, I guess
invlogit_deriv <- function(z) {
    derivative <- z - 2 * log(1 + exp(z))
    inverse_derivative <- exp(derivative)
    return(inverse_derivative)
}

#' Take the second derivative of an inverse logit
#'
#' hmm yeah whatever
#'
#' @param z a number to be second derivatived logitted
#' @return an inverse second derivative, maybe
invlogit_second_deriv <- function(z) {
    first_derivative <- z - 3 * log(1 + exp(z))
    second_derivative <- 2 * z - 3 * log(1 + exp(z))
    first_minus_second <- exp(first_derivative) - exp(second_derivative)
    return(first_minus_second)
}

#' wtf is lik0 supposed to mean
#'
#' Given 2 vectors of the same length, add the first to
#' ((1-itself) * e^(-1 * second vector))
#' I am sure this makes sense
#'
#' @param p a vector
#' @param mu another vector
#' @return lik, no lik0 man
lik0 <- function(p, mu) {
    ## p, mu: vector with the same length.
    lik <- p + ((1-p) * exp(-mu))
    return(lik)
}

#' perform a maximization!
#'
#' this does something I am sure which includes unicorns and ponies
#'
#' @param count counts are fun
#' @param z_temp I am pretty sure you aren't supposed to name your variables z_temp
#' @return maximized pain
mstep_mu <- function(count, z_temp) {
    at <- tapply(count * (1 - z_temp), dat$ord, sum)
    am <- tapply((1 - z_temp), dat$ord, sum)
    mu_temp <- at / am
    return(mu_temp)
}

#' Expectation Step in the EM algorithm
#'
#' do some stuff, how the hell did z_temp get here?
#'
#' @param count counts!
#' @param muG oh come on really
#' @param beta I am getting no love
#' @param X oh a capital x, I am sure that makes sense
#' @param covariance is this a cov model?
#' @return z_temp, because I said so
estep <- function(count, muG, beta, X, covariance=TRUE) {
    z_temp <- 0
    if (isTRUE(covariance)) {
        z_temp <- ifelse(count == 0, invlogit(muG + as.numeric(X %*% beta)), 0)
    } else {
        z_temp <- ifelse(count == 0, invlogit(muG + beta), 0)
    }
    return(z_temp)
}

#' Maximization Step---beta in the EM algorithm
#'
#' oh so z_temp is a parameter for this function but not the previous!?
#'
#' @param X oh noes the dredded capital x
#' @param z_temp sigh
#' @param covariance is this a cov model?
#' @return the coefficients from the general linear model fit of this crap
mstep_beta <- function(X, z_temp, covariance=TRUE) {
    oldw <- getOption("warn")
    options(warn = -1)
    fit <- NULL
    if (isTRUE(covariance)) {
        fit <- glm(z_temp ~ as.factor(dat$match) + dat$GC_percent, family=binomial)
    } else {
        fit <- glm(z_temp ~ 1, family=binomial)
    }

    options(warn = oldw)
    return(fit$coef)
}

#' the main function of EM Algorithm
#'
#' This does the actual work, and for a moment I thought I understood what was happening, but that was woefully not true.
#'
#' @param mu_inits mu inits are awesome!
#' @param beta_inits I want to init some betas
#' @param maxit 500 iterations sound good to me
#' @param tolerance 1 in a million is awesome
#' @param covariance is this a covariance model?
#' @return stuff
#' @export
em_algorithm <- function(mu_inits, beta_inits, maxit=500, tolerance=1e-6, covariance=FALSE) {
    ## Initial parameter estimates
    flag <- 0
    mu_cur <- mu_inits
    beta_cur <- beta_inits
    ## Iterate between expectation and maximization steps
    for (i in 1:maxit) {
	cur <- c(mu_cur, beta_cur)
	mu_cur_G <- mu_cur[dat$ord]
	X <- rep(1,600)
	z <- estep(y, mu_cur_G, as.numeric(beta_cur), X, covariance=covariance)
	mu_new <- mstep_mu(y, z)
	beta_new <- mstep_beta(X, z, covariance=covariance)
	new_step <- c(mu_new, beta_new)

        ## Stop iteration if the difference between the current and new estimates is less than a tolerance level
        if (all(abs(cur - new_step) < tol)) {
            flag <- 1
            break
        }
        ## Otherwise continue iteration
        mu_cur <- mu_new
        beta_cur <- beta_new
    }
    if (!flag) {
        warning("not converge\n")
    }
    ret <- list(mu=mu_cur, beta=beta_cur)
    return(ret)
}

#' obtain standardard deviation for the difference between the mean for
#' each of the normal genes with the common mean for pseudogenes
#' covariate is included, the dimension of coefficient regression is 1.
#' the estimated second derivative of loglikelihood w.r.t beta
#'
#' ok
#'
#' @param betaest this reads like beta test to me
#' @param muest this is the muest mu ever!
#' @param nbet nbets are yummy
#' @param covariance is this a covariance model?
#' @return stuff
deriv_beta2 <- function(betaest, muest, nbet=1, covariance=TRUE) {
    ## betaest: the maximum likelihood estimate of regression coefficient beta
    ## muest:   the maximum likelihood estimate of mean read count for each normal genes
    ## and common mean for pseudogenes.
    ## return the estimated second derivative of loglikelihood w.r.t beta
    beta <- betaest
    mu <- muest
    mu0 <- mu[N+1]
    muN <- mu[1:N]
    ## mu for each site
    mu <- rep(mu0, nrow(dat))
    mu[dat$gene_type==0] <- rep(muN, lengthN)

    p <- NULL
    if (isTRUE(covariance)) {
        X <- model.matrix(~ as.factor(dat$match))
        p <- invlogit(X %*% beta)
    } else {
        p <- invlogit(beta)
    }

    ##  derivp<- ifelse(y == 0, derivp1, derivp2) * invlogit_deriv2(beta )
    derivp <- NULL
    if (isTRUE(covariance)) {
        derivp <- ifelse(y == 0, derivp1, derivp2) * invlogit_second_deriv(X %*% beta)
    } else {
        derivp <- ifelse(y == 0, derivp1, derivp2) * invlogit_second_deriv(beta)
    }
    derivp <- sum(as.numeric(derivp))

    derivp2 < -1 / (p - 1)
    derivp21 <- -derivp1 ^ 2
    derivp22 <- -derivp2 ^ 2
    ret <- NULL
    if (isTRUE(covariance)) {
        derivp2 <- ifelse(y == 0, derivp21, derivp22) * (invlogit_deriv(X %*% beta)) ^ 2
        ret <- derivp + derivp2
    } else {
        ## derivp21<- - (1- exp(-mu))^2 / (lik0(p,mu))^2
        ## derivp22<-  -1/(p-1)^2
        derivp2 <- sum(ifelse(y == 0, derivp21, derivp22)) * (invlogit_deriv(beta)) ^ 2
        for(i in 1:nbet) {
            for (j in 1:nbet) {
                mat1[i, j] <- t(X[, i] * derivp) %*% X[, j]
                mat2[i, j] <- t(X[, i] * derivp2) %*% X[, j]
            }
        }
        ret <- mat1 + mat2
    }
    return(ret)
}

#' the estimated second derivative of loglikelihood w.r.t mu
#'
#' yep I agree
#'
#' @param betaest more betas
#' @param muest more mus!
#' @param covariance is this a covariance model?
#' @return stuff
deriv_mu2 <- function(betaest, muest, covariance=TRUE) {
    ## betaest: the maximum likelihood estimate of regression coefficient beta
    ## muest:   the maximum likelihood estimate of mean read count for each normal genes
    ## and common mean for pseudogenes.
    ## return (N+1)*(N+1) matrix,
    ## the estimated second derivative of loglikelihood w.r.t mu.
    beta <- betaest
    mu <- muest
    mu0 <- mu[N + 1]
    muN <- mu[1:N]
    ## mu for each site
    mu <- rep(mu0,nrow(dat))
    mu[dat$gene_type==0] <- rep(muN, lengthN)

    p <- NULL
    if (isTRUE(covariance)) {
        p <- invlogit(X%*%beta)
    } else {
        p <- invlogit(beta)
    }

    deriv_mus2_1 <- p * (1 - p) * exp(-mu) / (lik0(p,mu)) ^ 2
    deriv_mus2_2 <- -y / (mu ^ 2)
    deriv_mus2 <- ifelse(y== 0, deriv_mus2_1, deriv_mus2_2)

    deriv2 <- rep(NA, N+1)
    deriv2[1] <- sum(deriv_mus2[pos_pseudo])
    deriv2[-1] <- tapply(deriv_mus2[-pos_pseudo], as.numeric(dat[dat$gene_type==0,1]), sum)

    return(diag(as.numeric(deriv2)))
}

#' the estimated second derivative of loglikelihood w.r.t beta and mu
#'
#' hmm ok
#'
#' @param betaest more betas please
#' @param muest and more mues
#' @param covariance more covariances also
#' @return stuff
deriv_betamu <- function(betaest, muest, covariance=TRUE) {
    ## betaest: the maximum likelihood estimate of regression coefficient beta
    ## muest:   the maximum likelihood estimate of mean read count for each normal genes
    ## and common mean for pseudogenes.
    ## return nbet*N matrix
    ## the estimated second derivative of loglikelihood w.r.t beta and mu
    beta <- betaest
    mu <- muest
    mu0 <- mu[N + 1]
    muN <- mu[1:N]
    ## mu for each site
    mu <- rep(mu0, nrow(dat))
    mu[dat$gene_type==0] <- rep(muN, lengthN)

    p <- NULL
    if (isTRUE(covariance)) {
        X <- model.matrix(~ as.factor(dat$match))
        p <- invlogit(X %*% beta)
    } else {
        p <- invlogit(beta)
    }

    derivp1 <- exp(-mu) / (lik0(p, mu)) ^ 2
    derivp2 <- 0

    derivp <- NULL
    if (isTRUE(covariance)) {
        derivp <- ifelse(y == 0, derivp1, derivp2) * ivlogit_deriv(X %*% beta)
    } else {
        derivp <- ifelse(y == 0, derivp1, derivp2) * invlogit_deriv(beta)
    }

    mat <- matrix(0, nbet, length(y))
    derivp <- as.numeric(derivp)
    for (i in 1:nbet) {
        if (isTRUE(covariance)) {
            mat[i, ] <- t(X[, i] * derivp)
        } else {
            mat[i, ] <- derivp
        }
  }

    res <- aggregate(t(mat), list(group=dat$ord), sum)
    return(t(res[,-1]))
}

#' derive the variance of the difference b.w. between the mean for each of
#' the normal genes with common mean for pseudogenes
#'
#' I agree!
#'
#' @param bet almost like bert
#' @param mu and nu
#' @return stuff
deriv_SD_diff<-function(bet, mu) {
    ## bet: the regression coefficients
    ## mu:   the mean read count for each normal genes
    ## and common mean for pseudogenes.
    tem_mb <- deriv_betamu(bet, mu)
    temp_bet2 <- deriv_beta2(bet, mu)
    temp_mu2 <- deriv_mu2(bet, mu)
    ## This is the hessian matrix for mu's, not for diff's
    hessian <- matrix(0, N+1+nbet, N+1+nbet)
    hessian[1:nbet, 1:nbet] <- temp_bet2
    hessian[(nbet+1):(N+1+nbet), (nbet+1):(N+1+nbet)] <- temp_mu2
    hessian[1:nbet, (nbet+1):(N+1+nbet)] <- tem_mb
    hessian[(nbet+1):(N+1+nbet), 1:nbet] <- t(tem_mb)
    fisher <- solve(-hessian)
    var <- diag(solve(-hessian))
    ##beta, diff
    ##var[1:nbet]
    var_diff <- rep(0,N)
    for (i in 1:N) {
        var_diff[i] <- fisher[i+1+nbet, i+1+nbet] + fisher[1+nbet, 1+nbet] - (2 * fisher[i+1+nbet, 1+nbet])
    }
    return(var_diff)
}

#' derive the variance for the mean offor each of
#' the normal genes
#'
#' hmm ok
#'
#' @param bet is this bert?
#' @param mu and ernie!?
#' @param covariance covar?
deriv_SD_mu <- function(bet, mu, covariance=TRUE) {
    ## bet: the regression coefficients
    ## mu:   the mean read count for each normal genes
    ## and common mean for pseudogenes.
    tem_mb <- deriv_betamu(bet, mu)
    temp_bet2 <- deriv_beta2(bet, mu)
    temp_mu2 <- deriv_mu2(bet, mu)
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
    ret <- NULL
    if (isTRUE(covariance)) {
        ret <- var
    } else {
        var_diff <- rep(0,N)
        for(i in 1:N) {
            var_diff[i] <- fisher[i+1+nbet, i+1+nbet]
        }
        ret <- var_diff
    }
    return(var_diff)
}
