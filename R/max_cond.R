lasso_max_d <- function(X, Y, a, b, s)
{
    Xs <- scale(X)
    Ys <- scale(Y)
    q <- (a %*% Xs) * (b %*% Ys)

    .lasso_max_d(q, s)
}

#' Maximize the d vector in biclustering
#'
#' Find the solution to the quadratic program:
#' f(d) = d^2 * q
#' s.t. |d| <= s, 0 <= d_i <= 1
#'
#' @param q a vector of length k, where k is the number of conditions
#' @param s a regularization parameter such that s > 0 and t(d) *%* q %* d  <= s
#' @return a valid maximum
#' @export
.lasso_max_d <- function(q, s)
{
    # XXX: Not checking if s > 0... a program above should check

    d <- rep.int(0, length(q))

    n_gt0 <- length(which(q > 0))

    # no terms are positive, thus the only way to maximize is to choose all
    # values to be zero
    if (n_gt0 == 0 )
        return(d)

    q_order <- order(q, decreasing = T)
    cur_sum <- 0
    i <- 1
    for (i in 1:n_gt0)
    {
        if (cur_sum + 1 <= s)
        {
            cur_sum <- cur_sum + 1
            d[q_order[i]] <- 1
        }
        else {
            d[q_order[i]] <- s - cur_sum
            break
        }
    }

    d
}

old_lasso_max_d <- function(X, Y, a, b, s)
{
    # XXX: Possibly remove... Unsure if I should scale in this fashion
    xDf <- scale(X)
    yDf <- scale(Y)

    # Now, maximize D using quadratic programming
    # q <- - 2 * as.numeric((a %*% xDf) * (b %*% yDf))
    q <- -as.numeric((a %*% xDf) * (b %*% yDf))

    if (verbose)
        cat("Solving optimization of d\n")

    # lower bound is zero
    constrMat<- diag(length(q))
    # upper bound is one
    constrMat <- rbind(constrMat, -diag(length(q)))
    # regularize (sum(d) <= lam)
    constrMat <- rbind(constrMat, rep.int(-1, length(q)))
    constrLimit <- c(rep(0, length(q)), 
                     rep(-1, length(q)),
                     -lam)

    # FIXME: temporary until I figure out what is wrong with this function...
    # regularize (sum(d) >= lam.lwr)
    # constrMat <- rbind(constrMat, rep.int(1, length(q)))
    # constrLimit <- c(rep(0, length(q)), 
    #                  rep(-1, length(q)),
    #                  -lam, lam.lwr)

    objectiveFn <- function(d, qVec)
    {
        # sum(0.5 * d^2 * qVec)
        sum(d^2 * qVec)
    }

    optimIt <- 0
    repeat {
        optimRes <- constrOptim(d.start, objectiveFn, grad = NULL,
                                control = list(maxit = as.integer(1000000000)),
                                # control = list(maxit = 10000),
                                ui = constrMat, ci = constrLimit,
                                qVec = q
                                )

        if (optimRes$convergence != 0)
        {
            warning("Error with convergence.\n", "\tError code: ", 
                    optimRes$convergence, "\tMessage: ", optimRes$message,
                    immediate. = TRUE)
            if (optimIt >= optim.max)
            {
                warning("*****Reached maximum number of iterations of constrOptim. Exiting",
                        immediate. = TRUE)
                break
            }
            optimIt <- optimIt + 1
            d.start <- optimRes$par
            cat("\t\tRestarting optimization at new point.\n")
            cat("\t\tconstrOptim iteration: ", optimIt + 1, "\n")
        }
        else
        {
            break
        }

    }

    
    return(optimRes$par)
}
