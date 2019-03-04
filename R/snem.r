#' @title EM algorithm for multivariate skew normal distribution.
#'
#' @description The maximum likelihood estimate for the multivariate skew-normal distribution is obtained via an EM algorithm.
#'
#' @param x A data matrix. Each row is an observation vector.
#' @param eps Weight parameter with \eqn{0 \le eps < 1}. Default is 0.9.
#' @param maxit Maximum number of iterations.
#' @param iter.eps Convergence threshold. Default is 10^-6.
#' @param stop.rule \code{"parameter"}: The difference of the parameter is used as a stopping rule. \code{"log-likelihood"} The ratio of the log-likelihood is used as a stopping rule.
#'
#' @return Location parameter (\code{mu}), covariance matrix (\code{omega}), skewness parameter (\code{delta}), and another expression of skewness parameter (\code{lambda}).
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @importFrom graphics plot
#'
#' @references Abe, T., Fujisawa, H., and Kawashima, T. (2019) \emph{EM algorithm using overparametrization for multivariate skew-normal distribution,} \emph{in preparation.}
#'
#' @details The parameter \code{eps} is a tuning parameter which ensures that an initial covariance matrix is positive semi-definite.
#'
#' @examples
#' library(DAAG)
#' x <- ais[c("bmi")]
#' snem(x, stop.rule ="log-likelihood")
#' @export
snem <- function(x, eps = 0.9, maxit = 6000, iter.eps = 10^-6, stop.rule = c("parameter","log-likelihood") ){

  if(eps < 0 || eps > 1 ){
    stop(message="[Warning] 0 < eps < 1.")
  }

  if(iter.eps < 0 ){
    stop(message="[Warning] iter.eps must be a positive real value.")
  }

  stop <- match.arg(stop.rule)

  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]

  x.colmean <- colMeans(x)
  x.center <- x-t(matrix(x.colmean, p, n))

  gam <- colMeans(x.center^3)
  gam <- sign(gam)*(abs(gam)/((4/pi-1)*sqrt(2/pi)))^(1/3)

  mu0 <- x.colmean - sqrt(2/pi)*gam

  omega0 <- (t(x.center)%*%x.center)/n + (2/pi)*gam%*%t(gam)

  tmp <- eigen(omega0)
  omega0.inv.half <- tmp$vectors %*% diag( (tmp$values)^(-1/2), p ) %*% t(tmp$vectors)
  omega0.half <- tmp$vectors %*% diag( (tmp$values)^(1/2), p ) %*% t(tmp$vectors)
  delta0 <- omega0.inv.half %*% gam

  if (sum(delta0^2) > 1) {
    delta0 <- eps*delta0 / sqrt(sum(delta0^2))
  }

  tau0 <- 1
  mu1 <- mu0
  omega1 <- omega0
  delta1 <- delta0

  if(stop=="parameter"){

    ll.val <- snll(x-t(matrix(mu0, p, n)), n, omega0.inv.half, delta0)

    for (i in 1:maxit){

      mu1 <- x.colmean - ( omega0.half  %*% delta0 )/ tau0 * mean(yexabs( x-t(matrix(mu0, p, n)), omega0.inv.half, delta0, tau0 ))

      tau1 <-  mean(yex2( x-t(matrix(mu0, p, n)), omega0.inv.half, delta0, tau0 ))^(1/2)

      omega1 <- (t(x-t(matrix(mu1, p, n)))%*%( x-t(matrix(mu1, p, n))))/n

      tmp <- eigen(omega1)
      omega1.inv.half <- tmp$vectors %*% diag( (tmp$values)^(-1/2), p) %*% t(tmp$vectors)
      omega1.half <- tmp$vectors %*% diag( (tmp$values)^(1/2), p) %*% t(tmp$vectors)

      werr <- delta1

      delta1 <- omega1.inv.half%*%t((t(yexabs( x - t(matrix(mu0, p, n)), omega0.inv.half, delta0, tau0 ))  %*% (x - t(matrix(mu1, p, n))))/n)/tau1

      werr <- werr - delta1
      werr <- sqrt(sum(werr^2))

      ll.val <- append(ll.val, snll(x-t(matrix(mu1, p, n)), n, omega1.inv.half, delta1))

      if (werr <= iter.eps){
        break
      }

      mu0 <- mu1
      omega0.inv.half <- omega1.inv.half
      omega0.half <- omega1.half
      delta0 <- delta1
      tau0 <- tau1

    }

    lambda1 <- delta1 / sqrt( 1- sum( delta1^2 ) )
    cat("\n")
    cat("stopping rule: ", stop,"\n")
    cat("iteration: ", i, "\n")
    cat("log-likelihood: ", ll.val[i+1], "\n")

    cat("mu \n")
    print(mu1)

    cat("Omega \n")
    print(unname(omega1))

    cat("delta \n")
    print(delta1)

    cat("lambda \n")
    print(lambda1)

  }else if(stop=="log-likelihood"){

    ll.val <- snll(x-t(matrix(mu0, p, n)), n, omega0.inv.half, delta0)

    for(i in 1:maxit){

      mu1 <- x.colmean - ( omega0.half  %*% delta0 )/ tau0 * mean(yexabs( x-t(matrix(mu0, p, n)), omega0.inv.half, delta0, tau0 ))

      tau1 <-  mean(yex2( x-t(matrix(mu0, p, n)), omega0.inv.half, delta0, tau0 ))^(1/2)

      omega1 <- (t(x-t(matrix(mu1, p, n)))%*%( x-t(matrix(mu1, p, n))))/n

      tmp <- eigen(omega1)
      omega1.inv.half <- tmp$vectors %*% diag( (tmp$values)^(-1/2), p) %*% t(tmp$vectors)
      omega1.half <- tmp$vectors %*% diag( (tmp$values)^(1/2), p) %*% t(tmp$vectors)

      delta1 <- omega1.inv.half%*%t((t(yexabs( x - t(matrix(mu0, p, n)), omega0.inv.half, delta0, tau0 ))  %*% (x - t(matrix(mu1, p, n))))/n)/tau1

      ll.val <- append(ll.val, snll(x-t(matrix(mu1, p, n)), n, omega1.inv.half, delta1))

      werr <- abs(ll.val[length(ll.val)]/ll.val[length(ll.val)-1]-1)

      if (werr <= iter.eps){
        break
      }

      mu0 <- mu1
      omega0.inv.half <- omega1.inv.half
      omega0.half <- omega1.half
      delta0 <- delta1
      tau0 <- tau1

    }

    lambda1 <- delta1 / sqrt( 1- sum( delta1^2 ) )
    cat("\n")
    cat("stopping rule: ", stop,"\n")
    cat("iteration: ", i,"\n")
    cat("log-likelihood: ", ll.val[i+1], "\n")

    cat("mu \n")
    print(mu1)

    cat("Omega \n")
    print(unname(omega1))

    cat("delta \n")
    print(delta1)

    cat("lambda \n")
    print(lambda1)

  }

  plot(ll.val , main="EM-algorithm for SN", xlab="step", ylab="log-likelihood")

  return( invisible(list(mu = mu1, omega = omega1, delta = delta1, lambda = lambda1)) )
}

yexabs <- function(xi, omega.half.inv, delta, tau){

  lambda <- delta/sqrt(1-sum( delta^2 ))

  w1 <- 1/sqrt(1+sum( lambda^2 ))

  wa <-  xi%*%(omega.half.inv%*%lambda)

  ww <- tau*w1*( dnorm(wa)/pnorm(wa) + wa )

  return( ww )
}

yex2 <- function(xi, omega.half.inv, delta, tau){

  lambda <- delta/sqrt(1-sum(delta^2 ))

  w1 <- 1/sqrt(1+sum(lambda^2 ))

  wa <- xi%*%(omega.half.inv%*%lambda)

  ww <- tau^2*w1^2*wa*( dnorm(wa)/pnorm(wa) + wa + 1/wa )

  return(ww)
}

snll <- function( xi, dim.row.xi, omega.half.inv, delta ){

  lambda <- delta/sqrt(1-sum(delta^2))

  wx <- xi %*% omega.half.inv

  ww<- dim.row.xi*log(2) + sum( log( pnorm( wx %*% lambda ) ) ) + sum( log( dmvnorm(wx) ) ) + dim.row.xi*log(det(omega.half.inv))

  return(ww)
}
