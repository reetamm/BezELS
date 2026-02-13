#' Generate data with unit exponential margins from 4 different bivariate copulas
#'
#' @param n Number of observations
#' @param theta Dependence parameter (between 0 and 1)
#' @param copula \code{l} for logistic, \code{il} for inverted-logist, \code{al} for asymmetric logistic, \code{g} for Gaussian
#' @param tau Quantile threshold using marginal thresholding
#'
#' @return x n x 2 matrix of data with unit exponential marginal distributions
#' @return r Vector of radii
#' @return w Vector of angles
#' @return r_0_marg truncation point for truncated Gamma distribution for all n data points
#' @return above_thresh_marg TRUE/FALSE vector of whether points are above the marginal threshold
#' @return data_marg_r_0 List with n, r, w, and r_0 for data points where above_thresh_marg=T
#' @export
#'
#' @examples simdata <- gen_data_exp(n = 600, theta = 0.8, tau=0.75, copula = 'g')
#' summary(simdata$x)
gen_data_exp <- function(n,theta,copula=NULL, tau = 0.75){
    if(copula=='g'){
        x <- mvtnorm::rmvnorm(n, sigma=matrix(c(1, theta, theta, 1), 2, 2))  # generate bivariate Gaussian data
        x <- qexp(pnorm(x))
    }
    if(copula=='l'){
        x   <- evd::rbvevd(n,dep=theta, mar1=c(0, 1, 0))    # generate bivariate logistic data in Gumbel margins
        x   <- qexp(evd::pgumbel(x))                        # convert to exponential margins
    }
    if(copula=='il'){
        x <- 1 / evd::rbvevd(n, dep=theta, mar1=c(1,1,1))  # generate inverted bivariate logistic data in exponential margins
    }
    if(copula=='al'){
        x   <- evd::rbvevd(n,dep=theta, asy = c(0.6,0.8), model = 'alog',mar1=c(0, 1, 0))    # generate bivariate logistic data in Gumbel margins
        x   <- qexp(evd::pgumbel(x))                        # convert to exponential margins
    }
    r <- x[ ,1] + x[ ,2]
    w <- x[ ,1] / r

    rm <- get_r0_marg(x, tau)
    r_0_marg <- rm$r_0_marg
    above_thresh_marg <- rm$above_thresh_marg

    N_marg <- sum(above_thresh_marg)

    ## Make the data objects to use for fitting
    data_marg_r_0 <- list(N = N_marg,
                          r = r[above_thresh_marg],
                          w = w[above_thresh_marg],
                          r_0 = r_0_marg[above_thresh_marg])

    return(list(x=x, r=r, w=w, r_0_marg = r_0_marg, above_thresh_marg=above_thresh_marg,data_marg_r_0=data_marg_r_0))
}
