#' Evaluate the analytical value of \eqn{\eta} for bivariate Gaussian, logistic, asymmetric logistic, and inverted logistic copulas
#'
#' @param theta Strength of dependence, between (0,1)
#' @param dep Copula; \code{g} for Gaussian, \code{l} for logistic, \code{il} for inverted logistic, \code{al} for asymmetric logistic
#'
#' @return scalar value of \eqn{\eta}
#' @export
#'
#' @examples
eta_analytical <- function(theta,dep=NULL){
    if(dep=='g')
        eta_a <- (1+theta)/2
    if(dep=='l')
        eta_a <- 1
    if(dep == 'il')
        eta_a <- 2^(-theta)
    if(dep=='al')
        eta_a <- 1
    return(eta_a)
}

#' Compute \eqn{\eta} based on Bezier spline control points, usually a posterior draw obtained from \code{fit_mcmc_bezier}
#'
#' @param control_points A 9x1 vector of Bezier spline control points
#'
#' @return scalar value of \eqn{\eta}
#' @export
#'
#' @examples
eta_empirical <- function(control_points){
    P <- matrix(c(
        0,       control_points[1],
        control_points[2],control_points[3],
        control_points[4],       1,
        control_points[5],control_points[5],
        1,       control_points[6],
        control_points[7],control_points[8],
        control_points[9],      0),
        nrow = 7,byrow = T)

    tt = bezier_intersection(p0 = P[3,],
                             p1 = P[4,],
                             p2 = P[5,], m = 1)
    if(tt==0 | tt == 1)
        eta=1
    if(tt>0 & tt <1)
        eta = bezier_eval(p0 = P[3,],
                          p1 = P[4,],
                          p2 = P[5,], t = tt)

    return(eta)
}

#' Evaluate the analytical value of \eqn{\lambda(\omega)} for bivariate Gaussian, logistic, asymmetric logistic, and inverted logistic copulas
#'
#' @param omega scalar between (0,1)
#' @param theta Strength of dependence, between (0,1)
#' @param dep Copula; \code{g} for Gaussian, \code{l} for logistic, \code{il} for inverted logistic, \code{al} for asymmetric logistic
#'
#'
#' @return scalar value of \eqn{\lambda(\omega)}
#' @export
#'
#' @examples
lambda_analytical <- function(theta, omega,dep=NULL){
    if(dep == 'g'){
        t = min(omega,1-omega)/max(omega,1-omega)
        if(t<theta^2)
            lambda = max(omega,1-omega)
        if(t >= theta^2)
            lambda  = (1-2*theta*sqrt(omega*(1-omega)))/(1-theta^2)
    }
    if(dep=='l'){
        lambda = max(omega,1-omega)
    }
    if(dep=='il'){
        lambda = (omega^(1/theta) + (1-omega)^(1/theta))^theta
    }
    if(dep=='al'){
        lambda = max(omega,1-omega)
    }
    return(lambda)
}

#' Compute \eqn{\lambda(\omega)} based on Bezier spline control points, usually a posterior draw obtained from \code{fit_mcmc_bezier}
#'
#' @param omega scalar between (0,1)
#' @param control_points A 9x1 vector of Bezier spline control points
#'
#'
#' @return scalar value of \eqn{\lambda(\omega)}
#' @export
#'
#' @examples
lambda_empirical <- function(control_points, omega){
    m = (1-omega)/omega
    if(m<control_points[6] | m > 1/control_points[4])
        s_omega = 1
    if(m >= control_points[6] & m <= 1/control_points[4]){
        P <- matrix(c(
            0,       control_points[1],
            control_points[2],control_points[3],
            control_points[4],       1,
            control_points[5],control_points[5],
            1,       control_points[6],
            control_points[7],control_points[8],
            control_points[9],      0),
            nrow = 7,byrow = T)
        tt = bezier_intersection(p0 = P[3,],
                                 p1 = P[4,],
                                 p2 = P[5,], m = m)
        if(omega>0.5)
            s_omega = bezier_eval(p0 = P[3,1],
                                  p1 = P[4,1],
                                  p2 = P[5,1], t = tt)
        if(omega<=0.5)
            s_omega = bezier_eval(p0 = P[3,2],
                                  p1 = P[4,2],
                                  p2 = P[5,2], t = tt)
    }
    lambda = max(omega,1-omega)/s_omega
    # return(list(s_omega = s_omega,lambda = lambda))
    return(lambda)
}

#' Evaluate the analytical value of \eqn{\tau_1(\delta)} for bivariate Gaussian, logistic, asymmetric logistic, and inverted logistic copulas
#'
#' @param delta scalar between (0,1)
#' @param theta Strength of dependence, between (0,1)
#' @param dep Copula; \code{g} for Gaussian, \code{l} for logistic, \code{il} for inverted logistic, \code{al} for asymmetric logistic
#'
#'
#' @return scalar value of \eqn{\tau_1(\delta)}
#' @export
#'
#' @examples
tau_analytical = function(theta, delta,dep=NULL){
    if(dep=='g'){
        if(delta >= theta^2)
            tau = 1
        if(delta < theta^2)
            tau = (1-theta^2)/(1+delta-2*theta*sqrt(delta))
    }
    if(dep=='l'){
        tau = theta/(1+theta*delta-delta)
    }
    if(dep=='il'){
        tau=1
    }
    if(dep=='al'){
        tau=1
    }
    return(tau)
}

#' Compute \eqn{\tau_1(\delta)} based on Bezier spline control points, usually a posterior draw obtained from \code{fit_mcmc_bezier}
#'
#' @param delta scalar between (0,1)
#' @param control_points A 9x1 vector of Bezier spline control points
#'
#' @return scalar value of \eqn{\tau_1(\delta)}
#' @export
#'
#' @examples
tau_empirical <- function(control_points, delta){
    m = delta
    if(m>=control_points[6])
        tau = 1
    if(m < control_points[6]){
        P <- matrix(c(
            0,       control_points[1],
            control_points[2],control_points[3],
            control_points[4],       1,
            control_points[5],control_points[5],
            1,       control_points[6],
            control_points[7],control_points[7],
            control_points[9],      0),
            nrow = 7,byrow = T)
        tt = bezier_intersection(p0 = P[5,],
                                 p1 = P[6,],
                                 p2 = P[7,], m = m)
        tau = bezier_eval(p0 = P[5,1],
                          p1 = P[6,1],
                          p2 = P[7,1], t = tt)
    }
    if(tau<control_points[9])
        tau <- control_points[9]
    return(tau)
}
