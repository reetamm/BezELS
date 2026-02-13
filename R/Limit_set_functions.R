#' Evaluate the rate parameter of the truncated Gamma distribution
#'
#' @param N Number of data points that the truncated Gamma distribution is to be fitted to
#' @param p 7x2 matrix of points that define the Bezier spline
#' @param m Vector of slopes of length N; entries are slopes of lines connecting the origin to each point
#' @param p_l x-coordinate of each point
#'
#' @return Vector of rate parameters of the truncated Gamma distribution
#' @keywords internal
#' @export
#'
#' @examples
gx = function(N=N,p=P, m = x[,2]/x[,1], p_l = x[,1]) # p has all points, m is the slope, p_l is the x-coord of the point on the line
{
    g_b <- rep(1,N);
    t_3 <- matrix(NA,length(m),3)

    p0 <- p[1,]
    p1 <- p[2,]
    p2 <- p[3,]
    t_3[,1] <- bezier_intersection(p0,p1,p2,m)


    eval = which(t_3[,1]>0)
    g_b[eval] <- bezier_eval(p0,p1,p2,t_3[eval,1])


    p0 <- p[3,]
    p1 <- p[4,]
    p2 <- p[5,]
    t_3[,2] <- bezier_intersection(p0,p1,p2,m)


    eval = which(t_3[,2]>0)
    g_b[eval] <- bezier_eval(p0,p1,p2,t_3[eval,2])


    p0 <- p[5,]
    p1 <- p[6,]
    p2 <- p[7,]
    t_3[,3] <- bezier_intersection(p0,p1,p2,m)

    eval = which(t_3[,3]>0)
    g_b[eval] <- bezier_eval(p0,p1,p2,t_3[eval,3])

    g_ratio <- p_l/g_b
    return(g_ratio)
}

#' Make a 7x2 matrix that defines a Bezier spline, based on parameters and constraints
#'
#' @param theta Vector of 9 control points that define the shape of the spline
#' @param pmix0 Unit mass probability of p0y and p6x
#' @param pmix1 Unit mass probability for p1x and p5y
#' @param pmix2 Unit mass probability for p2x and p4y
#' @param pmix3 Unit mass probability for p3
#'
#' @return P 7x2 matrix that defines a Bezier spline comprised of 3 Bezier curves
#' @keywords internal
#' @export
#'
#' @examples
makeP = function(theta,pmix0,pmix1,pmix2,pmix3){
    theta <- pnorm(theta)

    # asymmetric logistic
    theta[1] <- ifelse(theta[1]>pmix0,1,theta[1]/pmix0)
    theta[9] <- ifelse(theta[9]>pmix0,1,theta[9]/pmix0)
    # inverted logistic
    theta[2] <- ifelse(theta[2]<pmix1,0,(theta[2]-pmix1)/(1-pmix1))
    theta[8] <- ifelse(theta[8]<pmix1,0,(theta[8]-pmix1)/(1-pmix1))
    # 2 point masses for p2 and p4
    tmp <- (c(theta[4],theta[6]) - pmix2[1])/(pmix2[2] - pmix2[1])
    if(theta[4] < pmix2[1])
        tmp[1] <- 0
    if(theta[4] > pmix2[2])
        tmp[1] <- 1

    if(theta[6] < pmix2[1])
        tmp[2] <- 0
    if(theta[6] > pmix2[2])
        tmp[2] <- 1

    theta[4] <- tmp[1]
    theta[6] <- tmp[2]
    # point mass along the diagonal
    if(max(theta[4],theta[6])==1)
        theta[5] <- ifelse(theta[5]>pmix3,1,theta[5]/pmix3)


    P <- matrix(c(
        0,       theta[1],
        theta[2],theta[3],
        theta[4],       1,
        theta[5],theta[5],
        1,       theta[6],
        theta[7],theta[8],
        theta[9],      0),
        nrow = 7,byrow = T)
    return(P)
}

#' Check whether the parameter vector satisfies the constraints on the limit set
#'
#' @param theta Vector of (9) control points and (unused) truncated Gamma shape parameter
#' @param pmix0 Unit mass probability of p0y and p6x
#' @param pmix1 Unit mass probability for p1x and p5y
#' @param pmix2 Unit mass probability for p2x and p4y
#' @param pmix3 Unit mass probability for p3
#'
#' @return TRUE/FALSE if constraints are satisfied/not-satisfied
#' @keywords internal
#' @export
#'
#' @examples
thetaOK <- function(theta,pmix0,pmix1,pmix2,pmix3){
    theta <- pnorm(theta)
    # asymmetric logistic
    theta[1] <- ifelse(theta[1]>pmix0,1,theta[1]/pmix0)
    theta[9] <- ifelse(theta[9]>pmix0,1,theta[9]/pmix0)
    # inverted logistic
    theta[2] <- ifelse(theta[2]<pmix1,0,(theta[2]-pmix1)/(1-pmix1))
    theta[8] <- ifelse(theta[8]<pmix1,0,(theta[8]-pmix1)/(1-pmix1))
    # 2 point masses for p2 and p4
    tmp <- (c(theta[4],theta[6]) - pmix2[1])/(pmix2[2] - pmix2[1])
    if(theta[4] < pmix2[1])
        tmp[1] <- 0
    if(theta[4] > pmix2[2])
        tmp[1] <- 1

    if(theta[6] < pmix2[1])
        tmp[2] <- 0
    if(theta[6] > pmix2[2])
        tmp[2] <- 1

    theta[4] <- tmp[1]
    theta[6] <- tmp[2]
    # point mass along the diagonal
    if(max(theta[4],theta[6])==1)
        theta[5] <- ifelse(theta[5]>pmix3,1,theta[5]/pmix3)

    check1 <- theta[2] <= theta[4]
    check2 <- theta[8] <= theta[6]
    check3 <- theta[2]/theta[3] <= theta[4]
    check4 <- theta[6] >= theta[8]/theta[7]
    check5 <- theta[5] >= max(theta[4],theta[6])
    check6 <- 1
    # if(max(theta[2],theta[8]) == 0)
    #     check6 <- theta[5] < 1
    # print(paste(check1,check2,check3,check4,check5,check6))
    return(check1*check2*check3*check4*check5*check6)
}

#' Evaluate the truncated Gamma likelihood for a vector of radii
#'
#' @param r Vector of radii
#' @param r_0 Vector of truncation points for each radii
#' @param shape Common shape parameter for the truncated Gamma distribution
#' @param rates Vector of rate parameters for each truncated Gamma distribution
#'
#' @return Sum of the truncated Gamma log-likelihoods
#' @keywords internal
#' @export
#'
#' @examples log_lik(r = 10, r_0 = 3, shape = 2, rates = 1)
log_lik = function(r,r_0,shape,rates){
    loglik <- dgamma(r,exp(shape),rates,log = T) - log(1-pgamma(r_0,exp(shape),rates))
    return(sum(loglik))
}

#' Truncated Gamma calculations
#'
#' @param x bivariate data with exponential margins, in a matrix with 2 colums
#' @param tau quantile level for marginal threshold
#'
#' @return r_0_marg truncation point for truncated Gamma distribution
#' @return above_thresh_marg TRUE/FALSE vector of whether points are above the marginal threshold
#' @export
#'
#' @examples
get_r0_marg <- function(x, tau) {
    thresh_marg_1 <- quantile(x[ ,1], tau)
    thresh_marg_2 <- quantile(x[ ,2], tau)
    thresh_marg <- cbind(thresh_marg_1, thresh_marg_2)
    above_thresh_marg <- (x[ ,1] > thresh_marg[ ,1]) | (x[ ,2] > thresh_marg[ ,2])
    x_max_coord <- apply(x, 1, which.max)
    x_max <- apply(x, 1, max)
    x_0_marg <- x / x_max * thresh_marg[x_max_coord]
    r_0_marg <- x_0_marg[ ,1] + x_0_marg[ ,2]
    return(list(r_0_marg = r_0_marg, above_thresh_marg = above_thresh_marg))
}
