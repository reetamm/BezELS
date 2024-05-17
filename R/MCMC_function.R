#' @title Main MCMC fit function
#' @description This function is the main MCMC function to fit a 2D Bezier spline to the limit set boundary of bivariate data in exponential margins
#' @param N Number of data points that the truncated Gamma distribution is to be fitted to
#' @param r Vector of radii of length N
#' @param w Vector of angles of length N
#' @param r_0 Vector of truncation thresholds of length N
#' @param priors list of priors
#' @param inits list of initial values
#' @param pmix list of mixture components
#' @param iters Total number of iterations
#' @param burn Number of burn-in iterations
#' @param print.every How many iterations between output
#' @param traceplot Plot traceplot (TRUE/FALSE)
#'
#' @return \code{samples} Matrix of MCMC samples of Bezier control points, shape parameter, and eta
#' @export
#'
#' @examples set.seed(1)
#' simdata             <- gen_data_exp(n = 500, theta = 0.3, tau=0.75, copula = 'l')
#' x                   <- simdata$x
#' r                   <- simdata$r
#' w                   <- simdata$w
#' data_marg_r_0       <- simdata$data_marg_r_0
#' samples  <- fit_mcmc_bezier( N = data_marg_r_0$N,
#'                          r = data_marg_r_0$r,
#'                          w = data_marg_r_0$w,
#'                          r_0 = data_marg_r_0$r_0,
#'                          iters = 1100, burn = 100,
#'                          traceplot=T, print.every = 100)
#' median(samples[101:1100,11]) # posterior median of eta
fit_mcmc_bezier = function(N,r,w,r_0,
                       priors = list(p0y_mean = 0, p0y_sd = 1, p1x_mean = 0, p1x_sd = 1,
                                     p1y_mean = 0, p1y_sd = 1, p2x_mean = 0, p2x_sd = 1,
                                     p3_mean = 0,  p3_sd = 1, p4y_mean = 0, p4y_sd = 1,
                                     p5x_mean = 0, p5x_sd = 1, p5y_mean = 0, p5y_sd = 1,
                                     p6x_mean = 0, p6x_sd = 1, alpha_mean = 0, alpha_sd = 1),
                       inits = list(p0y = qnorm(0.5), p1x = qnorm(0.01), p1y = qnorm(0.99),
                                    p2x = qnorm(0.5), p3  = qnorm(0.8), p4y = qnorm(0.5),
                                    p5y = qnorm(0.01), p5x = qnorm(0.99), p6x = qnorm(0.5)),
                       pmix = list(pmix0 = 0.1, pmix1 = 0.1, pmix2 = c(0.1,0.1),pmix3 = 0.4),
                       iters = 11000, burn = 1000, print.every = 1000, traceplot=T){
    ## Priors
    p0y_mean    <- priors$p0y_mean
    p0y_sd      <- priors$p0y_sd
    p1x_mean    <- priors$p1x_mean
    p1x_sd      <- priors$p1x_sd
    p1y_mean    <- priors$p1y_mean
    p1y_sd      <- priors$p1y_sd
    p2x_mean    <- priors$p2x_mean
    p2x_sd      <- priors$p2x_sd
    p3_mean     <- priors$p3_mean
    p3_sd       <- priors$p3_sd
    p4y_mean    <- priors$p4y_mean
    p4y_sd      <- priors$p4y_sd
    p5x_mean    <- priors$p5x_mean
    p5x_sd      <- priors$p5x_sd
    p5y_mean    <- priors$p5y_mean
    p5y_sd      <- priors$p5y_sd
    p6x_mean    <- priors$p6x_mean
    p6x_sd      <- priors$p6x_sd
    alpha_mean  <- priors$alpha_mean
    alpha_sd    <- priors$alpha_sd

    ## inits
    p0y <- inits$p0y
    p1x <- inits$p1x
    p1y <- inits$p1y
    p2x <- inits$p2x
    p3  <- inits$p3
    p4y <- inits$p4y
    p5y <- inits$p5y
    p5x <- inits$p5x
    p6x <- inits$p6x

    ## Mixture probabilityes
    pmix0   <- 1 - pmix$pmix0
    pmix1   <- pmix$pmix1
    pmix2   <- pmix$pmix2
    pmix2[2]<- 1 - pmix2[2]
    pmix3   <- 1 - pmix$pmix3

    # this will be the 'data' for the MCMC
    x <- matrix(NA,N,2)
    x[,1] <- w
    x[,2] <- 1-w

    theta_mn    <- c(p0y_mean,p1x_mean,p1y_mean,p2x_mean,p3_mean,p4y_mean,p5x_mean,p5y_mean,p6x_mean,alpha_mean)
    theta_sd    <- c(p0y_sd,  p1x_sd,  p1y_sd,  p2x_sd,  p3_sd,  p4y_sd,  p5x_sd,  p5y_sd,  p6x_sd,  alpha_sd)
    num_p       <- length(theta_mn)
    var_names   <- c('p0y','p1x','p1y','p2x','p3','p4y','p5x','p5y','p6x','alpha')
    can_sd      <- rep(.1,num_p)
    theta       <- c(p0y,p1x,p1y,p2x,p3,p4y,p5x,p5y,p6x,alpha_mean)
    gamma_rate  <- rep(NA,N)
    gamma_rate_can  <- rep(NA,N)
    P           <- makeP(theta,pmix0,pmix1,pmix2,pmix3)
    # for(i in 1:N){
    #     gamma_rate[i] <- gx(p = P,m = x[i,2]/x[i,1],p_l = x[i,1])
    # }
    gamma_rate <- gx(N=N,p = P,m = x[,2]/x[,1],p_l = x[,1])

    check   <- 50       # Iterations between checks of the acceptance rate
    att     <- rep(0,num_p)  # Keep track of the number of MH attempts
    acc     <- rep(0,num_p)  # Keep track of the number of MH accepts
    samples <- matrix(NA,iters,num_p+1)

    # begin iterations
    for(iter in 1:iters){
        for(j in 1:num_p){
            att[j] <- att[j] + 1
            can    <- theta
            can[j] <- rnorm(1,theta[j],can_sd[j])
            gamma_rate_can <- gamma_rate

            if(thetaOK(can,pmix0,pmix1,pmix2,pmix3)){
                if(j!=10){ #update the P matrix for every parameter except the trunc-gamma shape
                    P_can <- makeP(can,pmix0,pmix1,pmix2,pmix3)
                    gamma_rate_can <- gx(N=N,P_can,x[,2]/x[,1],x[,1])

                }
                if(j<=10){
                    cur_ll <- log_lik(r,r_0,theta[num_p],    gamma_rate)
                    can_ll <- log_lik(r,r_0,  can[num_p],gamma_rate_can)
                }
                can_prior   <- dnorm(  can[j],theta_mn[j],theta_sd[j],log = T)
                cur_prior   <- dnorm(theta[j],theta_mn[j],theta_sd[j],log = T)

                mhprob   <-  can_ll - cur_ll + can_prior - cur_prior

                if(!is.na(mhprob) & !is.nan(mhprob) & !is.infinite(mhprob))
                    if(log(runif(1))<mhprob){
                        theta       <- can
                        P           <- P_can
                        gamma_rate  <- gamma_rate_can
                        acc[j]      <- acc[j] + 1
                    }
            }
        }

        # tuning
        for(j in 1:length(att)){
            if(iter<burn & att[j]==check){
                if(acc[j]/att[j]<0.2){can_sd[j]<-can_sd[j]*0.8}
                if(acc[j]/att[j]>0.6){can_sd[j]<-can_sd[j]*1.2}
                acc[j] <- att[j] <- 0
            }
        }

        # save outputs
        samples[iter,1:10]  <- c(theta)
        samples[iter,10]    <- exp(samples[iter,10])    # shape parameter
        samples[iter,1:9]   <- pnorm(samples[iter,1:9]) # other parameters

        # asymmetric logistic
        samples[iter,1] <- ifelse(samples[iter,1]>pmix0,1,samples[iter,1]/pmix0)
        samples[iter,9] <- ifelse(samples[iter,9]>pmix0,1,samples[iter,9]/pmix0)
        # inverted logistic
        samples[iter,2] <- ifelse(samples[iter,2]<pmix1,0,(samples[iter,2]-pmix1)/(1-pmix1))
        samples[iter,8] <- ifelse(samples[iter,8]<pmix1,0,(samples[iter,8]-pmix1)/(1-pmix1))
        # 2 point masses for p2 and p4
        tmp <- (c(samples[iter,4],samples[iter,6]) - pmix2[1])/(pmix2[2] - pmix2[1])
        if(samples[iter,4] < pmix2[1])
            tmp[1] <- 0
        if(samples[iter,4] > pmix2[2])
            tmp[1] <- 1

        if(samples[iter,6] < pmix2[1])
            tmp[2] <- 0
        if(samples[iter,6] > pmix2[2])
            tmp[2] <- 1

        samples[iter,4] <- tmp[1]
        samples[iter,6] <- tmp[2]
        # point mass along the diagonal
        if(max(samples[iter,4],samples[iter,6])==1)
            samples[iter,5] <- ifelse(samples[iter,5]>pmix3,1,samples[iter,5]/pmix3)

        # calculate eta
        tt = bezier_intersection(p0 = c(samples[iter,4],1),
                                 p1 = c(samples[iter,5],samples[iter,5]),
                                 p2 = c(1,samples[iter,6]), m = 1)
        if(tt==0 | tt == 1)
            eta=1
        if(tt>0 & tt <1)
            eta = bezier_eval(p0 = c(samples[iter,4],1),
                              p1 = c(samples[iter,5],samples[iter,5]),
                              p2 = c(1,samples[iter,6]), t = tt)

        samples[iter,11] <- eta

        # print results
        if(iter == 1 | iter %% print.every == 0){
            cat(iter,round(samples[iter,1:10],2),'// eta = ',eta,'\n')
            if(traceplot==T){
                par(mfrow=c(3,3))
                for(variable in c(1:9))
                    plot(samples[1:iter,variable],type='l',
                         xlab = var_names[variable],ylab='')
                par(mfrow=c(1,1))
            }
        }
    }
    return(samples)
}
