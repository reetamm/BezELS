#' Plot the limit set based on MCMC output in polar coordinates and export median curve
#'
#' @param mcmc_samples Output object from mcmc_bezier function
#' @param x Original data in exponential margins
#' @param r_0_marg Vector of truncation thresholds for all points
#' @param above_thresh_marg Vector of TRUE/FALSE values, where TRUE indicates points selected for fitting the truncated Gamma model
#' @param theta True value of copula dependence parameter (if known)
#' @param copula True parametric copula form (if known)
#' @param thin.by Thin the posterior for plotting (and median curve calculation)
#' @param plot.fig TRUE/FALSE for whether the limit set should be plotted
#' @param plottitle Title of plot (only used if plot.fig = T)
#' @param plot.truth TRUE/FALSE for whether the true limit set should be plotted (if known)
#'
#' @return Row of thinned mcmc_samples which corresponds to the median curve
#' @export
#'
#' @examples
plot_bezier_polar       <- function(mcmc_samples, x, r_0_marg, above_thresh_marg,
                                    theta = NULL, copula = NULL, thin.by = 10, plot.fig = T,
                                    plottitle = 'Polar coordinates', plot.truth = T){

    r <- x[ ,1] + x[ ,2]
    w <- x[ ,1] / r

    n.data = nrow(x)
    n.samples = nrow(mcmc_samples)

    t <- seq(0, 3, length=100)

    sss = seq(1,n.samples,by=thin.by)
    n.post = length(sss)
    p0y <- (mcmc_samples[sss,1])
    p1x <- (mcmc_samples[sss,2])
    p1y <- (mcmc_samples[sss,3])
    p2x <- (mcmc_samples[sss,4])
    p3 <- (mcmc_samples[sss,5])
    p4y <- (mcmc_samples[sss,6])
    p5x <- (mcmc_samples[sss,7])
    p5y <- (mcmc_samples[sss,8])
    p6x <- (mcmc_samples[sss,9])

    bezier_points_polar = array(dim = c(n.post,100,2))
    mean_rbez = rep(NA,n.post)
    for(i in 1:n.post){
        ## BEZIER CONTROL POINTS
        p <- matrix(c(0,p0y[i],
                      p1x[i],  p1y[i],
                      p2x[i],1,
                      p3[i],  p3[i],
                      1,p4y[i],
                      p5x[i],  p5y[i],
                      p6x[i],0), nrow=7, ncol=2, byrow=TRUE)

        ## CREATE A 2D BEZIER SPLINE WITH 3, 2-DEGREE BEZIER CURVES
        bezier_points <- bezier::bezier(t=t, p=p[, 1:2], deg=2)
        rbez <- bezier_points[ ,1] + bezier_points[ ,2]
        wbez <- bezier_points[ ,1] / rbez
        w_sorted <- sort(wbez, index.return=TRUE)$ix
        bezier_points_polar[i,,] = cbind(wbez,rbez)
    }
    xout = seq(0,1,length=100)
    bezier_points_gridded = matrix(NA,n.post,100)
    for(i in 1:n.post)
        bezier_points_gridded[i,] = approx(bezier_points_polar[i,,1],bezier_points_polar[i,,2],xout)$y

    w_sorted <- sort(w, index.return=TRUE)$ix

    fbplots = fda::fbplot(t(bezier_points_gridded),xout,xlim=c(0,1),ylim=c(0,2),color=scales::alpha(rgb(0,0,0), 0.05),
                              xlab='w',ylab='r/log(n)',cex.lab=1.5,main=plottitle,plot = plot.fig,outliercol = 0)

    medcurve = fbplots$medcurve
    if(plot.fig==T){
        points(w, r/log(n.data), pch=20, main="Pseudo-polar Coordinates", col=above_thresh_marg+1,cex=.25)
        lines(w[w_sorted], r_0_marg[w_sorted]/log(n.data), col="gray80", lty=2, lwd=3)

        lines(xout,bezier_points_gridded[medcurve,],col='blue',lwd=2)
        lines(1, 1, xlim=c(0, 1), ylim=c(0, 1.75),
              xlab="W", ylab=expression(r/log(n)), main="", type="n")

        if(plot.truth == T){
            w1 <- seq(0, 1, length=200)
            if(copula=='g'){
                unscaled_gauge <- (w1 + (1-w1) - 2 * theta * sqrt(w1*(1-w1))) / (1 - theta^2)
            }
            if(copula=='l'){
                max_w <- apply(cbind(w1, 1-w1), 1, max)
                min_w <- apply(cbind(w1, 1-w1), 1, min)
                unscaled_gauge <- 1/theta * max_w + (1 - 1/theta) * min_w
            }
            if(copula=='il'){
                unscaled_gauge <- (w1^(1/theta) +
                                       (1-w1)^(1/theta))^(theta)
            }
            if(copula=='al'){
                max_w <- apply(cbind(w1, 1-w1), 1, max)
                min_w <- apply(cbind(w1, 1-w1), 1, min)
                unscaled_gauge <- max_w/theta + min_w*(1-1/theta)
                unscaled_gauge = apply(cbind(unscaled_gauge,1),1,min)
            }
            lines(w1, 1/unscaled_gauge, lwd=3, col="red")
        }
    }
    return(medcurve)
}

#' Plot the limit set based on MCMC output in Euclidean coordinates and export median curve
#' @importFrom scales alpha
#' @importFrom bezier bezier
#'
#' @param mcmc_samples Output object from mcmc_bezier function
#' @param x Original data in exponential margins
#' @param above_thresh_marg Vector of TRUE/FALSE values, where TRUE indicates points selected for fitting the truncated Gamma model
#' @param theta True value of copula dependence parameter (if known)
#' @param copula True parametric copula form (if known)
#' @param thin.by Thin the posterior for plotting - should ideally correspond to the thin.by parameter of the plot_bezier_polar function
#' @param medcurve Median curve output from the plot.bezier.polar function
#' @param plottitle Title of the plot
#' @param plot.truth TRUE/FALSE for whether the true limit set should be plotted (if known)
#'
#' @return
#' @export
#'
#' @examples
plot_bezier_euclidean   <- function(mcmc_samples, x, above_thresh_marg ,
                                    theta=NULL,copula=NULL, thin.by = 10, medcurve = NULL,
                                    plottitle = 'Euclidean coordinates', plot.truth = T){
    r <- x[ ,1] + x[ ,2]
    w <- x[ ,1] / r

    n.data = nrow(x)
    n.samples = nrow(mcmc_samples)
    t <- seq(0, 3, length=100)
    sss = seq(1,n.samples,by=thin.by)
    n.post = length(sss)
    p0y <- (mcmc_samples[sss,1])
    p1x <- (mcmc_samples[sss,2])
    p1y <- (mcmc_samples[sss,3])
    p2x <- (mcmc_samples[sss,4])
    p3 <- (mcmc_samples[sss,5])
    p4y <- (mcmc_samples[sss,6])
    p5x <- (mcmc_samples[sss,7])
    p5y <- (mcmc_samples[sss,8])
    p6x <- (mcmc_samples[sss,9])

    plot(x[, 1]/log(n.data), x[, 2]/log(n.data), xlim=c(0, 1.2), ylim=c(0, 1.2), pch=20,cex=0.3,
         asp = 1,col=above_thresh_marg+1, cex.lab=1.5,
         xlab=expression(X[1]/log(n)), ylab=expression(X[2]/log(n)), main=plottitle)
    segments(x0=0, x1=1, y0=1, y1=1, lwd=3, lty=2, col="gray80")
    segments(x0=1, x1=1, y0=0, y1=1, lwd=3, lty=2, col="gray80")

    for(l in 1:n.post){
        ## BEZIER CONTROL POINTS
        p <- matrix(c(0,p0y[l],
                      p1x[l],  p1y[l],
                      p2x[l],1,
                      p3[l],  p3[l],
                      1,p4y[l],
                      p5x[l],  p5y[l],
                      p6x[l],0), nrow=7, ncol=2, byrow=TRUE)

        ## CREATE A 2D BEZIER SPLINE WITH 3, 2-DEGREE BEZIER CURVES
        bezier_points <- bezier::bezier(t=t, p=p[, 1:2], deg=2)
        ## PLOT BEZIER SPLINE
        lines(bezier_points,col=scales::alpha(rgb(0,0,0), 0.1))
    }
    boxplot(mcmc_samples[,4], add = T,horizontal = T,at = 1.1,boxwex = .2, col='white', outline=F)
    boxplot(mcmc_samples[,6], add = T,horizontal = F,at = 1.1,boxwex = .2, col='white', outline=F)

    if(plot.truth==T){
        w_order <- sort(w, index.return=TRUE)$ix
        if(copula=='g'){
            unscaled_gauge <- (w + (1-w) - 2 * theta * sqrt(w*(1-w))) / (1 - theta^2)
        }
        if(copula=='l'){
            max_w <- apply(cbind(w, 1-w), 1, max)
            min_w <- apply(cbind(w, 1-w), 1, min)
            unscaled_gauge <- 1/theta * max_w + (1 - 1/theta) * min_w
        }
        if(copula=='il'){
            unscaled_gauge <- (w^(1/theta) +
                                   (1-w)^(1/theta))^(theta)
        }
        if(copula=='al'){
            max_w <- apply(cbind(w, 1-w), 1, max)
            min_w <- apply(cbind(w, 1-w), 1, min)
            unscaled_gauge <- max_w/theta + min_w*(1-1/theta)
            unscaled_gauge = apply(cbind(unscaled_gauge,1),1,min)
        }
        lines((w /unscaled_gauge)[w_order],
              ((1-w)/unscaled_gauge)[w_order], lwd=3, col="red")
    }
    if(!is.null(medcurve)){
        p <- matrix(c(0,p0y[medcurve],
                      p1x[medcurve],  p1y[medcurve],
                      p2x[medcurve],1,
                      p3[medcurve],  p3[medcurve],
                      1,p4y[medcurve],
                      p5x[medcurve],  p5y[medcurve],
                      p6x[medcurve],0), nrow=7, ncol=2, byrow=TRUE)

        ## CREATE A 2D BEZIER SPLINE WITH 3, 2-DEGREE BEZIER CURVES
        bezier_points <- bezier::bezier(t=t, p=p[, 1:2], deg=2)
        ## PLOT BEZIER SPLINE
        lines(bezier_points,col='blue',lwd=2)
    }
}
