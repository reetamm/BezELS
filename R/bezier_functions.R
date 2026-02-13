#' Find a point on a 2D Bezier curve B(t) at a particular value of t
#'
#' @param p0 Start control point
#' @param p1 Intermediate control point
#' @param p2 End control point
#' @param t: Location on the curve
#' @returns The x-coordinate of B(t)
#'
#' @examples bezier_eval(p0 = c(0.4,0.1),p1 = c(0.9,0.9),p2 = c(0.1,0.4),t = 0)
#'
bezier_eval = function(p0, p1, p2, t)
{
    x0  <- p0[1]
    x1  <- p1[1]
    x2  <- p2[1]
    bez <- t^2*(x0-2*x1+x2) - 2*t*(x0-x1) + x0
    return(bez)
}

#' Find the intersection of a 2D Bezier curve B(t) and a line passing through the origin with slope m
#'
#' @param p0 Start control point
#' @param p1 Intermediate control point
#' @param p2 End control point
#' @param m Slope of straight line passing through the origin
#'
#' @returns Value of t where the line intersects the curve B(t)
#' @export
#'
#' @examples bezier_intersection(p0 = c(0.4,0.1),p1 = c(0.9,0.9),p2 = c(0.1,0.4),m = 1)
bezier_intersection = function(p0, p1, p2, m)
{
    t <- matrix(0,length(m),2)
    if(p0[1] == p2[1] && p0[2] == p2[2]){
        return(t[,1])
    }
    aa   <- p0[2] - 2*p1[2] + p2[2] - m*(p0[1] - 2*p1[1] + p2[1])
    bb   <- -2*(p0[2]-p1[2] - m*(p0[1] - p1[1]))
    cc   <- p0[2] - m*p0[1]

    disc <- (bb^2 - 4*aa*cc)
    update = which(disc >= 0 & aa != 0)

    t[update,1] <- (-bb[update] - sqrt(disc[update]))/(2*aa[update])
    t[update,2] <- (-bb[update] + sqrt(disc[update]))/(2*aa[update])
    update2 = which(aa==0)
    t[update2,2] <- t[update2,1] <- c(-cc/bb)[update2]
    t = round(t,8)  ## newly added line to deal with really small values
    makezero = which(t[,1]<=0 | t[,1]>= 1)
    t[makezero,1] = 0
    makezero = which(t[,2]<=0 | t[,2]>= 1)
    t[makezero,2] = 0
    keept = apply(t,1,max)
    return(keept)
}
