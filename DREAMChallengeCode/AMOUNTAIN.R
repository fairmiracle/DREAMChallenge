#' Module Identification
#' 
#' Algorithm for Module Identification on single network
#' 
#' @param W edge score matrix of the network, n x n matrix
#' @param z node score vector of the network, n-length vector
#' @param x0 initial solution, n-length vector
#' @param a parameter in elastic net the same as in \code{\link{EuclideanProjectionENNORM}}
#' @param lambda parameter in objective, coefficient of node score part
#' @param maxiter maximal interation of whole procedure
#' 
#' @return a list containing function objective vector and the solution 
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @references MOUNTAIN
#' @seealso \code{\link{EuclideanProjectionENNORM}}
#' @keywords module identification
#' 
#' @examples
#' n = 100
#' k = 20
#' theta = 0.5
#' pp <- networkSimulation(n,k,theta)
#' moduleid <- pp[[3]]
#' ## use default parameters here
#' x <- moduleIdentificationGPFixSS(pp[[1]],pp[[2]],rep(1/n,n))
#' predictedid<-which(x[[2]]!=0)
#' recall <- length(intersect(predictedid,moduleid))/length(moduleid)
#' precise <- length(intersect(predictedid,moduleid))/length(predictedid)
#' Fscore <- (2*precise*recall/(precise+recall))
#' @export
#' 
moduleIdentificationGPFixSS <- function(W,z,x0,a=0.5,lambda=1,maxiter=1000){
    x = x0
    epsilon = 1e-6
    grad = -W%*%x-lambda*z
    f_x = -0.5*t(x)%*%(W%*%x)-lambda*(t(z)%*%x)
    func = numeric(maxiter)

    for (iteration in 1:maxiter){
            #y = x-1*grad
            #print(sum(y)+0.5*gamma*sum(y*y))
            func[iteration] = f_x
            #x_cand = EuclideanProjectionEN(x-1*grad,t=1,alpha = a)
            #x_cand = EuclideanProjection(x-1*grad,t=radius)
            x_cand = EuclideanProjectionENNORM (x-1*grad,t=1,alpha = a)
            if(sum(abs(x_cand-x)^2)^(1/2) < epsilon){break}
            x=x_cand
            grad = -W%*%x-lambda*z
            f_x = -0.5*t(x)%*%(W%*%x)-lambda*(t(z)%*%x)
            
    }
    return (list(func[1:(iteration-1)],x))
}
#' Module Identification
#' 
#' Algorithm for Module Identification on single network without node score
#' 
moduleIdentificationGPFixSSnoNs <- function(W,x0,a=0.5,maxiter=1000){
    x = x0
    epsilon = 1e-6
    grad = -W%*%x
    f_x = -0.5*t(x)%*%(W%*%x)
    func = numeric(maxiter)
    
    xlist=list()
    x_candlist=list()
    for (iteration in 1:maxiter){
        func[iteration] = f_x[1,1]
        xlist[[iteration]] = x
        x_cand = EuclideanProjectionENNORM (x-1*grad,t=1,alpha = a)
        x_candlist[[iteration]] = x_cand
        if(sum(abs(x_cand-x)^2)^(1/2) < epsilon){break}
        x=x_cand
        grad = -W%*%x
        f_x = -0.5*t(x)%*%(W%*%x)
        
    }
    return (list(func[1:(iteration-1)],x))
}

#' Euclidean projection on elastic net
#' 
#' Piecewise root finding algorithm for Euclidean projection on elastic net
#' 
#' @param y constant vector
#' @param t radius of elastic net ball
#' @param alpha parameter in elastic net: alpha x_1 + (1-alpha)*x_2^2=t
#' 
#' @return a list containing network adjacency matrix, node score and module membership 
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @references Gong, Pinghua, Kun Gai, and Changshui Zhang. "Efficient euclidean projections via piecewise root finding and its application in gradient projection." Neurocomputing 74.17 (2011): 2754-2766.
#' @keywords Euclidean projection
#' 
#' @examples
#' y=rnorm(100)
#' x=EuclideanProjectionENNORM(y,1,0.5)
#' sparistyx = sum(x==0)/100
#' @export
#' 
EuclideanProjectionENNORM <- function(y,t,alpha = 0.5){
    #f(theta) = \sum max(0,y_i-theta) - t
    maxiter = 100
    epsilon = 1e-9
    n = length(y)
    #theta_old = y[n]
    theta = 0
    x = numeric(n)
    if(alpha*sum(y)+(1-alpha)*sum(y*y) <= t){
        x = y
        return (x)
    } else {
        tmpx = (y-theta*alpha)/(1+(2-2*alpha)*theta)
        S = which(tmpx >= 0)
        a = (alpha^3-alpha^2)*length(S)-4*t*(1-alpha)^2
        b = -alpha^2*length(S)-4*t*(1-alpha)
        c = alpha*sum(y[S])+(1-alpha)*sum(y[S]*y[S])-t
        for (i in 1:maxiter) {
            theta = (-b-sqrt(b^2-4*a*c))/(2*a)
            # doesn't work for (-b+sqrt(b^2-4*a*c))/(2*a) why
            tmpx = (y-theta*alpha)/(1+(2-2*alpha)*theta)
            S = which(tmpx >= 0)
            a = (alpha^3-alpha^2)*length(S)-4*t*(1-alpha)^2
            b = -alpha^2*length(S)-4*t*(1-alpha)
            c = alpha*sum(y[S])+(1-alpha)*sum(y[S]*y[S])-t        
            if (a*theta^2+b*theta+c < epsilon)
                break
            #theta_old = theta
        }
    }
    for (i in 1:n) {
        x[i] = max(0,(y[i] - theta*alpha)/(1+(2-2*alpha)*theta))
    }
    return (x)
}
#' Module Identification for multi-layer network
#' 
#' Algorithm for Module Identification on multi-layer network sharing the same 
#' set of nodes (aligned), without inter-layer interactions. The result is one module 
#' across these layers.
#' 
#' @param listW a list of edge score matrices, each reprents edge score for one layer
#' @param listz a list of node score vectors, each reprents node score for one layer
#' @param lambda parameter in objective, coefficient of node score part
#' @param x0 initial solution, n-length vector
#' @param a parameter in elastic net the same as in \code{\link{EuclideanProjectionENNORM}}
#' @param maxiter maximal interation of whole procedure
#' 
#' @return a vector represents membership of module nodes
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @references AMOUNTAIN
#' @seealso \code{\link{vecconsensus}}
#' @keywords module identification, many-layer
#' 
#' @examples
#' Module Identification for multi-layer network without node scores
#' 
#' Algorithm for Module Identification on multi-layer network sharing the same 
#' set of nodes (aligned), without inter-layer interactions. The result is one module 
#' across these layers. Another simplified version of \code{\link{moduleIdentificationGPFixSSManylayer}}.
#' 
#' @param listW a list of edge score matrices, each reprents edge score for one layer
#' @param x0 initial solution, n-length vector
#' @param a parameter in elastic net the same as in \code{\link{EuclideanProjectionENNORM}}
#' @param maxiter maximal interation of whole procedure
#' 
#' @return a vector represents membership of module nodes
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @references AMOUNTAIN
#' @seealso \code{\link{moduleIdentificationGPFixSSManylayer}}
#' @keywords module identification, many-layer
#' 
#' @examples
moduleIdentificationGPFixSSManylayer3 <- function(listW,x0,a=0.5, maxiter=1000){
    numlayer = length(listW)
    x = x0
    epsilon = 1e-6
    W = listW[[1]]

    grad = -W%*%x
    f_x = -0.5*t(x)%*%(W%*%x)
    for (i in 2:numlayer){
        W = listW[[i]]
        grad = grad-W%*%x
        f_x = f_x-0.5*t(x)%*%(W%*%x)
    }
    
    func = numeric(maxiter)
    
    for (iteration in 1:maxiter){
        func[iteration] = f_x[1,1]
        x_cand = EuclideanProjectionENNORM (x-1*grad,t=1,alpha = a)
        if(sum(abs(x_cand-x)^2)^(1/2) < epsilon){break}
        x=x_cand
        
        W = listW[[1]]

        grad = -W%*%x
        f_x = -0.5*t(x)%*%(W%*%x)
        
        for (i in 2:numlayer){
            W = listW[[i]]
            grad = grad -W%*%x
            f_x = f_x-0.5*t(x)%*%(W%*%x)
        }
    }
    return (list(func[1:(iteration-1)],x))
}
