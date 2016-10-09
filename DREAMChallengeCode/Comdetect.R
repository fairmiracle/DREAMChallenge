#' Using exisiting algorithms on Disease Module Identification DREAM Challenge subchallenge1
#' Using recursive programming to find mode size to be 3~100
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
library(igraph)
library(Matrix)
recursiveigraph <- function(g, savefile, method = c('fastgreedy','louvain')){
    
    if (method == "fastgreedy")
        fc <- cluster_fast_greedy(g)
    else if (method == "louvain")
        fc <- cluster_louvain(g)
    #fc <- cluster_label_prop(g)         #label propagation algorithm 1min for large net
    #fc <- cluster_fast_greedy(g)        #fast greedy propagation algorithm    10min for large net
    #fc <- cluster_louvain(g)             #multi-level algorithm    1min for large net
    #fc <- cluster_infomap(g)            #InfoMAP algorithm    10min for large net
    #fc <- cluster_walktrap(g)

    memfc <- membership(fc)
    msize <- sizes(fc)
    
    if(length(msize) > 1){
        
    
    for (i in 1:length(msize)) {
        if(msize[i] < 100 & msize[i] >= 3){
            mgeneids <- V(g)$name[which(memfc==i)]
            mgeneids <- as.numeric(mgeneids) - 1
            cp = c()
            for (k in 1:length(mgeneids)){
                cp = paste(cp,mgeneids[k],sep='\t')
            }
            write(paste(1.00,cp,sep=''),file = savefile,append = TRUE)
        } else if (msize[i] > 100) {
            #large modules, in recursive way
            ids = which(memfc==i)
            g2 <- induced.subgraph(graph=g,vids=ids)
            recursiveigraph(g2,savefile,method)
        } else {
            next
        }
    }
    }
    else{
        print(paste('Atom! with size',length(V(g)),sep=' '))
        edges <- get.edgelist(g)
        write.table(edges,paste(savefile,'Atomsize',length(V(g)),sep='_'),sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE,append = TRUE)      
    }
}

excuteigraph <- function(filename,method = c('fastgreedy','louvain')){
    savefile=paste('result/',filename,sep='')
    net = read.delim(filename,header = FALSE)
    x = net[,1]+1
    y = net[,2]+1
    # W = sparseMatrix(i = x, j = y, x = net[,3], symmetric = TRUE)
    # W = sparseMatrix(i = x, j = y, x = net[,3],dims = c(max(max(x),max(y)),max(max(x),max(y))))
    # weighted network construction in igraph
    g = graph_from_edgelist(cbind(x,y),directed=F)
    V(g)$name=1:length(V(g))
    E(g)$weight = net[,3]
    recursiveigraph(g,savefile,method)
}

# Star net, only pick the first 40 nodes as the module
# normally larger than 100, no worries about size 40
forstargraph <- function(gt,savefile1){
    mgeneids <- V(gt)$name[1:40]
    mgeneids <- as.numeric(mgeneids) - 1
    cp = c()
    for (k in 1:length(mgeneids)){
            cp = paste(cp,mgeneids[k],sep='\t')
    }
    write(paste(1.00,cp,sep=''),file = savefile1,append = TRUE)
}

# W: the whole adjacency matrix
# modulefile: Atom file
# savefile1: new module ids file
ForAtom <- function(W,savefile1,modulefile){
    el = read.table(modulefile)
    gt = graph_from_edgelist(cbind(as.character(el[,1]),as.character(el[,2])),directed=F)
    predictedid=as.numeric(V(gt)$name)
    Wt = W[predictedid,predictedid]
    rownames(Wt) = predictedid
    N = dim(Wt)[1]

    if(ecount(gt)<vcount(gt)*1.2){
        print('for star net!')
        forstargraph(gt,savefile1)
    } else {
        source("AMOUNTAIN.R")
        for (k in 1:10) {
            abegin = 0.01
            aend = 0.9
            maxsize = 100
            minsize = 3
            for (i in 1:20) {
                x <- moduleIdentificationGPFixSSnoNs(Wt,rep(1/N,N),
                            a=(abegin+aend)/2,maxiter = 50)
                predictedid <- which(x[[2]]!=0)
                #predictedid <- getLargestComp(W,predictedid)
                if(length(predictedid) > maxsize){
                    abegin = (abegin+aend)/2
                }else if (length(predictedid) < minsize){
                    aend = (abegin+aend)/2
                }else
                    break
            }
            if(length(predictedid) > 100 ){
                print('for totally completed!')
                forTotalcompletegraph(gt,savefile1)
                break
            }
            if(length(predictedid) <= maxsize){
                tmpstr = as.character(as.numeric(rownames(Wt)[predictedid])-1)
                cp = c()
                for (j in 1:length(tmpstr)){
                    cp = paste(cp,tmpstr[j],sep='\t')
                }
                write(paste(1,cp,sep=''),file = savefile1,append = TRUE)
            }
            Wt = Wt[-predictedid,-predictedid]
            N = N - length(predictedid)

            if( N <3 || sum(Wt)==0)
                break
        }
    }
}
