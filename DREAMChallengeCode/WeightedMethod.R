####################################################
##	from weighted edgelist file to modules 
##	require package Matrix, igraph and AMOUNTAIN
####################################################
source('AMOUNTAIN.R')
ModulesAmoutain <- function(W,Nmodule,savefile,minsize,maxsize){
	saveAtomfile = paste(savefile,'Atom',sep = '')
	N = dim(W)[1]
	GeneNames = 1:N

	for (ii in 1:Nmodule) {
		abegin = 0.01
		aend = 0.9
		for (i in 1:20) {
			x <- moduleIdentificationGPFixSSnoNs(W,rep(1/N,N),
		                a=(abegin+aend)/2,maxiter = 50)
			#x <-moduleIdentificationGPNMFixSS(W,rep(1/N,N),
		    #            a=(abegin+aend)/2,maxiter = 50)
			predictedid <- which(x[[2]]!=0)
			if(length(predictedid) > maxsize){
				abegin = (abegin+aend)/2
			}else if (length(predictedid) < minsize){
				aend = (abegin+aend)/2
			}else
				break
		}

		if(length(predictedid) <= maxsize){
			modulescore = sum(W[predictedid,predictedid])
		
			tmpstr = as.character(as.numeric(GeneNames[predictedid])-1)
			cp = c()
			for (j in 1:length(tmpstr)){
				cp = paste(cp,tmpstr[j],sep='\t')
			}
			write(paste(modulescore,cp,sep=''),file = savefile,append = TRUE)
		}
		else {
			modulescoreW = W[predictedid,predictedid]
			print(paste('Atom! with size',length(predictedid),sep=' '))
			tmpstr = as.numeric(GeneNames[predictedid])-1
			forTotalcompletegraph(tmpstr,modulescoreW,saveAtomfile)
        }
		W = W[-predictedid,-predictedid]
		GeneNames = GeneNames[-predictedid]
		N = length(GeneNames)
		print(paste('Finishing module ',ii,sep=''))

		if(N < 3 | sum(W)==0)
			break
	}
}

ModulesAmoutain2 <- function(W,Nmodule,previousid,savefile){
	#savefile = paste('SubR1/',filename,sep = '')
	saveAtomfile = paste(savefile,'Atom',sep = '')
	N = dim(W)[1]
	
	GeneNames = 1:N

	W = W[-previousid,-previousid]
	GeneNames = GeneNames[-previousid]
	N = N - length(previousid)

	for (ii in 1:Nmodule) {
		abegin = 0.01
		aend = 0.9
		maxsize = 100
		minsize = 3
		for (i in 1:20) {
			x <- moduleIdentificationGPNMFixSS(W,rep(1/N,N),
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
		if(length(predictedid) <= maxsize){
			modulescore = sum(W[predictedid,predictedid])
			tmpstr = as.character(as.numeric(GeneNames[predictedid])-1)
			cp = c()
			for (j in 1:length(tmpstr)){
				cp = paste(cp,tmpstr[j],sep='\t')
			}
			write(paste(modulescore,cp,sep=''),file = savefile,append = TRUE)
		} else {
			modulescoreW = W[predictedid,predictedid]
			print(paste('Atom! with size',length(predictedid),sep=' '))
			tmpstr = as.numeric(GeneNames[predictedid])-1
			forTotalcompletegraph(tmpstr,modulescoreW,saveAtomfile)
        }

		W = W[-predictedid,-predictedid]
		GeneNames = GeneNames[-predictedid]
		N = length(GeneNames)
		print(paste('Finishing module ',ii,sep=''))

		if(N < 3 | sum(W)==0)
			break
	}
}

# a large clique, even AMOUNTAIN cannot devided it
forTotalcompletegraph <- function(predictedid,modulescoreW,savefile1){
    #el = read.table(modulefile)
    #g = graph_from_edgelist(cbind(as.character(el[,1]),as.character(el[,2])),directed=F)
    #predictedid=as.numeric(V(g)$name) - 1
    fixedlen = 60
    for (i in 1:(floor(length(predictedid)/fixedlen))) {
        cp = c()
        for (k in (fixedlen*(i-1)+1):(fixedlen*i)){
            cp = paste(cp,predictedid[k],sep='\t')
        }
        mscore = sum(modulescoreW[(fixedlen*(i-1)+1):(fixedlen*i),(fixedlen*(i-1)+1):(fixedlen*i)])
        write(paste(mscore,cp,sep=''),file = savefile1,append = TRUE)
    }

    if( (length(predictedid) - fixedlen*i)>=3){
    	cp = c()
    	for (k in (fixedlen*i+1):length(predictedid)){
            cp = paste(cp,predictedid[k],sep='\t')
    	}
    	mscore = sum(modulescoreW[(fixedlen*i+1):length(predictedid),(fixedlen*i+1):length(predictedid)])
    	write(paste(mscore,cp,sep=''),file = savefile1,append = TRUE)
	}
}

getLargestComp <- function (W, predictedid){
		require(igraph)
		g <- graph.adjacency(W[predictedid,predictedid],mode="undirected", weighted=TRUE)	
		comps <- decompose.graph(g,min.vertices=3)
		return (as.numeric(V(comps[[1]])$name))
}

ModulesAmoutainML <- function(savefile,listW,Nmodule, minsize=15, maxsize=100,amtiter=50){
	saveAtomfile = paste(savefile,'Atom',sep = '')
	W = listW[[1]]
	N = dim(W)[1]
	GeneNames = 1:N
	numlayer = length(listW)
	for (ii in 1:Nmodule) {
		abegin = 0.01
		aend = 0.9
		#maxsize = 100
		#minsize = 15
		for (i in 1:20) {
			x <- moduleIdentificationGPFixSSManylayer3(listW, rep(1/N,N), 
				a=(abegin+aend)/2, maxiter = amtiter)

			predictedid <- which(x[[2]]!=0)
			#predictedid <- getLargestComp(W,predictedid)
			if(length(predictedid) > maxsize){
				abegin = (abegin+aend)/2
			}else if (length(predictedid) < minsize){
				aend = (abegin+aend)/2
			}else
				break
		}
		if(length(predictedid) <= maxsize){
		modulescore = 0
		for (i in 1:numlayer){
            W = listW[[i]]
			modulescore = modulescore + sum(W[predictedid,predictedid])
		}
		tmpstr = as.character(as.numeric(GeneNames[predictedid])-1)
		cp = c()
		for (j in 1:length(tmpstr)){
			cp = paste(cp,tmpstr[j],sep='\t')
		}
		write(paste(modulescore,cp,sep=''),file = savefile,append = TRUE)
		} else {
			modulescorelistW = list()
			for (i in 1:numlayer){
            	W = listW[[i]]
				modulescorelistW[[i]] = W[predictedid,predictedid]
			}
			print(paste('Atom! with size',length(predictedid),sep=' '))
			tmpstr = as.numeric(GeneNames[predictedid])-1
			forTotalcompletegraph2(tmpstr,modulescorelistW,saveAtomfile)
        }

		for (i in 1:numlayer){
            W = listW[[i]]
			W = W[-predictedid,-predictedid]
			listW[[i]] = W
		}
		GeneNames = GeneNames[-predictedid]
		N = N - length(predictedid)
		print(paste('Finishing module ',ii,sep=''))
		#save(listW,GeneNames,N,file='sub2.RData')
	}
}

ModulesAmoutainML2 <- function(savefile,listW,Nmodule,previousid){
	saveAtomfile = paste(savefile,'Atom',sep = '')
	W = listW[[1]]
	N = dim(W)[1]
	GeneNames = 1:N
	numlayer = length(listW)

	for (i in 1:numlayer){
            W = listW[[i]]
			W = W[-previousid,-previousid]
			listW[[i]] = W
	}
	GeneNames = GeneNames[-previousid]
	N = N - length(previousid)

	for (ii in 1:Nmodule) {
		abegin = 0.01
		aend = 0.9
		maxsize = 100
		minsize = 20
		for (i in 1:20) {
			x <- moduleIdentificationGPFixSSManylayer3(listW, rep(1/N,N), 
				a=(abegin+aend)/2, maxiter = 20)

			predictedid <- which(x[[2]]!=0)
			#predictedid <- getLargestComp(W,predictedid)
			if(length(predictedid) > maxsize){
				abegin = (abegin+aend)/2
			}else if (length(predictedid) < minsize){
				aend = (abegin+aend)/2
			}else
				break
		}
		if(length(predictedid) <= maxsize){
			modulescore = 0
			for (i in 1:numlayer){
            	W = listW[[i]]
				modulescore = modulescore + sum(W[predictedid,predictedid])
			}
			tmpstr = as.character(as.numeric(GeneNames[predictedid])-1)
			cp = c()
			for (j in 1:length(tmpstr)){
				cp = paste(cp,tmpstr[j],sep='\t')
			}
			write(paste(modulescore,cp,sep=''),file = savefile,append = TRUE)
		} else {
			modulescorelistW = list()
			for (i in 1:numlayer){
            	W = listW[[i]]
				modulescorelistW[[i]] = W[predictedid,predictedid]
			}
			print(paste('Atom! with size',length(predictedid),sep=' '))
			tmpstr = as.numeric(GeneNames[predictedid])-1
			forTotalcompletegraph2(tmpstr,modulescorelistW,saveAtomfile)
        }
		for (i in 1:numlayer){
            W = listW[[i]]
			W = W[-predictedid,-predictedid]
			listW[[i]] = W
		}
		GeneNames = GeneNames[-predictedid]
		N = N - length(predictedid)
		print(paste('Finishing module ',ii,sep=''))
		#save(listW,GeneNames,N,file='sub2.RData')
	}
}

# a large clique, even AMOUNTAIN cannot devided it
forTotalcompletegraph2 <- function(predictedid,modulescorelistW,savefile1){
    #el = read.table(modulefile)
    #g = graph_from_edgelist(cbind(as.character(el[,1]),as.character(el[,2])),directed=F)
    #predictedid=as.numeric(V(g)$name) - 1
    fixedlen = 60
    numlayer = length(modulescorelistW)
    for (i in 1:(floor(length(predictedid)/fixedlen))) {
        cp = c()
        for (k in (fixedlen*(i-1)+1):(fixedlen*i)){
            cp = paste(cp,predictedid[k],sep='\t')
        }
        modulescore = 0
		for (j in 1:numlayer){
            modulescoreW = modulescorelistW[[j]]
			modulescore = modulescore + sum(modulescoreW[(fixedlen*(i-1)+1):(fixedlen*i),(fixedlen*(i-1)+1):(fixedlen*i)])
		}
        write(paste(modulescore,cp,sep=''),file = savefile1,append = TRUE)
    }

    if( (length(predictedid) - fixedlen*i)>=3){
    	cp = c()
    	for (k in (fixedlen*i+1):length(predictedid)){
            cp = paste(cp,predictedid[k],sep='\t')
    	}
    	modulescore = 0
		for (j in 1:numlayer){
            modulescoreW = modulescorelistW[[j]]
			modulescore = modulescore + sum(modulescoreW[(fixedlen*i+1):length(predictedid),(fixedlen*i+1):length(predictedid)])
		}
       	write(paste(modulescore,cp,sep=''),file = savefile1,append = TRUE)
	}
}

library(Matrix)
excutescrpt <- function(filename,N,savefile,minsize,maxsize){
	net = read.delim(filename,header = FALSE)
	x = net[,1]+1
	y = net[,2]+1
	W = sparseMatrix(i = x, j = y, x = net[,3], symmetric = TRUE)
	#savefile = paste('SubR1/',filename,sep = '')
	ModulesAmoutain(W,N,savefile,minsize,maxsize)
}
#knock down existing modules and extract new ones
excutescrpt2 <- function(filename,savefile,N,exisitingmodulefile){
	net = read.delim(filename,header = FALSE)
	x = net[,1]+1
	y = net[,2]+1
	adj = sparseMatrix(i = x, j = y, x = net[,3], symmetric = TRUE)

	rl = readLines(exisitingmodulefile)
	previousid = c()
	for (i in 1:length(rl)) {
		tmpid = as.numeric(strsplit(rl[i],'\t')[[1]])
		previousid = c(previousid,tmpid[3:length(tmpid)])
	}
	previousid = previousid+1
	stardmid = length(rl)
	#savefile = paste('SubR1/new',filename,sep = '')
	ModulesAmoutain2(adj,N,previousid,savefile)

	rl = readLines(savefile)
    file.remove(savefile)
    for (i in 1:length(rl)) {
        write(paste((stardmid+i),rl[i],sep='\t'),file = savefile,append = TRUE)
    }
}
