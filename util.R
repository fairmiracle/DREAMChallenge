#putting module id ahead
modulesid <- function(filename) {
	rl=readLines(filename)
	file.remove(filename)
	for (i in 1:length(rl)) {
		write(paste(i,rl[i],sep='\t'),file = filename,append = TRUE)
	}
}
#reading adjacency matrix
returnW <- function(filename,dims){
	net = read.delim(filename,header = FALSE)
	x = net[,1]+1
	y = net[,2]+1
	adj = sparseMatrix(i = x, j = y, x = net[,3], symmetric = TRUE, dims = dims)
}

#checking overlap and too large or too small modules, NULL means good
checking <- function(filename){
	rlines=readLines(filename)
	cp=list()
	for (i in 1:length(rlines)) {
	ap=strsplit(rlines[i],'\t')[[1]]
	ap=ap[3:length(ap)]
	cp[[i]] = ap
	}
	for (i in 1:(length(rlines)-1)) {
		for (j in (i+1):length(rlines)) {
			acp = intersect(cp[[i]],cp[[j]])
			if(length(acp) > 0)
			{
				print(c(i,j,min(length(cp[[i]]),length(cp[[j]]))))
				print(length(acp))
				print('\n')
			}
		}
	}

	minsize = 3
	rlines=readLines(filename)
	delid=c()
	for (i in 1:length(rlines)) {
		ap=strsplit(rlines[i],'\t')[[1]]
		ap=ap[3:length(ap)]
		if(length(ap)>100 | length(ap)<minsize)
			delid=c(delid,i)
	}
	print(delid)
}

# re-assign module socre
reassignScore <- function(filename,modulefile){
	net = read.delim(filename,header = FALSE)
	x = net[,1]+1
	y = net[,2]+1
	W = sparseMatrix(i = x, j = y, x = net[,3], symmetric = TRUE)

	rlines=readLines(modulefile)
	file.remove(modulefile)
	cp1=c()
	for (i in 1:length(rlines)) {
		ap=strsplit(rlines[i],'\t')[[1]]
		ap=ap[3:length(ap)]
		predictedid = as.numeric(ap) + 1
		mscore = sum(W[predictedid,predictedid])
		strtmp=c()
		for (j in 1:length(ap)) {
			strtmp = paste(strtmp,ap[j],sep='\t')
		}
		write(paste(paste(i,mscore,sep='\t'),strtmp,sep=''),file = modulefile,append = TRUE)
	}

	# Ordering the modules in decreasing
	cp1=c()
	rlines = readLines(modulefile)
	file.remove(modulefile)
	for(i in 1:length(rlines))
	{
		sp=strsplit(rlines[i],'\t')[[1]]
		cp1[i] = as.numeric(sp[2])
	}
	sortedsc=sort(cp1,decreasing=T,index.return =T)

 	for (i in 1:length(rlines)) {
		ap=strsplit(rlines[sortedsc$ix[i]],'\t')[[1]]
		ap=ap[2:length(ap)]
		strtmp = i
		for (j in 1:length(ap)) {
			strtmp = paste(strtmp,ap[j],sep='\t')
		}
		write(strtmp,file = modulefile,append = TRUE)
	}

}

# re-assign module socre when given adjacency matrix W
reassignScore2 <- function(W,modulefile){
	rlines=readLines(modulefile)
	file.remove(modulefile)
	cp1=c()
	for (i in 1:length(rlines)) {
		ap=strsplit(rlines[i],'\t')[[1]]
		ap=ap[3:length(ap)]
		predictedid = as.numeric(ap) + 1
		mscore = sum(W[predictedid,predictedid])
		strtmp=c()
		for (j in 1:length(ap)) {
			strtmp = paste(strtmp,ap[j],sep='\t')
		}
		write(paste(paste(i,mscore,sep='\t'),strtmp,sep=''),file = modulefile,append = TRUE)
	}

	# Ordering the modules in decreasing
	cp1=c()
	rlines = readLines(modulefile)
	file.remove(modulefile)
	for(i in 1:length(rlines))
	{
		sp=strsplit(rlines[i],'\t')[[1]]
		cp1[i] = as.numeric(sp[2])
	}
	sortedsc=sort(cp1,decreasing=T,index.return =T)

 	for (i in 1:length(rlines)) {
		ap=strsplit(rlines[sortedsc$ix[i]],'\t')[[1]]
		ap=ap[2:length(ap)]
		strtmp = i
		for (j in 1:length(ap)) {
			strtmp = paste(strtmp,ap[j],sep='\t')
		}
		write(strtmp,file = modulefile,append = TRUE)
	}

}

# re-assign module socre for network 3
reassignScore3 <- function(filename,modulefile){
	net = read.delim(filename,header = FALSE)
	x = net[,1]+1
	y = net[,2]+1
	W = sparseMatrix(i = x, j = y, x = net[,3],dims = c(max(max(x),max(y)),max(max(x),max(y))))

	rlines=readLines(modulefile)
	file.remove(modulefile)
	for (i in 1:length(rlines)) {
		ap=strsplit(rlines[i],'\t')[[1]]
		ap=ap[3:length(ap)]
		predictedid = as.numeric(ap) + 1
		mscore = sum(W[predictedid,predictedid])
		strtmp=c()
		for (j in 1:length(ap)) {
			strtmp = paste(strtmp,ap[j],sep='\t')
		}
		write(paste(paste(i,mscore,sep='\t'),strtmp,sep=''),file = modulefile,append = TRUE)
	}

	# Ordering the modules in decreasing
	cp1=c()
	rlines = readLines(modulefile)
	file.remove(modulefile)
	for(i in 1:length(rlines))
	{
		sp=strsplit(rlines[i],'\t')[[1]]
		cp1[i] = as.numeric(sp[2])
	}
	sortedsc=sort(cp1,decreasing=T,index.return =T)

 	for (i in 1:length(rlines)) {
		ap=strsplit(rlines[sortedsc$ix[i]],'\t')[[1]]
		ap=ap[2:length(ap)]
		strtmp = i
		for (j in 1:length(ap)) {
			strtmp = paste(strtmp,ap[j],sep='\t')
		}
		write(strtmp,file = modulefile,append = TRUE)
	}

}

# re-assign module socre when module score is already there
reassignScore4 <- function(modulefile){
	# Ordering the modules in decreasing
	cp1=c()
	rlines = readLines(modulefile)
	file.remove(modulefile)
	for(i in 1:length(rlines))
	{
		sp=strsplit(rlines[i],'\t')[[1]]
		cp1[i] = as.numeric(sp[2])
	}
	sortedsc=sort(cp1,decreasing=T,index.return =T)

 	for (i in 1:length(rlines)) {
		ap=strsplit(rlines[sortedsc$ix[i]],'\t')[[1]]
		ap=ap[2:length(ap)]
		strtmp = i
		for (j in 1:length(ap)) {
			strtmp = paste(strtmp,ap[j],sep='\t')
		}
		write(strtmp,file = modulefile,append = TRUE)
	}

}


# avg size
avgsize <-function(filename){
	cp=c()
	rlines=readLines(filename)
	for (i in 1:length(rlines)) {
		ap=strsplit(rlines[i],'\t')[[1]]
		cp[i] = length(ap)-2
	}

	return(mean(cp))
}

#module size
getsize <-function(filename){
	cp=c()
	rlines=readLines(filename)
	for (i in 1:length(rlines)) {
		ap=strsplit(rlines[i],'\t')[[1]]
		cp[i] = length(ap)-2
	}

	return(cp)
}


# module score
getscore <-function(filename){
	cp=c()
	rlines=readLines(filename)
	for (i in 1:length(rlines)) {
		ap=strsplit(rlines[i],'\t')[[1]]
		cp[i] = as.numeric(ap[2])
	}
	return(cp)
}
