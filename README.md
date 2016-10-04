# Disease Module Identification DREAM Challenge

The R Code for ShanHeLab, see https://www.synapse.org/#!Synapse:syn7136340/wiki/403063

Depends: R (>= 3.1.0)
Files
==============
- Comdetect.R: Recursive community detection based on igraph package.
- AMOUNTAIN.R: Modified from AMOUNTAIN package.
- WeightedMethod.R: Multiple modules extraction using AMOUNTAIN.
- util.R: Some additional functions such as rank modules

## Usage
==============
Put the folder "DREAMChallengeCode" at the same level of data folder for both subchallenges. Make a new folder "result" to save the modules, and two subfolders named "result/subchallenge1", "result/subchallenge2".
### Subchallenge 1
```{r}
source('DREAMChallengeCode/Comdetect.R')
filename='subchallenge1/1_ppi_anonym_v2.txt'
excuteigraph(filename,'fastgreedy')
```
It taks a while for network 1 using fast-greedy. And also make sure the memory is large enough for holding large matrix and deep recursive procedure. Then the preliminary result would save under "result/subchallenge1". We can see four modules are larger then expected, which needs further division.
```{r}
net = read.delim(filename,header = FALSE)
x = net[,1]+1
y = net[,2]+1
W = sparseMatrix(i = x, j = y, x = net[,3], symmetric = TRUE)
savefile=paste('result/',filename,sep='')
ForAtom(W,savefile,'result/subchallenge1/1_ppi_anonym_v2.txt_Atomsize_109')
ForAtom(W,savefile,'result/subchallenge1/1_ppi_anonym_v2.txt_Atomsize_159')
ForAtom(W,savefile,'result/subchallenge1/1_ppi_anonym_v2.txt_Atomsize_338')
ForAtom(W,savefile,'result/subchallenge1/1_ppi_anonym_v2.txt_Atomsize_467')
```
Now all clusters(modules) are saved in "result/subchallenge1/1_ppi_anonym_v2.txt". We need to assign the module score (defined as the sum of edges weights) to each module and rank them from highest to lowest. In previous submission we selected part of them to avoid too many false postive ones.
```{r}
source('DREAMChallengeCode/util.R')
modulesid('result/subchallenge1/1_ppi_anonym_v2.txt')
reassignScore(filename,'result/subchallenge1/1_ppi_anonym_v2.txt')
```
### Subchallenge 2
```{r}
library(Matrix)
source('DREAMChallengeCode/Comdetect.R')
returnW <- function(filename,dims){
	net = read.delim(filename,header = FALSE)
	x = net[,1]+1
	y = net[,2]+1
	adj = sparseMatrix(i = x, j = y, x = net[,3], symmetric = TRUE, dims = dims)
	#adj = as.matrix(adj)
}
W1 = returnW('subchallenge2/1_ppi_anonym_aligned_v2.txt',c(21115,21115))
W2 = returnW('subchallenge2/2_ppi_anonym_aligned_v2.txt',c(21115,21115))
W4 = returnW('subchallenge2/4_coexpr_anonym_aligned_v2.txt',c(21115,21115))
W5 = returnW('subchallenge2/5_cancer_anonym_aligned_v2.txt',c(21115,21115))
W6 = returnW('subchallenge2/6_homology_anonym_aligned_v2.txt',c(21115,21115))
W6 = W6/max(W6)
filename='subchallenge2/3_signal_anonym_aligned_directed_v3.txt'
net = read.delim(filename,header = FALSE)
x = net[,1]+1
y = net[,2]+1
W3 = sparseMatrix(i = x, j = y, x = net[,3],dims = c(21115,21115))
W3[lower.tri(W3)] = t(W3)[lower.tri(W3)]
W3 = W3/max(W3)

W = W1+W2+W3+W4+W5+W6
W[W<0.5] = 0
g = graph_from_adjacency_matrix(W,mode='undirected',weighted=TRUE)
V(g)$name=1:length(V(g))
recursiveigraph(g,'result/subchallenge2/sub2.txt',method)
```
