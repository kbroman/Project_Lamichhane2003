######################################################################
# OPERONS

# read data
gaps <- read.csv("../Data/gaps.csv") 
gaps$gap[gaps$gap < 0] <- 0  # overlap -> 0 gap
geneloc <- read.csv("../Data/mygeneloc.csv")

findop <-
function(gaps, mingap=0)
{
  n <- nrow(gaps)
  a <- (1:n)[gaps[,3] == 1 & gaps[,4] <= mingap]
  b <- diff(a)
  v <- list(unlist(gaps[a[1],1:2]))
  for(i in 2:length(a)) {
    if(b[i-1] == 1) v[[length(v)]] <- c(v[[length(v)]],unlist(gaps[a[i],2]))
    else v <- c(v,list(unlist(gaps[a[i],1:2])))
  }
  v
}

ops <- vector("list",21)
names(ops) <- paste(0:20)
for(i in 1:21) { ops[[i]] <- findop(gaps,i-1); cat(i,"\n") }

oplength <- cbind("1"=sapply(ops,function(a) 4250-length(unlist(a))),
                  t(sapply(ops,function(a)
                           table(factor(sapply(a,length),levels=2:10)))))



convertops <-
function(ops, genes=geneloc[,1])
{
  x <- rep(0,length(genes))
  for(i in 1:length(ops)) 
    x[match(ops[[i]], genes)] <- i
    
  x
}

operons <- convertops(ops[[1]])

# Create an "80% rule" data set, containing *all* 4250 genes.
attach("../R/.RData")
source("../../gibbs.R")
data80 <- taloc[taloc[,2]/taloc[,3] <= 0.8,]
genes80 <- unique(data80[,1])
data80 <- cbind(nsites=tapply(data80[,4],data80[,1],length),
                 counts=tapply(data80[,4],data80[,1],sum))
data80 <- data80[data80[,1] != 0,]
temp <- matrix(0,nrow=length(operons),ncol=2)
temp[match(genes80,geneloc[,1]),] <- data80
data80 <- temp;rm(genes80,temp)
detach(2)

simOPdata <-
function(n.sites=data80[,1], ops=operons, orient=geneloc[,4],
         n.mutants=sum(data80[,2]), prop.ess=0.4)
{
  n.genes <- length(n.sites)
  if(length(ops) != n.genes || length(orient) != n.genes)
    stop("n.sites, ops, and orient should all be the same length.")

  # simulate essential genes
  ess <- sample(1:n.genes)[1:ceiling(prop.ess*n.genes)]
  theta <- rep(1,n.genes) # theta = 1 if *not* essential
  theta[ess] <- 0 

  # use operon information to fill in additional theta=0
  n.ops <- max(ops)
  for(i in 1:n.ops) {
    wh <- (1:n.genes)[ops==i]
    if(orient[wh[1]]) wh <- rev(wh)

    for(j in 2:length(wh))
      if(theta[wh[j]] == 0)
        theta[wh[1:(j-1)]] <- 0
  }
    
  theta <- theta[n.sites>0]
  n.sites <- n.sites[n.sites>0]

  
  list(sum(theta),sim.mutants(n.sites,theta,n.mutants))
}

n.sim <- 10
outA <- outB <- outC <- matrix(nrow=n.sim,ncol=6)
colnames(outA) <- colnames(outB) <- colnames(outC) <- 
  c("sum.theta","n.unique","mean","sd","2.5%","97.5%")
n.sites <- data80[data80[,1]>0,1]
for(i in 1:n.sim) {
  dat <- simOPdata(prop.ess=0.5)
  outA[i,1:2] <- c(dat[[1]],sum(dat[[2]]>0))
  outA[i,-(1:2)] <- gibbs(n.sites,dat[[2]],n.mcmc=1000)$summary

  dat <- simOPdata(prop.ess=0.2)
  outB[i,1:2] <- c(dat[[1]],sum(dat[[2]]>0))
  outB[i,-(1:2)] <- gibbs(n.sites,dat[[2]],n.mcmc=1000)$summary

  dat <- simOPdata(prop.ess=0.8)
  outC[i,1:2] <- c(dat[[1]],sum(dat[[2]]>0))
  outC[i,-(1:2)] <- gibbs(n.sites,dat[[2]],n.mcmc=1000)$summary

  cat(i,"\n")
}
  
plot.out <-
function(out,add=FALSE)
{
  n <- nrow(out)
  par(las=1)
  if(!add) 
    plot(1:n,out[,1]/4207,ylim=c(0,1),xlim=c(0,n.sim)+0.5,
         xlab="Simulation replicate",ylab="Proportion non-essential",type="n")
  segments(1:n-0.2,out[,1]/4207,1:n+0.2,out[,1]/4207,col="red")
  segments(1:n-0.1,out[,3]/4207,1:n+0.1,out[,3]/4207)
  segments(1:n,out[,5]/4207,1:n,out[,6]/4207)
}

plot.out(outA)
plot.out(outB,add=TRUE)
plot.out(outC,add=TRUE)


#hist(gaps$gap[gaps$same.orient==1 & gaps$gap<30 & gaps$gap>0],breaks=seq(0.5,30.5,by=1))

#r <- seq(0,2*pi,length=1000)
#par(pty="s")
#plot(cos(r),sin(r),type="l",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
#segments(cos(pi/2),sin(pi/2),cos(pi/2),sin(pi/2)+0.05,lwd=2)
#phi <- 2*pi/4403837
#for(i in 2:nrow(geneloc)) {
#  if(geneloc[i,2] > geneloc[i-1,3]) {
#    r <- seq((geneloc[i,2]-1)*phi,(geneloc[i-1,3]-1)*phi,by = -0.5*phi)
#    lines(cos(r),sin(r),lwd=1,col="white")
#  }
#}
  
