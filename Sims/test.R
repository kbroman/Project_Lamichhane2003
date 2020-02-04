n <- nrow(dat)
sites <- dat[,2]
viable <- rep(0,n);viable[1:sample(1:n,1)] <- 1;viable <- sample(viable)
counts <- table(factor(sample(1:n,sum(dat[,3]>0),repl=TRUE,prob=viable*sites),levels=1:n))

temp <- gibbs(sites,counts,n.mcmc=1000,burnin=100,skip=9,start="none",trace=FALSE)
hist(temp$total,breaks=seq(0,n+50,by=50),xlab="No. viable",main="")
abline(v=sum(viable),lwd=2,col="blue")
abline(v=sum(counts>0),lwd=2,col="red")

