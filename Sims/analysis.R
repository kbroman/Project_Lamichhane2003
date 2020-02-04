#source("gibbs.R")

dat <- read.table("speciesex.txt", header=TRUE)

unix.time(out0 <- gibbs(dat[,2],dat[,3],trace=FALSE))
unix.time(out1 <- gibbs(dat[,2],dat[,3],trace=FALSE,start="all"))
unix.time(out2 <- gibbs(dat[,2],dat[,3],trace=FALSE,start="random"))

for(i in 1:10) { out[[i]] <- gibbs(dat[,2],dat[,3],trace=FALSE); cat(i,"\n") }
me <- sapply(out,function(a) mean(a$tot))
qu <- sapply(out,function(a) quantile(a$tot,c(0.025,0.975)))
plot(1,0,xlim=c(0.5,10.5),ylim=c(0,nrow(dat)),type="n",
     xlab="Run",ylab="No. viable genes")
abline(h=mean(me),lty=2,col="gray40")
segments(1:10,qu[1,],1:10,qu[2,],lwd=2)
points(1:10,me,lwd=2)

# inferred viable genes
out <- gibbs(dat[,2],dat[,3],return=TRUE)
probs <- apply(out$output,2,mean)
N <- round(out$summary[1])-sum(dat[,3]>0)
viable <- as.numeric(dat[,3]>0)
notobs <- (1:length(viable))[viable==0]
viable[notobs[rev(order(probs))[1:N]]] <- 1
rm(out,probs,N,notobs)

simcounts <- sim.data(dat[,2],(dat[,3]>0),sum(dat[,3]))

# simulation results, assuming 639 viable genes
n.sim <- 200
simres1 <- matrix(ncol=5,nrow=n.sim)
colnames(simres1) <- c("no.obs","mean","sd","2.5%","97.5%")
for(i in 1:n.sim) {
  simd <- sim.data(dat[,2],(dat[,3]>0),sum(dat[,3]))
  simres1[i,] <- c(sum(simd>0),gibbs(dat[,2],simd,trace=FALSE)$summary)
  cat(i,"\n")
}

# simulation results, assuming 2787 viable genes
simres2 <- matrix(ncol=5,nrow=n.sim)
colnames(simres2) <- c("no.obs","mean","sd","2.5%","97.5%")
for(i in 1:n.sim) {
  simd <- sim.data(dat[,2],viable,sum(dat[,3]))
  simres2[i,] <- c(sum(simd>0),gibbs(dat[,2],simd,trace=FALSE)$summary)
  cat(i,"\n")
}


# simulation results, assuming 2787 viable genes, with 2000 mutants observed
simres3 <- matrix(ncol=5,nrow=n.sim)
colnames(simres3) <- c("no.obs","mean","sd","2.5%","97.5%")
for(i in 1:n.sim) {
  simd <- sim.data(dat[,2],viable,2000)
  simres3[i,] <- c(sum(simd>0),gibbs(dat[,2],simd,trace=FALSE)$summary)
  cat(i,"\n")
}

plot(0,0,xlim=c(0.5,200.5),ylim=c(0,nrow(dat)),type="n",
     xlab="Simulation replicate",ylab="No. viable genes")
abline(h=2787)
segments(1:200,simres2[,4],1:200,simres2[,5])
segments(1:200,simres3[,4],1:200,simres3[,5],col="blue")
#points(1:200,simres2[,1],col="black")
#points(1:200,simres3[,1],col="blue")

plot(0,0,xlim=c(0.5,200.5),ylim=c(0,nrow(dat)),type="n",
     xlab="Simulation replicate",ylab="No. viable genes")
abline(h=639)
segments(1:200,simres1[,4],1:200,simres1[,5])
#points(1:200,simres1[,1],col="black")
