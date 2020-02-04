#attach("../R/.RData")
#data80 <- data80
#detach(2)
#save.image()

library(coda)

# run mcmc
#out1 <- negenes(data80[,1],data80[,2],data80[,3],data80[,4],n.mcmc=500000,skip=49,
#                burnin=500,startp=1)
#out2 <- negenes(data80[,1],data80[,2],data80[,3],data80[,4],n.mcmc=500000,skip=49,
#                burnin=500,startp=0)

# write to BUGS-like files
#write.mcmc(out1,"out1");write.mcmc(out2,"out2")

### Long runs from five starting points
# set.seed(2831342)
#load("../data80.RData")
#
#outA <- negenes(data80[,1],data80[,2],data80[,3],data80[,4],burnin=500,skip=49,
#                n.mcmc=50000,return=TRUE,startp=0)
#outB <- negenes(data80[,1],data80[,2],data80[,3],data80[,4],burnin=500,skip=49,
#                n.mcmc=500000,startp=0)

# combine information
outputA <- outputB <- vector("list",5)
names(outputA) <- names(outputB) <- paste(seq(0,1,by=0.25))
for(i in 1:5) {
  attach(paste("Work", i, "/.RData",sep=""))
  outputA[[i]] <- outA
  outputB[[i]] <- outB
  detach(2)
}
rm(i)

# gene-specific probs
prob <- sapply(outputA,function(a) 1-apply(a$output,2,mean))
prob <- prob[prob[,1] != 0,]
d <- matrix(0,5,5)
for(i in 1:4) for(j in (i+1):5) d[i,j] <- d[j,i] <- quantile(abs(prob[,i]-prob[,j]),0.9)
rm(i,j)
mp <- apply(prob,1,mean)
plot(tapply(prob[,1]-mp,cut(mp,quantile(mp,seq(0,1,by=0.01))),sd),ylim=c(0.005,0.025))
for(i in 2:5)
  points(tapply(prob[,i]-mp,cut(mp,quantile(mp,seq(0,1,by=0.01))),sd),
         col=c("black","blue","red","green","orange")[i])

# comparison of major results
round(t(sapply(outputB,function(a) a$summary))[,-2]/42.04,1)

# write in gibbs format
for(i in 1:5) write.mcmc(outputB[[i]],paste("out",i,sep=""))
