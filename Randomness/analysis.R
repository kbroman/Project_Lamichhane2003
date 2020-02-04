attach("../R/.RData")
taloc <- taloc
doubleta <- doubleta
detach(2)

temp1 <- paste(taloc[,1],taloc[,2],sep=":")
temp2 <- paste(c(doubleta[,1],doubleta[,3]),
               c(doubleta[,2],doubleta[,4]),sep=":")
taloc <- taloc[-match(temp2,temp1),]
rm(temp1,temp2,doubleta)

taloc <- cbind(taloc,prop=taloc[,2]/taloc[,3])


postscript("randomness.ps",horiz=TRUE,height=6.5,width=9)
plot.randomness(1,breaks=seq(0,100,by=5))
plot.randomness(0.8,breaks=seq(0,100,by=5))
dev.off()
