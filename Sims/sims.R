n.viable <- seq(500,3500,by=500)
n.mut <- seq(500,2500,by=500)
n.sim <- 100
n.genes <- nrow(dat)

results <- array(dim=c(n.sim,5,length(n.viable),length(n.mut)))
dimnames(results) <- list(NULL,c("count>0","mean","sd","2.5%","97.5%"),
                          as.character(n.viable),as.character(n.mut))

for(i in seq(along=n.viable)) {
  for(j in seq(along=n.mut)) {
    for(k in 1:n.sim) {
      simdat <- sim.data(dat[,2], sample(c(rep(1,n.viable[i]),rep(0,n.genes-n.viable[i]))),
                         n.mut[j])
      out <- gibbs(dat[,2], simdat, trace=FALSE)
      results[k,,i,j] <- c(sum(simdat>0),out$summary)
    }
    print(c(i,j))
  }
}

      
      
