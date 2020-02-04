coverage <- matrix(nrow=length(n.viable),ncol=length(n.mut))
dimnames(coverage) <- dimnames(results)[3:4]
for(i in seq(along=n.viable))
  for(j in seq(along=n.mut))
    coverage[i,j] <- sum(results[,4,i,j] <= n.viable[i] &
                         results[,5,i,j] >= n.viable[i])

rmserr <- coverage
for(i in seq(along=n.viable))
  for(j in seq(along=n.mut))
    rmserr[i,j] <- sqrt(mean((results[,2,i,j] -n.viable[i])^2))

aveintlen <- coverage
for(i in seq(along=n.viable))
  for(j in seq(along=n.mut))
    aveintlen[i,j] <- mean(results[,5,i,j] - results[,4,i,j])
