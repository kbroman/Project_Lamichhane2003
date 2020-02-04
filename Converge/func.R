
write.mcmc <-
function(object,filestem="myfile",all=FALSE)
{
  file1 <- paste(filestem,".out",sep="")
  file2 <- paste(filestem,".ind",sep="")

  n <- length(object$n.essential)

  if(all && !is.na(match("output",names(object)))) {
  }
  else {
    write.table(cbind(1:n,object$n.essential),file=file1, quote=FALSE,
                row.names=FALSE,col.names=FALSE)
    write.table(cbind("thetaplus","1",paste(n)),file=file2,quote=FALSE,
                row.names=FALSE,col.names=FALSE)
  }
}
