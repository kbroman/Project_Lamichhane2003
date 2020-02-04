plot.randomness <-
function(rule=0.8,tas=taloc,breaks=seq(0,100,by=10))
{
#  require(modreg)
  temp <- tas[tas[,5]<=rule,]
  noness <- tapply(temp[,4],temp[,1],sum)
  tas <- tas[!is.na(match(tas[,1],names(noness[noness>0]))),]
  tas[tas[,4]==2,4] <- 1
  
  tas <- tas[tas[,5]<=rule,]
  breaks <- breaks[breaks <= rule*100]

  p <- cut(tas[,5]*100,breaks=breaks)
  n <- tapply(tas[,4],p,length)
  val <- tapply(tas[,4],p,sum)/n*100
#  return(list(val,n,tas))
  
  oldlas <- par("las")
  oldmar <- par("mar")
  on.exit(par(las=oldlas,mar=oldmar))
  par(las=1,mar=c(5.1,4.1,6.1,2.1))

  x <- round(n*val/100)
  y <- n-x
  pval <- chisq.test(rbind(x,y))$p.value

  barplot(val,xaxt="n",ylab="Percent hits",xlab="TA location (% gene length)",
          col="white",width=1,space=0,
          main = paste(rule*100,"% Rule",sep=""))
  mtext(paste("P-value = ", round(pval*100), "%",sep=""),side=3,line=1.5)
  u <- par("usr")
  nbr <- length(breaks)-1
  segments(0,0,nbr,0)
  segments(0:nbr,u[3]-diff(u[3:4])*0.02,
           0:nbr,0,xpd=TRUE)
  text(0:nbr,u[3]-diff(u[3:4])*0.07,as.character(breaks),xpd=TRUE)
#  text((1:nbr)-0.5,val-0.2,n,xpd=TRUE)
#  text((1:nbr)-0.5,val+0.2,round(n*val/100),xpd=TRUE)

  barplot(n,xaxt="n",ylab="No. TA sites",xlab="TA location (% gene length)",
          col="white",width=1,space=0,
          main = paste(rule*100,"% Rule",sep=""))
  u <- par("usr")
  nbr <- length(breaks)-1
  segments(0,0,nbr,0)
  segments(0:nbr,u[3]-diff(u[3:4])*0.02,
           0:nbr,0,xpd=TRUE)
  text(0:nbr,u[3]-diff(u[3:4])*0.07,as.character(breaks),xpd=TRUE)
#  text((1:nbr)-0.5,val-0.2,n,xpd=TRUE)
#  text((1:nbr)-0.5,val+0.2,round(n*val/100),xpd=TRUE)

  print(c(sum(noness>0),sum(n),sum(val*n/100)))
  
}
