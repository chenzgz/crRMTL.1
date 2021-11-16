\name{Simdata}
\alias{Simdata}
\docType{data}
\title{
 The simulated dataset under competing risks with 1000 observations.
}
\description{
A dataset constructed by simulation with 1000 observations on the following 3 variables.
}
\usage{data(Simdata)}
\format{time:the followed up time of patients. status: the status indicator, 0= censored, 1= event of interest, 2= competing risks; group: The group indicator for comparison, 0= control group, 1= treatment group.}
\examples{
library(crRMTL)
data(Simdata)
#test <- RMTLd_test(Simdata$time, Simdata$status, Simdata$group, alpha = 0.05, tau = 5)
# The Figure corresponding to the RMTLs of Simdata at tau = 5.
#par(mfrow=c(1,1),mai=c(1,1,.4,.4))
#fit<-cuminc(Simdata$time,Simdata$status,Simdata$group,cencode=0);round(fit$Tests,3)
#x1<-data.frame(t=fit$'0 1'$time,e=fit$'0 1'$est)
#x2<-data.frame(t=fit$'1 1'$time,e=fit$'1 1'$est)
#par(mfrow=c(1,1),mai=c(1,1,.4,.4))
#plot(x1$t,x1$e,lwd=2,'S',xlim=c(0,5),ylim=c(0,1),col='red',lab=c(12,5,1),cex.lab=1.3,lty=2,cex.axis=1,xlab='',ylab='',xaxt='n',yaxt='n',mgp=c(0.2,1,0),bty='l')
#segments(x0=5,y0=0,x1=5,y1=1.0,lty=2,lwd=2,col='gray')
#lines(x2$t,x2$e,col='black',lwd=2,'S',lty=1)
#axis(side=1,c(0,1,2,3,4,5),cex.axis=1.2,las=1.2,lwd=1.2,line=0,mgp=c(1,0.7,0),tck=-0.015)
#axis(side=2,c(0.0,0.2,0.4,0.6,0.8,1.0),cex.axis=1.2,las=1.2,lwd=1.2,mgp=c(1,0.7,0),tck=-0.015)
#tau<-expression(paste(tau, ' = 5'))
#S1<-expression(S[1])
#S2<-expression(S[2])
#RMTL1<-expression(paste(Treatment,' (', RMTL[1]: S[1] + S[2], ')'))
#RMTL2<-expression(paste(Placebo,' (', RMTL[0]: S[1], ')'))
#RMTLd<-expression(paste(RMTLd,' = 2.427 (95% CI: 1.242, 3.612) (',S[2],')'))
#star<-expression(paste('*RMTLd : the difference of RMTLs between groups'))
#title<-expression('The RMTL of event of interest in 'Simdata' dataset')
#mtext('Time',side=1,line=2.2,cex=1.3)
#mtext('Cumulative incidence',side=2,line=2.5,cex=1.3)
#mtext(tau,cex=1.3,side=1,line=-2,bty='n',at = 4.75)
#mtext(S1,cex=1.3,side=1,line=-6,bty='n',at = 1)
#mtext(S2,cex=1.3,side=1,line=-9,bty='n',at = 3)
#mtext(RMTLd,cex=1.3,side=1,line=-17.1,bty='n',at=7.5)
#mtext(star, cex=1,side=1,line=3.7,at = 1.2)
#mtext(title,cex=1.5,side=3,line=0,bty='n',at = 14)
#legend('topleft',c(RMTL1,RMTL2),lty=c(1,2),col=c(1,2),lwd=2,cex=1.3,bty='n')
}

