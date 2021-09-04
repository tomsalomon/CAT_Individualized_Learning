#### FROM AKRAM (not everything works yet) ######
library(heplots)
labcex=1.5
axiscex=1.2
namecex=1.2
titlecex=1.8
pluscex=1.2
starcex=1.5
data2analyze$timechosen=NA
data2analyze$timeunchosen=NA
data2analyze$timechosen[data2analyze$Response=='b']=data2analyze$gaze_stim_l[data2analyze$Response=='b']
data2analyze$timechosen[data2analyze$Response=='y']=data2analyze$gaze_stim_r[data2analyze$Response=='y']
data2analyze$timeunchosen[data2analyze$Response=='b']=data2analyze$gaze_stim_r[data2analyze$Response=='b']
data2analyze$timeunchosen[data2analyze$Response=='y']=data2analyze$gaze_stim_l[data2analyze$Response=='y']

# data2analyze$numfixchosen[data2analyze$Response=='b']=data2analyze$numLeftFix[data2analyze$Response=='b']
# data2analyze$numfixchosen[data2analyze$Response=='y']=data2analyze$numRightFix[data2analyze$Response=='y']
# data2analyze$numfixunchosen[data2analyze$Response=='b']=data2analyze$numRightFix[data2analyze$Response=='b']
# data2analyze$numfixunchosen[data2analyze$Response=='y']=data2analyze$numLeftFix[data2analyze$Response=='y']

data2analyze$chose[data2analyze$IsleftGo==1&data2analyze$Response=='b']="Go"
data2analyze$chose[data2analyze$IsleftGo==0&data2analyze$Response=='y']="Go"
data2analyze$chose[data2analyze$IsleftGo==1&data2analyze$Response=='y']="NoGo"
data2analyze$chose[data2analyze$IsleftGo==0&data2analyze$Response=='b']="NoGo"
data2analyze$notchosen[data2analyze$IsleftGo==1&data2analyze$Response=='b']="NoGo"
data2analyze$notchosen[data2analyze$IsleftGo==0&data2analyze$Response=='y']="NoGo"
data2analyze$notchosen[data2analyze$IsleftGo==1&data2analyze$Response=='y']="Go"
data2analyze$notchosen[data2analyze$IsleftGo==0&data2analyze$Response=='b']="Go"


d=data.frame(subjectID=rep(levels(data2analyze$subjectID),2), PairType2=rep(levels(data2analyze$PairType2),2),time=c(data2analyze$timechosen,data2analyze$timeunchosen),ischosen=c(rep("chosen",length(data2analyze$timechosen)),rep("unchosen",length(data2analyze$timeunchosen))),whichchosen=c(data2analyze$chose,data2analyze$notchosen))

mod=lmer(time ~ ischosen * whichchosen + (1|subjectID), data=d, na.action=na.omit)

d2=subset(d,d$ischosen=="unchosen")
mod2=lmer(time ~ whichchosen + (1|subjectID), data=d2, na.action=na.omit) #simple effect within unchosen

m=cbind(tapply(data2analyze$timechosen,data2analyze$chose,mean),tapply(data2analyze$timeunchosen,data2analyze$notchosen,mean))
mchosen=aggregate(timechosen~subjectID+chose,data=data2analyze,mean)
munchosen=aggregate(timeunchosen~subjectID+notchosen,data=data2analyze,mean)
munchosen1=munchosen
mchosen1=mchosen
se=cbind(c(sd(mchosen1$timechosen[mchosen1$chose=="Go"]),sd(mchosen1$timechosen[mchosen1$chose=="NoGo"])),c(sd(munchosen1$timeunchosen[munchosen1$notchosen=="Go"]),sd(munchosen1$timeunchosen[munchosen1$notchosen=="NoGo"])))
colnames(m) = c('Chosen','Not Chosen')
par(mar=c(4,3,5,0.5), mgp=c(1.5,0.4,0),oma=c(.5,.5,.5,.5))
xvals=barplot(m,beside=T,space=c(.1,.1,.5,.1),ylim=range(0,.7),col=rep(c("forestgreen","firebrick3"),2),border=NA,ylab="Proportion total choice time eye on item",main=paste("Proportion of total choice time eye on item \n during Go vs. NoGo probe trials", sep=""), cex.lab=labcex-.1,cex.main=titlecex,axes=F)
Axis(side=2, labels=c("","0.1","0.2","0.3","0.4","0.5","0.6","1.0"),at=c(0,.1,.2,.3,.4,.5,.6,.7),cex.axis=axiscex)
axis.break(2,style="zigzag",breakpos=.65)
#errbar(xvals,m,m+se,m-se,pch="",add=T,lwd=2)
# lines(xvals[3:4],rep(.2,2),lwd=2)
# text(mean(xvals[3:4]),.23,'***',cex=starcex)
mtext("Eyes on", at=-1,side=1,line=.4,cex=1.2)
mtext("Go", at=xvals[1],side=1,line=.4,cex=1.2)
mtext("NoGo", at=xvals[2],side=1,line=.4,cex=1.2)
mtext("Go", at=xvals[3],side=1,line=.4,cex=1.2)
mtext("NoGo", at=xvals[4],side=1,line=.4,cex=1.2)
mtext("When", at=-1,side=1,line=1.75,cex=1.5)
mtext("Chosen", at=mean(xvals[1:2]),side=1,line=1.75,cex=1.5)
mtext("Not Chosen", at=mean(xvals[3:4]),side=1,line=1.75,cex=1.5)