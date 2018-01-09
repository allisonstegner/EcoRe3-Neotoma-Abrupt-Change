# FUNCTIONS________________________________________________________________

# Variance Index_______________________________________________________
#calculates the variance index for a moving window
# i.e., maximum eigenvalue of the covariance matrix from the dataset (Sundstrom)
# comm.matrix=transformed pollen data for the community with time bins as rows and taxa as cols
# window.size=number of samples of amount of time encompased by the moving window
# step.size=how far to slide the window in each iteration
# type=should window size be based on the the amount of time or the number of samples?

slidingVI<-function(comm.matrix,ages,window.size,step.size,type){
	bin.start<-c(ages[length(ages)],(ages[length(ages)]-window.size))
		max.bins<-(ceiling(ages[length(ages)]-ages[1])/step.size)-2
		VI<-c()
		
		midpoints<-c()
		nsamp<-c()
	if (type=="time") {
		bins<-matrix(NA,nrow=max.bins,ncol=2)
		for (i in 1:max.bins){
			bin.i<-bin.start-(i*step.size)
			bin.rows<-which(ages<=bin.i[1] & ages>=bin.i[2])
			T.pollen.wind<-comm.matrix[bin.rows,]
			nsamp[i]<-nrow(T.pollen.wind)
			VI[i]<-max(eigen(cov(T.pollen.wind))$values)	
			bins[i,]<-bin.i
			}	 
			midpoints<-bins[,2]+((bins[,1]-bins[,2])/2)
			out<-list(VIs=VI,bins=bins,midpoints=midpoints,nsamp=nsamp)
		} else if (type=="sample") {
			bins<-matrix(NA,nrow=(nrow(comm.matrix)-window.size),ncol=2)
			for (i in (1:(nrow(comm.matrix)-window.size))){
				T.pollen.wind<-comm.matrix[c(i:(i+(window.size-step.size))),]
				nsamp[i]<-nrow(T.pollen.wind)
				bins[i,]<-ages[c(i,(i+(window.size-step.size)))]
				VI[i]<-max(eigen(cov(T.pollen.wind))$values)
			}	
		} else {
			print("invalid type")
		}
	midpoints<-bins[,2]+((bins[,1]-bins[,2])/2)
	out<-list(VIs=VI,bins=bins,midpoints=midpoints,nsamp=nsamp)
	return(out)
}


#GAMs___________________________________________________________

plot.gam<-function(taxa,ages2, taxa.names){
	source("/Users/allisonstegner/Dropbox (Personal)/ACES/R files/Deriv.R")
	require(mgcv)

	dev.new(width=4,height=6)
	par(mfrow=c(length(taxa),1),mar=c(0,2,0,1),oma=c(4,4,2,1))

	for(p in 1:length(taxa)){
		sppm<-taxa[[p]]
		ages<-ages2[[p]]
		spp.data<-as.data.frame(cbind(sppm,ages))
		names(spp.data)<-c("sppm","ages")
		#m<-gam(sppm~s(ages,k=10,bs="ad"),method="REML",select=TRUE)
		m0<-gamm(sppm~s(ages,k=50),data=spp.data)
		#m1<-gamm(sppm~s(ages),data=spp.data,correlation=corARMA(form=~ages,p=1))
		#m2<-gamm(sppm~s(ages),data=spp.data,correlation=corARMA(form=~ages,p=2))
		#anova(m1$lme,m2$lme,m3$lme)
		mod<-m0$gam

		pdat<-with(spp.data,data.frame(ages=seq(min(ages),max(ages),length=200)))
		p0<-predict(mod,newdata=pdat)

		pred<-p0
		m.d<-Deriv(mod,n=200)

		pred1 <- with(spp.data, data.frame(ages = seq(min(ages), max(ages), length = 200)))
		pred1 <- cbind(pred1, as.data.frame(predict(mod, pred1, se.fit = TRUE, unconditional = TRUE)))
		pred1 <- transform(pred1,Fitted = fit,Upper = fit + (2 * se.fit),Lower = fit - (2 * se.fit),ages = ages)

		CI<-confint(m.d,alpha=0.01)
		S<-signifD(pred,m.d$ages$deriv,CI$ages$upper,CI$ages$lower,eval=0)

		if (p==length(taxa)){
			plot(sppm~ages,data=spp.data,type="p",pch=16,ylab="%pollen",xlim=c(10000,0),las=1)
			mtext("% pollen",2,line=3,cex=0.75)
		} else {
			plot(sppm~ages,data=spp.data,type="p",pch=16,ylab="%pollen",xlim=c(10000,0),xaxt="n",las=1)
			mtext("% pollen",2,line=3,cex=0.75)
		}
	
		abline(v=seq(0,10000,1000),col="gray80",cex=0.2,lwd=0.5)
		points(ages,sppm,pch=16,col="gray50")
		lines(pred1$ages,pred1$Upper,lty=3,lwd=2)
		lines(pred1$ages,pred1$Lower,lty=3,lwd=2)
		lines(pred1$ages,pred1$fit,lty=1,lwd=2)
		lines(S$incr~ages,data=pdat,lwd=3,col="blue")
		lines(S$decr~ages,data=pdat,lwd=3,col="red")	
		text(10000,max(sppm,na.rm=TRUE)*0.9,taxa.names[p],pos=4)
	}
	mtext(data.name[j],3,outer=TRUE)
}

# END FUNCTIONS__________________________________________________

#Lac Brule, pollen (Neotoma dataset id 19842)
#Irwin Smith Bog, pollen (dataset id 13047)
#Path Lake, pollen (dataset id 15209)
#Crevice, diatoms (dataset id 22099)
#Morrison, diatoms (dataset id 22196)

library(neotoma)
library(princurve)
library(rpart)

# pull site data from Neotoma
# requires internet connection

ids<-c(19842,13047,15209,22099,22196)
i=1
neotoma_dataset<-get_download(ids[i])
comm_data<-neotoma_dataset[[1]]$counts
ages<-neotoma_dataset[[1]]$chronologies[[1]]$age


dev.new(width=5,height=6)
par(mfrow=c(4,1),mar=c(1,1,1,1),oma=c(3,3,1,1))

# Principal curve_____________________________________________________
comm_curve<-principal.curve(comm_data)
plot(ages,comm_curve$lambda,xlim=c(max(ages),min(ages)),pch=16, cex=0.7,type="b")
text(max(ages),max(comm_curve$lambda)*0.95,"principal curves",pos=4)

# Variance Index_____________________________________________________
# function allows window to be defined by amt. of time or by number of samples
# inclined to think that defining window by number of samples is best
# may make sense to standardize by amount of time represented in each window

window=20
step=5
commVI<-slidingVI(comm_data,ages,window,step,"sample")
amt.time<-commVI$bins[,2]-commVI$bins[,1]
plot(commVI$bins[,1],commVI$VIs,xlim=c(max(ages),min(ages)),ylim=c(min(commVI$VIs),max(commVI$VIs)*1.2),pch=16,xaxt="s",type="b")
points(ages,rep(max(commVI$VIs)*1.2,length(ages)),pch="|",col="red")
#text(max(ages),max(commVI$VIs)*0.95,paste("windows=",window,"steps=",step),pos=4)
text(max(ages),max(commVI$VIs)*0.95,"Variance index",pos=4)

# VI standardized by the amt of time represented by each window
plot(commVI$bins[,1],commVI$VIs/amt.time,xlim=c(max(ages),min(ages)),ylim=c(min(commVI$VIs/amt.time),max(commVI$VIs/amt.time)*1.2),pch=16,xaxt="s",type="b")
points(ages,rep(max(commVI$VIs/amt.time)*1.2,length(ages)),pch="|",col="red")
#text(max(ages),max(commVI$VIs/amt.time)*0.95,paste("windows=",window,"steps=",step),pos=4)
text(max(ages),max(commVI$VIs/amt.time)*0.95,"Variance index, standardized by time",pos=4)


# BCART_____________________________________________________
rpart.test <- rpart(comm_curve$lambda~ages, method="anova")

layout(matrix(c(1,2,3,4),nrow=2,ncol=2))
par(mar=c(2,2,2,2),oma=c(2,2,2,2))
printcp(rpart.test)
plotcp(rpart.test)
rsq.rpart(rpart.test)
plot(rpart.test)
text(rpart.test, use.n=TRUE, all=TRUE, cex=.8)





# DLMs_____________________________________________________
source('/Users/allisonstegner/Dropbox/ACES/Simulation Project/R files/ODLMAR_Stegner_5Oct2017.R')


taxa<-comm_curve$lambda

dev.new(width=8,height=4)
par(mfrow=c(2,1),mar=c(0,1,0,1),oma=c(4,4,3,4))
plot(ages, (comm_curve$lambda),type="b",pch=16,cex=0.7,xaxt="n",las=1)
#par(mfrow=c(3,1),mar=c(0,1,0,1),oma=c(4,4,3,4))

taxon<-cbind(rev(ages),rev(taxa))
delta=0.9
nl=1
ODL.out = ODLMAR(nl,delta,taxon[,2],taxon[,1],"",FALSE)
T.ar<-taxon[1:(nrow(taxon)-nl),1]
	
plot(T.ar,ODL.out$LamdaMat[,"lamda"],type="l",lwd=2,col='blue', xlab='Time Steps',ylab='Eigenvalue +/- SE',ylim=c(-0.5,1.5),xlim=c(min(T.ar),max(T.ar)),las=1)
mtext("Cal YBP",1,2,cex=0.8)
abline(v=seq(0,14000,1000),col="gray70",lwd=0.5)
mtext("Eigenvalue",2,2.5,cex=0.8)

polygon(c((T.ar),rev(T.ar)),c((ODL.out$LamdaMat[,"lamda.plus"]),rev(ODL.out$LamdaMat[,"lamda.minus"])),col="deepskyblue",border=NA)
lines(T.ar,ODL.out$LamdaMat[,"lamda"],type='l',lwd=2,col='blue')
abline(h=1,lty=2)
	#text(min(ages),1.45,taxa.names[i],pos=2)



#_____________________________________________________



#Lac Brule, pollen (Neotoma dataset id 19842)
#Irwin Smith Bog, pollen (dataset id 13047)
#Path Lake, pollen (dataset id 15209)
#Crevice, diatoms (dataset id 22099)
#Morrison, diatoms (dataset id 22196)

neotoma_dataset<-get_download(22196)

ids<-c(19842,13047,15209,22099,22196)
taxa<-list()
ages<-list()
for (i in 1:length(ids)){
	neotoma_dataset<-get_download(ids[i])
	comm_data<-neotoma_dataset[[1]]$counts
	ages[[i]]<-neotoma_dataset[[1]]$chronologies[[1]]$age
	comm_curve<-principal.curve(comm_data)
	taxa[[i]]<-comm_curve$lambda
}

plot.gam(taxa,ages, as.character(ids))