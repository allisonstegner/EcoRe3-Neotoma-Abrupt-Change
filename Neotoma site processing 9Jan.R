# Script to pull sites from the Neotoma Paleoecology Database
# And eliminate sites that don't meet minimum resolution and chronology requirements

library(neotoma)
library(maps)

# FUNCTIONS___________________________________________________
# has.chron______________________________________________________
# function to select ids for datasets that have a chonolorgy in Neotoma db
# tax_dec_dl is a Neotoma download object

has.chron<-function(tax_dec_dl){
	#chron.table<-c()
	#chron.type<-c()
	ind<-c()
	for (i in 1:length(tax_dec_dl)) {
		sitei<-tax_dec_dl[[i]]
		sitei.chron<-tax_dec_dl[[i]]$chronologies #all the chronologies for site i
		i.chrons<-names(sitei.chron)
		n.pollen.samples<-nrow(tax_dec_dl[[i]]$counts) #number of samples for site 
		if (i.chrons=="No chronology"){ #toss sites with no chronologies 
			ind[i]<-NA 
		} else {
			ind[i]<-sitei$dataset$dataset.meta$dataset.id
		}
	}
return(has.chron=ind)
}


#ts.min.length_____________________________________________________
# function to select ids for datasets that have number of data points >= minimum samples
# tax_dec_dl is a Neotoma download object
# min.samples is an integer: sites with number of samples less than min.samples are excluded

ts.min.length<-function(tax_dec_dl,min.samples){
	ind<-c()
	for (i in 1:length(tax_dec_dl)) {
		sitei<-tax_dec_dl[[i]]
		n.pollen.samples<-nrow(sitei$counts) #number of samples for site i
		if (n.pollen.samples<min.samples){ 
			ind[i]<-NA
		} else {
			ind[i]<-sitei$dataset$dataset.meta$dataset.id
		}
	}
		return(adequate.n=ind)
}
	
	
#select.high.res________________________________________________________
# function to select ids for datasets where number of years represented per pollen sample is less than max.grain
# tax_dec_dl is a Neotoma download object
# max.grain is the maximum allowable number of years represented per pollen sample

select.high.res<-function(tax_dec_dl,max.grain){ 
	chron.table<-c()
	chron.type<-c()
	for (i in 1:length(tax_dec_dl)) {
		sitei.chron<-tax_dec_dl[[i]]$chronologies #all the chronologies for site i
		i.chrons<-names(sitei.chron)
		n.pollen.samples<-nrow(tax_dec_dl[[i]]$counts) #number of samples for site i
			
		#develop a table summarizing various aspects including average resolution
		chron.data<-c()
		for (j in 1:length(i.chrons)){
			chronj<-sitei.chron[[i.chrons[j]]]
			if (nrow(chronj)<2){	next }
			i.names<-tax_dec_dl[[i]]$dataset$site.data$site.name
			n.dates<-nrow(chronj)
			chronj.dur<-max(chronj$age)-min(chronj$age)
			dates.per.time<-chronj.dur/n.dates
			yrs.per.pollen<-chronj.dur/n.pollen.samples
			chron.data1<-c(i.names,unique(chronj$dataset.id),names(sitei.chron)[j],n.dates,unique(chronj$age.type)[1],dates.per.time,yrs.per.pollen)
			chron.data<-rbind(chron.data,chron.data1)
			}
		chron.table<-rbind(chron.table,chron.data)
	}
	rownames(chron.table)<-c(1:nrow(chron.table))
	chron.table<-as.data.frame(chron.table)
	colnames(chron.table)<-c("dataset.name","dataset.id","chronology.type","n.dates","age.type","n.year.per.date","yrs.per.pollen")

	#choose sites where duration/n.samples is less than max.grain yrs
	resolution.include.exclude<-ifelse(as.numeric(as.vector(chron.table$yrs.per.pollen))>=max.grain,"exclude","include") 
	dated_tax<-cbind(chron.table,resolution.include.exclude)
	dated_tax<-dated_tax[dated_tax$resolution.include.exclude %in% "include",]
	dataset.list<-unique(as.character(dated_tax$dataset.id)) #list of dataset ids
	
	high.res.sites<-unique(as.numeric(as.character(dated_tax$dataset.id)))
	return(high.res.sites)
}


#min.chron.control________________________________________________
# function to select ids for datasets where number of chron controls is adequate
# tax_dec_dl is a Neotoma download object
# min.chron.grain is the minimum allowable number of chronolgy controls for the length of the record

min.chron.control<-function(tax_dec_dl,max.chron.grain){
	well.dated.sites<-c()
	tsuga.chrons<-c()
	for (i in 1:length(tax_dec_dl)){
		print(i)
		sitei<-tax_dec_dl[[i]]
		
		sitei.id<-sitei$dataset$dataset.meta$dataset.id
		controls<-get_chroncontrol(sitei) #this is also slow, but unavoidable (?)
		print(controls)
		
		#get duration of the time series
		duration<-max(sitei$sample.meta$age,na.rm=T)-min(sitei$sample.meta$age,na.rm=T)
		
		# 9Jan18 Note to Alistair and Trisha: these are the "bad" chron control types I found in the hemlock dataset--we will likely need to search for and idenitfy others
		bad.chrons<-c("Tsuga decline","Biostratigraphic, pollen","Sediment stratigraphic","Guess","Interpolated")
		controls.temp<-cbind(as.character(controls[[1]]$control.type),controls[[1]]$age)
		
		if (length(which(as.numeric(controls.temp[,2])>200))==1) {
			controls.temp2<-controls.temp[which(as.numeric((controls.temp[,2]))>200)]
			xx<-(controls.temp2 %in% bad.chrons)
		} else if (length(which(as.numeric(controls.temp[,2])>200))==0){
			xx<-FALSE
		} else {
			controls.temp2<-controls.temp[which(as.numeric((controls.temp[,2]))>200),]
			xx<-(controls.temp2[,1] %in% bad.chrons)
		}
		
		if (length(which(xx==TRUE))==0) {
			nchrons<-length(controls[[1]]$control.type)
		} else {
			nchrons<-length(controls[[1]]$control.type)-length(which(xx==TRUE))
		} 
		
		if (nchrons==1) { #if there is only 1 chron control, toss
			next
		} else if (nchrons==0) {
			next
		} else if (duration/nchrons>max.chron.grain) { #toss sites with a max grain of more than 1 chron control per XXX years (2000 is reasonable)
			next
		} else { #if N chron controls is adequate, send ids into a vector
			well.dated.sites<-c(well.dated.sites,sitei$dataset$dataset.meta$dataset.id)			
		}
	}
	return(well.dated.sites)
}

	
#min_pol_pct________________________________________________________
# function to select ids for datasets where a single taxon pollen reaches a minimum %
# tax_dec_dl is a Neotoma download object
# eco.group is a vector of Neotoma database ecological group codes.
# taxon is a species or group name. * indicates partial matching
# min.pct is the pollen % cut off; if the site never acheives min.pct of taxon pollen, it is excluded
# pct.zeros is the maximum allowable number of pollen samples where taxon is unsampled
# if eco.sort="pulished" function chooses a published taxon list to use for calculating pollen %, else eco.group is used

min_pol_pct<-function(tax_dec_dl,eco.group,taxon,min.pct,pct.zeros,eco.sort){
	dataset.list<-c()
	for (i in 1:length(tax_dec_dl)){
		print(i)
		sitei<-tax_dec_dl[[i]]
		dataset.list[i]<-sitei$dataset$dataset.meta$dataset.id
	}
	high.pct<-c()
	max.pct<-c()
	Hpct.ids<-c()
	Hpct.names<-c()
	for (k in 1:length(dataset.list)){
		sitei<-tax_dec_dl[[k]]
		sitei.counts<-sitei$counts
		all_taxa <- do.call(rbind.data.frame, lapply(tax_dec_dl, function(x)x$taxon.list[,1:6]))
		all_taxa <- all_taxa[!duplicated(all_taxa),]
	if (eco.sort=="published"){
		counts.subset<-compile_taxa(sitei,list.name="WS64")
		sitei.counts<-counts.subset$counts[,-grep("Lycopodium", colnames(sitei.counts))]
	} else if (eco.sort=="ecological.group") {
		good_cols<-c(which(colnames(sitei.counts) %in% all_taxa[all_taxa$ecological.group %in% eco.group,1]))	
		sitei.counts<-sitei.counts[,good_cols]
	} else {
		print("choose eco.sort = 'published' or 'ecological group'")
	}

		tax_pct<-sitei.counts[,1:ncol(sitei.counts)]/rowSums(sitei.counts[,1:ncol(sitei.counts)],na.rm=TRUE)
		sitei.taxon<-tax_pct[,grep(taxon, colnames(tax_pct))]
		max.pct[k]<-max(sitei.taxon[])
	
		if (max(sitei.taxon,na.rm=TRUE)<min.pct){
			next
		} else if ((sum(sitei.taxon==0)/length(sitei.taxon))>pct.zeros) {
			next
		} else {
			high.pct<-c(high.pct,k)
		}
	}
	Hpct.ids<-dataset.list[high.pct]
}


#trim_to________________________________________________________
# function to trim dataset to a minimum age and determine if there are enough remaining datapoints
# tax_dec_dl is a Neotoma download object
# min.age is the age cut off: datapoints younger than min.age will by trimmed
# cutoff is the minimum allowable number of remaining points 

trim_to<-function(tax_dec_dl,min.age,cutoff){
	ind<-c()
	for (i in 1:length(tax_dec_dl)) {
		sitei<-tax_dec_dl[[i]]
		age<-sitei$chronologies[[1]]$age
		length(age[which(age>min.age)])
		if (length(age[which(age>min.age)])<cutoff){ 
			ind[i]<-NA
		} else {
			ind[i]<-sitei$dataset$dataset.meta$dataset.id
		}
	}
		return(adequate.n=ind)
}



#map_dl________________________________________________________
# function to map location of sites
# tax_dec_dl is a Neotoma download object
# X is a vector with 2 elements: min and max longitude
# Y is a vector with 2 elements: min and max latitude
# add: should points be added to an existing map?
# color: color to use for points
# label.sites: should site names be added as text to the map

map_dl<-function(tax_dec_dl,X,Y,add,color,label.sites){
	if (add==FALSE){
		map("world",xlim=X,ylim=Y)
	} 
	lat<-c()
	long<-c()
	dataset.id<-c()
	site.name<-c()
	site.id<-c()
	for (i in 1:length(tax_dec_dl)){
		long[i]<-tax_dec_dl[[i]]$dataset$site.data$long
		lat[i]<-tax_dec_dl[[i]]$dataset$site.data$lat
		points(long[i],lat[i],pch=16,col=color)
		dataset.id[i]<-tax_dec_dl[[i]]$dataset$dataset.meta$dataset.id
		site.name[i]<-tax_dec_dl[[i]]$dataset$site.data$site.name
		site.id[i]<-tax_dec_dl[[i]]$dataset$site.data$site.id
	}
	if (label.sites==TRUE){
		text(long,lat,site.id,cex=0.5,pos=4,offset=0.2)
	}
	
	out<-cbind(site.name,dataset.id,site.id,lat,long)
	return(out)
}

#END FUNCTIONS________________________________________


library(neotoma)
#library(dplyr)
library(maps)

# query Neotoma database
# requires internet connection

#pull pollen datasets for North America______
tax_dec<-get_dataset(
	datasettype = "pollen",
	loc = c(-180, 10, -50, 89.9999)) #north america

summary(tax_dec)

#runs slowly!
pollen_dl<-get_download(tax_dec)	
pol_dl<-pollen_dl

# OR pull diatom data for North America______
#tax_dec<-get_dataset(datasettype = "pollen") #global
tax_dec<-get_dataset(
	datasettype = "diatom",
	loc = c(-150, 20, -60, 89.9999)) 
	
summary(tax_dec)
	
#runs slowly!
dia_dl<-get_download(tax_dec)	
pol_dl<-dia_dl


#limit to sites with chronologies
pol_chron<-has.chron(hem_dec_dl) #expect warnigns here. Not an issue
pol_chron<-pol_chron[complete.cases(pol_chron)]
pol_dl_sub1<-hem_dec_dl[as.character(pol_chron)]

#limit to sites with at least n samples
pol_n<-ts.min.length(pol_dl_sub1,20)
pol_n<-pol_n[complete.cases(pol_n)]
pol_dl_sub2<-pol_dl_sub1[as.character(pol_n)]

#limit to sites with minimum level of average resolution
pol_high<-select.high.res(pol_dl_sub2,200)
pol_dl_sub3<-pol_dl_sub2[as.character(pol_high)]
length(pol_dl_sub3)

#limit to sites with adequate chron controls
pol_chroncont<-min.chron.control(pol_dl_sub3,2000)
pol_dl_sub4<-pol_dl_sub3[as.character(pol_chroncont)]
length(pol_dl_sub4)


# The code below was specific to our hemlock project, but I'm leaving it in just in case it's useful for other purposes
# limit to sites that reach minimum % Tsuga pollen
# and for which Tusga pollen is sampled in at least 50% of samples
eco<-c("TRSH", "UPHE")
taxon<-"Tsuga*"
min.pct<-0.1
pct.zeros<-0.5

min.tsuga<-min_pol_pct(pol_dl_sub4,eco,"Tsuga*",0.1,0.5,"ecological.group")	
pol_dl_sub5<-pol_dl_sub4[as.character(min.tsuga)]
length(pol_dl_sub5)

# limit to sites with at least 20 datapoints older than 1000 years BP
ids<-trim_to(pol_dl_sub5,1000,20)
id.list<-ids[complete.cases(ids)]
pol_dl_sub6<-pol_dl_sub5[as.character(id.list)]

# map sites
X<-c(-180,-50)
Y<-c(10,90)
mapX<-map_dl(pol_dl_sub6,X,Y,add=FALSE,color="blue",label.sites=FALSE)