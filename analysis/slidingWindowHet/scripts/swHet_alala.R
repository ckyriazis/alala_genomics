# Script for plotting heterozygosity in sliding windows across the genome
# Uses text files containing the following columns (for example):
# chromo	window_start		sites_total	sites_unmasked	sites_passing	sites_variant	calls_sample1	calls_sample2	calls_sample3	hets_sample1	hets_sample2		hets_sample3

library(plyr)
library(gtools)

file <- "174AL_FilterB_Round2"
window_size=1000000

setwd(paste("~/Documents/SDZWA/alala_genomics/analysis/slidingWindowHet/data/", file, sep=""))

# All files have *step.txt naming convention. Get list of files:
winfiles=list.files(pattern="step.txt")
winfiles <- mixedsort(winfiles)
winfiles <- winfiles[1:20]

# This is here to reorder the chromosomes and exclude the X
nchr=length(winfiles)

# Read in chromosome 1
hetdata=read.table(winfiles[1], header=T, sep='\t')

# Add subsequent chromosomes
for (i in 2:nchr){ 
	temp=read.table(winfiles[i], header=T, sep='\t')
	hetdata=rbind(hetdata,temp)
}

#filter out windows with fewer than x% of sites called
#NOTE that the hetdata dataframe needs to be reset if changing this
threshold=0.2
hetdata= subset(hetdata, hetdata$sites_total>window_size*threshold) 
dim(hetdata)
hist(hetdata$sites_total)
rownames(hetdata) <- 1:nrow(hetdata) # need to reset rownames since these are used to plot



# Get chromosome position info needed for plotting later:
# Get the start positions of each chromosome
pos=as.numeric(rownames(unique(data.frame(hetdata $chrom)[1])))

# Add the end position
pos=append(pos,length(hetdata $chrom))

# Get the midpoints to center chromosome labels on x-axis
numpos=NULL
for (i in 1:length(pos)-1){numpos[i]=(pos[i]+pos[i+1])/2}


# Move out of data files folder
setwd("~/Documents/SDZWA/alala_genomics/analysis/slidingWindowHet/plots")

# Get the number of samples. Here, the number of columns is 6 + 2*number of samples. Change if necessary.
numsamps=174

# Which columns contain the numbers of calls and the numbers of hets?
callcols=seq(4,3+numsamps)
hetcols=seq(4+numsamps,length(names(hetdata)))



# Get the sample names
sampnames=gsub("calls_", "", names(hetdata)[callcols])

# Which individuals to plot?
samp1='AL101'
samp2='AL102'
samp3='AL103'
samp4='AL114'
samp5='AL119'
samp6='AL120'
samp7='AL121'
samp8='AL124'
samp9='AL125'
samp10='AL127'
samp11='AL129'
samp12='AL133'
samp13='AL134'
samp14='AL137'
samp15='AL138'
samp16='AL141'
samp17='AL142'
samp18='AL143'
samp19='AL146'
samp20='AL149'
samp21='AL153'
samp22='AL155'
samp23='AL156'
samp24='AL159'
samp25='AL160'
samp26='AL162'
samp27='AL166'
samp28='AL170'
samp29='AL173'
samp30='AL174'
samp31='AL175'
samp32='AL176'
samp33='AL177'
samp34='AL178'
samp35='AL180'
samp36='AL181'
samp37='AL182'
samp38='AL184'
samp39='AL185'
samp40='AL187'
samp41='AL189'
samp42='AL193'
samp43='AL196'
samp44='AL197'
samp45='AL199'
samp46='AL200'
samp47='AL201'
samp48='AL204'
samp49='AL207'
samp50='AL208'
samp51='AL209'
samp52='AL210'
samp53='AL213'
samp54='AL215'
samp55='AL216'
samp56='AL220'
samp57='AL223'
samp58='AL231'
samp59='AL238'
samp60='AL241'
samp61='AL242'
samp62='AL243'
samp63='AL244'
samp64='AL246'
samp65='AL254'
samp66='AL255'
samp67='AL256'
samp68='AL263'
samp69='AL266'
samp70='AL269'
samp71='AL272'
samp72='AL273'
samp73='AL274'
samp74='AL276'
samp75='AL278'
samp76='AL279'
samp77='AL283'
samp78='AL284'
samp79='AL288'
samp80='AL292'
samp81='AL300'
samp82='AL303'
samp83='AL304'
samp84='AL307'
samp85='AL309'
samp86='AL313'
samp87='AL32'
samp88='AL35'
samp89='AL75'
samp90='AL82'
samp91='AL86'
samp92='AL88'
samp93='AL91'
samp94='AL92'
samp95='AL95'
samp96='AL97'
samp97='AL99'
samp98='K1423'
samp99='K1461'
samp100='K1542'
samp101='K1544'
samp102='K1567'
samp103='K1579'
samp104='K1580'
samp105='K1583'
samp106='K1603'
samp107='K1606'
samp108='K1608'
samp109='K16100'
samp110='K1614'
samp111='K1617'
samp112='K1618'
samp113='K1623'
samp114='K1627'
samp115='K1632'
samp116='K1640'
samp117='K1650'
samp118='K1660'
samp119='K1663'
samp120='K1665'
samp121='K1669'
samp122='K1672'
samp123='K1680'
samp124='K1684'
samp125='K1693'
samp126='K1696'
samp127='K17007'
samp128='K17008'
samp129='K17010'
samp130='K17013'
samp131='K17014'
samp132='K17025'
samp133='K17028'
samp134='K17030'
samp135='K17032'
samp136='K17033'
samp137='K17034'
samp138='K17035'
samp139='K17037'
samp140='K17040'
samp141='K17042'
samp142='K17044'
samp143='K17045'
samp144='K17048'
samp145='K17052'
samp146='K17053'
samp147='K17056'
samp148='K17062'
samp149='K17064'
samp150='K17073'
samp151='K17077'
samp152='K17078'
samp153='K17081'
samp154='K17083'
samp155='K17087'
samp156='K17091'
samp157='K17094'
samp158='K17098'
samp159='K17099'
samp160='K17101'
samp161='K17102'
samp162='K17103'
samp163='K17106'
samp164='K18010'
samp165='K18017'
samp166='K18019'
samp167='K18022'
samp168='K18034'
samp169='K18078'
samp170='K18083'
samp171='K18084'
samp172='K18091'
samp173='K19018'
samp174='K19019'


# Use alternating colors to distinguish chromosomes
mycols=rep(c("chartreuse4","chartreuse4"), nchr/2)


hetplot=function(sampname, ymax, title, xlab, ylab){
	# Start an empty plot, change ylim as needed
	plot(0,0, type="n", xlim=c(0,pos[length(pos)]), ylim=c(0,ymax), axes=F, xlab="", ylab=ylab, main=title, cex.lab=1.2)
	
	# Add lines for chromosome data
	aa=which(sampnames==sampname)
	for (i in 1:nchr){
	  temp=hetdata[which(hetdata$chrom==unique(hetdata$chrom)[i]),]
	  x <- as.numeric(rownames(temp))
	  y <- temp[,hetcols[aa]]/temp[,callcols[aa]]

	  for(j in 1:length(x)){
	    lines(x=c(x[j],x[j]),c(0,y[j]), col=mycols[i], lwd=1.1)
	  }
	  
		#lines(as.numeric(rownames(temp)),temp[,hetcols[aa]]/temp[,callcols[aa]], col=mycols[i], lwd=1.5)

	}
	
	mean_het=mean(hetdata[,hetcols[aa]]/hetdata[,callcols[aa]], na.rm=T)
	print(mean_het)

	# Add y-axis
	axis(2)
	
	# Add x-axis and labels
	title(xlab=xlab, line=2, cex.lab=1.2)
	axis(side=1, at=pos, labels=F)
	axis(side=1, at=numpos, tick=F, labels=1:nchr, las=3, cex.axis=.8, line=-.2)
	return(mean_het)
}



# Save as file
pdf(paste(file,"_6samp.pdf",sep=""), width=8, height=4)

par(mfrow=c(3,2))
par(mar=c(4,4,2,1))

ymax=0.0015

# Define a list to store all the hetplot objects
het_plots <- list()

samples <- c(1,12,22,43,51,62)

# Loop to create hetplot objects
for (i in samples) {
  het_plots[[i]] <- hetplot(get(paste0("samp", i)), ymax, get(paste0("samp", i)), "Chromosome", "Heterozygosity")
}

# Naming the hetplot objects
names(het_plots) <- paste0("het_", samples)


# Close figure file
dev.off()






# Save as file
pdf(paste(file,"1.pdf",sep=""), width=12, height=24)

par(mfrow=c(20,3))
par(mar=c(4,4,2,1))

ymax=0.0015

# Define a list to store all the hetplot objects
het_plots1 <- list()

# Loop to create hetplot objects
for (i in 1:58) {
  het_plots1[[i]] <- hetplot(get(paste0("samp", i)), ymax, get(paste0("samp", i)), "Chromosome", "Heterozygosity")
}


# Close figure file
dev.off()


# Save as file
pdf(paste(file,"2.pdf",sep=""), width=12, height=24)

par(mfrow=c(20,3))
par(mar=c(4,4,2,1))

ymax=0.0015

# Define a list to store all the hetplot objects
het_plots2 <- list()

# Loop to create hetplot objects
for (i in 59:116) {
  het_plots2[[i]] <- hetplot(get(paste0("samp", i)), ymax, get(paste0("samp", i)), "Chromosome", "Heterozygosity")
}


# Close figure file
dev.off()



# Save as file
pdf(paste(file,"3.pdf",sep=""), width=12, height=24)

par(mfrow=c(20,3))
par(mar=c(4,4,2,1))

ymax=0.0015

# Define a list to store all the hetplot objects
het_plots3 <- list()

# Loop to create hetplot objects
for (i in 117:174) {
  het_plots3[[i]] <- hetplot(get(paste0("samp", i)), ymax, get(paste0("samp", i)), "Chromosome", "Heterozygosity")
}


# Close figure file
dev.off()




hets <- c(unlist(het_plots1), unlist(het_plots2), unlist(het_plots3))



### plot heterozygosity as a function of depth

pdf("alala_het_vs_depth.pdf" , width=6, height=5)
par(mar=c(5,5,2,2))

#depth from scaff1
depth <- c(18.2818,20.1779,19.7924,19.5881,13.8308,16.7056,14.5598,15.029,22.7435,21.5615,18.4655,16.8813,1.19248,2.76566,22.619,22.0028,13.6,15.6239,22.25,2.45287,2.92935,22.2958,2.05529,15.2898,27.8705,16.2215,13.7319,17.9327,19.2593,15.9075,14.6874,15.8819,15.9616,2.67406,15.2413,1.69128,15.4116,15.9558,16.3054,20.802,13.5612,15.0515,14.0095,17.4786,14.6281,15.5985,14.0914,18.1493,13.542,2.22328,14.6104,23.1014,15.9803,24.133,15.5821,17.1385,6.07508,15.5318,21.3363,2.71768,3.05631,16.0355,11.3681,13.7374,19.4131,2.02307,2.55042,29.445,16.5605,16.2438,21.3065,2.06069,16.71,16.9267,2.48198,34.8651,26.6125,25.951,21.0187,13.6631,13.663,3.95814,7.62491,16.6926,15.8024,23.735,18.5481,13.1859,17.1903,18.1246,22.7832,18.2285,18.3344,23.3395,17.0893,18.3666,18.8308,13.3181,16.643,15.7638,13.4236,22.0929,44.0667,14.6946,17.6096,12.0834,25.2417,19.3781,14.2979,12.719,22.0412,21.4627,3.07478,4.51878,18.0157,20.3542,15.9794,7.29203,9.17548,16.0727,28.1728,29.6051,16.0897,27.456,13.8452,16.2672,15.7181,1.85204,23.5423,7.1122,14.6584,14.9694,15.7557,15.5676,44.4572,31.0654,14.1245,12.7141,15.2413,16.761,21.2331,18.6228,69.9429,14.5474,31.1972,16.7123,18.7359,19.2906,14.9094,17.0001,17.0375,17.3769,20.0386,19.0325,24.0227,17.3502,16.5588,15.1134,16.1928,21.65,16.7264,43.9003,16.8157,14.3134,20.2245,17.7819,15.718,8.27694,16.5864,14.4709,14.3899,13.1962,15.3586,14.5636)
plot(depth, hets, ylab = 'Heterozygosity', xlab = 'Sequencing depth',  pch=19, cex=1.2, cex.lab=1.2)
loess_fit <- loess(hets ~ depth, span = 0.75)  # Adjust `span` for smoothness
x_new <- seq(min(depth), max(depth), length.out = 100)  # Fine grid for smooth curve
y_pred <- predict(loess_fit, newdata = x_new)

lines(x_new, y_pred, col = "gray", lwd = 3)


dev.off()


boxplot(hets[1:97], hets[98:174])



hets_highcov <- hets[depth>20]
mean(hets_highcov)
mean(hets_highcov)/(1-0.32)




# non-ROH het
mean(hetdata$hets_AL86[hetdata$hets_AL86>100]/hetdata$calls_AL86[hetdata$hets_AL86>100])
mean(hetdata$hets_AL88[hetdata$hets_AL88>100]/hetdata$calls_AL88[hetdata$hets_AL88>100])
mean(hetdata$hets_K18078[hetdata$hets_K18078>100]/hetdata$calls_K18078[hetdata$hets_K18078>100])



data <- read.csv("~/Documents/SDZWA/alala_genomics/analysis/slidingWindowHet/het_depth_alive.csv")

living_highcov <- data[data$depth.from.bam>20 & data$Alive==1,]
dim(living_highcov)

mean(living_highcov$het )

