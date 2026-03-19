## This script is for filtering and plotting bcftools roh output
## Input is (optionally gzipped) .out files from bcftools ROH
## subsetted down to have one file for each individual
## as the total file size can be quite large

# Author: Chris Kyriazis


#### Set working directory, file names, individuals, etc ####

setwd("~/Documents/SDZWA/alala_genomics/analysis/ROH/data/175AL_round1/")

# set filename excluding individual and ".out.gz"
file <- "_175AL_FilterB_Round1_autosomes_PASS"

## set array of individuals 

individuals <- c("AL101","AL102","AL103","AL114","AL119","AL120","AL121","AL124","AL125","AL127",
                 "AL129","AL133","AL134","AL137","AL138","AL141","AL142","AL143","AL146","AL149",
                 "AL153","AL155","AL156","AL159","AL160","AL162","AL166","AL170","AL173","AL174",
                 "AL175","AL176","AL177","AL178","AL180","AL181","AL182","AL184","AL185","AL187",
                 "AL189","AL193","AL196","AL197","AL199","AL200","AL201","AL204","AL207","AL208",
                 "AL209","AL210","AL213","AL215","AL216","AL220","AL223","AL231","AL238","AL241",
                 "AL242","AL243","AL244","AL246","AL254","AL255","AL256","AL263","AL266","AL269",
                 "AL272","AL273","AL274","AL276","AL278","AL279","AL283","AL284","AL288","AL292",
                 "AL300","AL303","AL304","AL307","AL309","AL313","AL32","AL35","AL75","AL82","AL86",
                 "AL88","AL91","AL92","AL95","AL97","AL99","K1423","K1461","K1542","K1544","K1567",
                 "K1579","K1580","K1583","K1603","K1606","K1608","K16100","K1614","K1617","K1618",
                 "K1623","K1627","K1632","K1640","K1650","K1660","K1663","K1665","K1669","K1672",
                 "K1680","K1684","K1693","K1696","K17007","K17008","K17010","K17013","K17014","K17025",
                 "K17028","K17030","K17032","K17033","K17034","K17035","K17037","K17040","K17042","K17044",
                 "K17045","K17048","K17052","K17053","K17056","K17062","K17064","K17073","K17077","K17078",
                 "K17081","K17083","K17087","K17091","K17094","K17098","K17099","K17101","K17102","K17103","K17104",
                 "K17106","K18010","K18017","K18019","K18022","K18034","K18078","K18083","K18084","K18091",
                 "K19018","K19019")

chroms <- c(99062180,123451405,121534940,80112800,76278832,66143299,40978052,38659334,33574751,28816175,
            23986915,23235100,22362767,23265108,21088201,17487702,16938386,15572589,12787533,13899114,
            13093161,8696975,7971722,11344391,7312786,46055897,6387999,7803721,3496993,22021689,4714721,
            736907,477022,104925,1321807,1784968,336161,251864,167149,337809,123330)
genome_length <- sum(chroms)/1e6

# min size allowable for ROHs - typically 100kb or 300kb
min_roh_length=100000


#### Define functions ####

## Define function that takes in data frame and divides it into three length classes
classify_roh <- function(roh_dataframe, min_roh_length){
  short_roh <- subset(roh_dataframe,length>min_roh_length & length<1000000) # roh_dataframe[length>100000 & length<1000000]
  med_roh <- subset(roh_dataframe, length>1000000 & length<10000000)
  long_roh <-  subset(roh_dataframe, length>10000000 & length<100000000)
  
  #sum each class and divide by 1000000 to convert to Mb
  sum_short_Mb <- sum(short_roh$length)/1000000
  sum_med_Mb <- sum(med_roh$length)/1000000
  sum_long_Mb <- sum(long_roh$length)/1000000
  
  print(paste("This individual has",dim(short_roh)[1],"short ROHs summing to",sum_short_Mb, "Mb",
              dim(med_roh)[1],"medium ROHs summing to",sum_med_Mb,"Mb, and",
              dim(long_roh)[1],"long ROHs summing to",sum_long_Mb, "Mb"))
  
  return(c(sum_short_Mb, sum_med_Mb, sum_long_Mb))
  #roh_matrix_Mb <- roh_matrix/1000000 #convert to Mb
  
}

## Define function to read in output files for each individual and filter out ROHs less than min_roh_length
read_filter_roh <- function(data, min_roh_length){
  output <- read.table(paste(data,".out.gz",sep=""), col.names=c("row_type","sample","chrom","start","end","length","num_markers","qual"), fill=T)
  output1 <- subset(output, row_type == "RG")
  output_class_sums <- classify_roh(output1, min_roh_length=min_roh_length)
  return(output_class_sums)
}


#### Read in data, filter, and classify ROHs ####

## initialize data frame
roh_size_df <- data.frame(matrix(nrow=3, ncol=length(individuals)))
colnames(roh_size_df) <- individuals
froh_1mb <- c()
froh_10mb <- c()

## read in data for each individual
## note that this can be VERY slow - takes several min per individual
roh_size_df <- read.csv(file = "roh_size_df.csv", header=T)

## read in data for each individual
## note that this can be VERY slow - takes several min per individual
for(i in 1:length(individuals)){
  #roh_size_df[,i] <- read_filter_roh(paste(individuals[i],file, sep=""), min_roh_length=min_roh_length)
  froh_1mb = c(froh_1mb,sum(roh_size_df[2:3,i])/genome_length) # sum ROHs > 1Mb and divide by genome length to estimate Froh
  froh_10mb = c(froh_10mb,sum(roh_size_df[3,i])/genome_length) # sum ROHs > 5Mb and divide by genome length to estimate Froh
}

#write.csv("roh_size_df.csv", x = roh_size_df, row.names=F)




#### Plot results ####

# stacked barplots of binned ROHs
setwd("~/Documents/SDZWA/alala_genomics/analysis/ROH/plots/")
pdf(paste("barplot",file,".pdf",sep=""), width=24, height=8)

par(mar=c(5,5,2,1))

barplot(as.matrix(roh_size_df), names.arg = individuals, col=c("seashell","seashell3", "seashell4"), ylab="Summed ROH length (Mb)", xlab="", ylim=c(0,650), cex.names=0.75, las=2, cex.lab = 1.5, cex.axis=1.2)
legend("topleft",legend=c(paste(min_roh_length/1000000,"-1 Mb",sep=""), "1-10 Mb", ">10Mb"), col=c("red","orange", "yellow"), fill=c("seashell","seashell3", "seashell4"), cex=1.3)

dev.off()




pdf("~/Documents/SDZWA/alala_genomics/analysis/ROH/Plots/froh_vs_depth_alala.pdf", width=6, height=5)
par(mar=c(5,5,2,2))
#depth from scaff1
depth <- c(18.942,21.1717,20.7955,22.0912,14.3123,17.444,15.0665,16.1825,24.7138,24.7377,22.0079,18.7245,1.03158,3.09115,23.661,24.0613,13.6141,16.2613,23.6201,3.1785,2.84681,25.0609,2.57365,15.8341,29.1006,16.7831,14.3039,18.6119,21.5584,18.928,15.0291,22.899,17.7855,2.62858,15.8183,1.52258,16.1097,16.5356,16.9394,21.6746,14.7755,15.6524,14.6574,18.2099,15.2478,16.1442,14.5662,21.3868,14.1061,2.86584,15.1224,23.8733,16.5558,21.8666,16.2248,21.4093,6.85374,17.3912,22.1901,3.48188,3.77107,17.0246,13.6137,15.0737,20.6925,2.45495,3.44111,28.4479,17.867,17.6361,22.4343,2.7875,18.0402,17.8618,3.1269,36.4858,27.302,27.135,21.5568,16.4171,14.4264,4.81038,10.5208,17.1831,16.8378,24.7924,19.5613,13.7649,18.6923,19.0218,23.8702,19.0251,19.1585,24.5887,20.2818,19.1839,19.6242,13.5889,19.4837,17.3535,13.6242,23.125,43.876,15.3453,18.4149,16.2226,26.3075,24.5687,19.118,16.2292,24.5577,22.015,3.66354,5.88908,18.7709,26.1311,18.1388,10.3518,12.3867,20.0498,29.4444,30.6787,18.2532,29.142,14.3787,17.0793,16.595,2.67335,24.8883,9.85173,15.6857,18.5882,16.6608,16.3087,41.3353,32.0196,15.0533,12.7469,16.2075,22.0387,22.9181,23.1671,87.1928,15.4949,33.5358,17.8981,24.5529,20.2775,16.398,17.8178,17.6112,18.4784,21.0429,19.7268,25.5016,18.2274,22.1027,18.0724,17.5675,26.4918,20.6466,31.6978,17.9841,21.9524,15.3673,24.4391,22.0461,17.4972,11.0178,20.1232,17.4424,16.5904,16.3963,17.4015,15.8311)


plot(depth,froh_1mb,ylab = expression('F'[ROH]), xlab = "Sequencing depth", pch=19, cex=1.2,cex.lab=1.2)

loess_fit <- loess(froh_1mb ~ depth, span = 0.75)  # Adjust `span` for smoothness
x_new <- seq(min(depth), max(depth), length.out = 100)  # Fine grid for smooth curve
y_pred <- predict(loess_fit, newdata = x_new)

lines(x_new, y_pred, col = "gray", lwd = 3)


dev.off()

