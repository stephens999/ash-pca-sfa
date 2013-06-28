cmd_args = commandArgs()
#Data1 <- t(read.table(cmd_args[5], header=FALSE))
Data1 <- read.table(cmd_args[5], header=FALSE)
	
Data1 <- Data1[2:length(Data1[,1]),] #Getting rid of top header (first row)

Data1[,4:length(Data1[1,])] <-  apply(Data1[,4:length(Data1[1,])], 2, as.numeric)

#Data2 <- Data1[,4:length(Data1[1,])] #If want to check median of all gene expression data per individual -- commonly done as QC for other dataset checks
#medianRslts <- apply(Data2, 1, median)
#hist(medianRslts)

#for (i in 4:length(Data1[1,])) {		#Quantile-Normalizing each gene
#	Data1[Data1[,3]=="CEPH",i] <- qqnorm(Data1[Data1[,3]=="CEPH",i], plot.it=FALSE)$x
#	Data1[Data1[,3]=="YRI",i] <- qqnorm(Data1[Data1[,3]=="YRI",i], plot.it=FALSE)$x
#	Data1[Data1[,3]=="CHB",i] <- qqnorm(Data1[Data1[,3]=="CHB",i], plot.it=FALSE)$x
#}

#Within-population quantile-normalizing each gene
#Data1[Data1[,3]=="CEPH",4:length(Data1[1,])] <- apply(Data1[Data1[,3]=="CEPH",4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)
#Data1[Data1[,3]=="YRI",4:length(Data1[1,])] <- apply(Data1[Data1[,3]=="YRI",4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)
#Data1[Data1[,3]=="CHB",4:length(Data1[1,])] <- apply(Data1[Data1[,3]=="CHB",4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)

#Quantile-normalizing each gene independent of sub-populations
#Data1[,4:length(Data1[1,])] <- apply(Data1[,4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)

#Quantile-normalizing a random subset of 90, 60 and 60
Seq1 <- sample(1:210)
Data1[Seq1[1:90],4:length(Data1[1,])] <- apply(Data1[Seq1[1:90],4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)
Data1[Seq1[91:150],4:length(Data1[1,])] <- apply(Data1[Seq1[91:150],4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)
Data1[Seq1[151:210],4:length(Data1[1,])] <- apply(Data1[Seq1[151:210],4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)

Data1 <- Data1[,4:length(Data1[1,])] #Getting rid of row headers (first two columns)

Means <- apply(Data1, 2, mean)	#Getting Mean and SD for each column
SDs <- apply(Data1, 2, sd)

#I <- matrix(rep(1, length(Data1[,1])*length(Data1[1,])), nrow=length(Data1[,1]))
#MeanMtrx <- matrix(rep(Means, length(Data1[,1])), byrow=TRUE, nrow=length(Data1[,1]))	#Creating Mean and SD matrices
#SDMtrx <- matrix(rep(SDs, length(Data1[,1])), byrow=TRUE, nrow=length(Data1[,1]))

#Data1 <- (Data1 - MeanMtrx) / SDMtrx	#Centralizing and Standardizing main data matrix

CovMtrx <- (1/length(Data1[1,])) * (as.matrix(Data1) %*% t(as.matrix(Data1)))		#Creating 1/p * covariance matrix

Results1 <- eigen(CovMtrx)		#Getting eigenvectors and eigenvalues

write.table(Results1$vectors, cmd_args[6], quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(Results1$values, cmd_args[7], quote=FALSE, row.names=FALSE, col.names=FALSE)


##############
###Data log###
##############

#> Data1 <- read.table("/mnt/lustre/home/mturchin20/Data/Stranger07/beforeCorrection_onlyExpGenes.txt.pedOrder.Transposed.wHPMPid.gz", header=FALSE)
#> Data1[1:5,1:5]
#   V1      V2  V3              V4              V5
#   1 FID     IID POP ENSG00000187634 ENSG00000187961
#   2   1 NA18524 CHB         6.42588         9.00134
#   3   2 NA18526 CHB         6.34857         8.68707
#   4   3 NA18529 CHB         6.46448         8.76156
#   5   4 NA18532 CHB         6.37794         9.17925
#   > Data1 <- Data1[2:length(Data1[,1]),] #Getting rid of top header (first row)
#   >
#   > Data1[,4:length(Data1[1,])] <-  apply(Data1[,4:length(Data1[1,])], 2, as.numeric)
#   >
#   > Data1[1:5,1:5]
#     V1      V2  V3      V4      V5
#     2  1 NA18524 CHB 6.42588 9.00134
#     3  2 NA18526 CHB 6.34857 8.68707
#     4  3 NA18529 CHB 6.46448 8.76156
#     5  4 NA18532 CHB 6.37794 9.17925
#     6  5 NA18537 CHB 6.39174 8.70063
#     > Data1[,4:length(Data1[1,])] <- apply(Data1[,4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)
#     > Data1[1:5,1:5]
#       V1      V2  V3         V4         V5
#       2  1 NA18524 CHB -0.4770404  0.5039654
#       3  2 NA18526 CHB -1.0781824 -0.6156959
#       4  3 NA18529 CHB -0.1497620 -0.3470265
#       5  4 NA18532 CHB -0.8501565  1.0570776
#       6  5 NA18537 CHB -0.7201566 -0.5589579
#       > Data1 <- Data1[,4:length(Data1[1,])] #Getting rid of row headers (first two columns)
#       >
#       > Means <- apply(Data1, 2, mean)  #Getting Mean and SD for each column
#       > SDs <- apply(Data1, 2, sd)
#       >
#       > #I <- matrix(rep(1, length(Data1[,1])*length(Data1[1,])), nrow=length(Data1[,1]))
#       > MeanMtrx <- matrix(rep(Means, length(Data1[,1])), byrow=TRUE, nrow=length(Data1[,1]))   #Creating Mean and SD matrices
#       > SDMtrx <- matrix(rep(SDs, length(Data1[,1])), byrow=TRUE, nrow=length(Data1[,1]))
#       >
#       > #Data1 <- (Data1 - MeanMtrx) / SDMtrx   #Centralizing and Standardizing main data matrix
#       >
#       > CovMtrx <- (1/length(Data1[1,])) * (as.matrix(Data1) %*% t(as.matrix(Data1)))           #Creating 1/p * covariance matrix
#       >
#       > Results1 <- eigen(CovMtrx)              #Getting eigenvectors and eigenvalues
#       > CovMtrx[1:5,1:5]
#                  2          3           4           5          6
#		  2 0.95511276  0.0540701  0.06472216  0.47045660 0.08949228
#		  3 0.05407010  1.3437347  0.60808475 -0.05819240 0.55403112
#		  4 0.06472216  0.6080847  1.26739125 -0.06139981 0.36951139
#		  5 0.47045660 -0.0581924 -0.06139981  0.79545792 0.01013631
#		  6 0.08949228  0.5540311  0.36951139  0.01013631 0.76484822
#		  > beforeCorrection_onlyExpGenes.txt.pedOrder.Transposed.wthnPopNorm
#		  > write.table(Results1$vectors, "beforeCorrection_onlyExpGenes.txt.pedOrder.Transposed.allPopNorm.eigenVectors.txt" , quote=FALSE, row.names=FALSE, col.names=FALSE)
#		  > write.table(Results1$values, "beforeCorrection_onlyExpGenes.txt.pedOrder.Transposed.allPopNorm.eigenValues.txt" , quote=FALSE, row.names=FALSE, col.names=FALSE)
#		  > write.table(CovMtrx, "beforeCorrection_onlyExpGenes.txt.pedOrder.Transposed.allPopNorm" , quote=FALSE, row.names=FALSE, col.names=FALSE)

#20130626
#> Data1 <- read.table("beforeCorrection_onlyExpGenes.txt.pedOrder.Transposed.wHPMPid.gz", header=FALSE)
#> Data1 <- Data1[2:length(Data1[,1]),] #Getting rid of top header (first row)
#>
#> Data1[,4:length(Data1[1,])] <-  apply(Data1[,4:length(Data1[1,])], 2, as.numeric)
#>
#> #Data2 <- Data1[,4:length(Data1[1,])] #If want to check median of all gene expression data per individual -- commonly done as QC for other dataset checks
#> #medianRslts <- apply(Data2, 1, median)
#> #hist(medianRslts)
#>
#> #for (i in 4:length(Data1[1,])) {               #Quantile-Normalizing each gene
#> #       Data1[Data1[,3]=="CEPH",i] <- qqnorm(Data1[Data1[,3]=="CEPH",i], plot.it=FALSE)$x
#> #       Data1[Data1[,3]=="YRI",i] <- qqnorm(Data1[Data1[,3]=="YRI",i], plot.it=FALSE)$x
#> #       Data1[Data1[,3]=="CHB",i] <- qqnorm(Data1[Data1[,3]=="CHB",i], plot.it=FALSE)$x
#> #}
#>
#> #Within-population quantile-normalizing each gene
#> #Data1[Data1[,3]=="CEPH",4:length(Data1[1,])] <- apply(Data1[Data1[,3]=="CEPH",4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)
#> #Data1[Data1[,3]=="YRI",4:length(Data1[1,])] <- apply(Data1[Data1[,3]=="YRI",4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)
#> #Data1[Data1[,3]=="CHB",4:length(Data1[1,])] <- apply(Data1[Data1[,3]=="CHB",4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)
#>
#> #Quantile-normalizing each gene independent of sub-populations
#> #Data1[,4:length(Data1[1,])] <- apply(Data1[,4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)
#>
#> #Quantile-normalizing a random subset of 90, 60 and 60
#> Seq1 <- sample(1:210)
#> Seq1
#  [1] 184  14 148 107   5 152 189 169   1 166 192  28 139 146  47  19  62  36
#   [19] 150  83  45 125 163  85  11 138 122  82 183 151  79  38 190 195 165  17
#       > Data1[Seq1[1:5],1:5]
#            V1      V2   V3      V4      V5
#	    185 184 NA19131  YRI 6.46336 8.75060
#	    15   14 NA18561  CHB 6.30135 8.96216
#	    149 148 NA12875 CEPH 6.72132 9.10995
#	    108 107 NA11881 CEPH 6.51859 9.06993
#	    6     5 NA18537  CHB 6.39174 8.70063
#	    > 
#> Data1[Seq1[1:90],4:length(Data1[1,])] <- apply(Data1[Seq1[1:90],4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)
#> Data1[Seq1[1:5],1:5]
#     V1      V2   V3         V4         V5
#     185 184 NA19131  YRI -0.1820346 -0.3853205
#     15   14 NA18561  CHB -1.9145058  0.2967378
#     149 148 NA12875 CEPH  2.5391848  0.8616341
#     108 107 NA11881 CEPH  0.2104284  0.7461852
#     6     5 NA18537  CHB -0.8219399 -0.5729675
#     > Data1[1:20,1:5]
#        V1      V2  V3         V4         V5
#	2   1 NA18524 CHB -0.5404469  0.5084881
#	3   2 NA18526 CHB  6.3485700  8.6870700
#	4   3 NA18529 CHB  6.4644800  8.7615600
#	5   4 NA18532 CHB  6.3779400  9.1792500
#> Data1[Seq1[101:105],1:5]
#     V1      V2   V3      V4      V5
#     148 147 NA12874 CEPH 6.88838 9.59689
#     106 105 NA11839 CEPH 6.71669 9.10257
#     23   22 NA18573  CHB 6.38641 8.95633
#     8     7 NA18542  CHB 6.31718 8.65836
#     3     2 NA18526  CHB 6.34857 8.68707
#     
#> Data1[Seq1[91:150],4:length(Data1[1,])] <- apply(Data1[Seq1[91:150],4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)
#> Data1[Seq1[151:210],4:length(Data1[1,])] <- apply(Data1[Seq1[151:210],4:length(Data1[1,])], 2, function(y) qqnorm(y, plot.it=FALSE)$x)
#> Data1[1:20,1:5]
#   V1      V2  V3          V4          V5
#   2   1 NA18524 CHB -0.54044687  0.50848806
#   3   2 NA18526 CHB -1.23544034 -0.75541503
#   4   3 NA18529 CHB -0.14674496 -0.23183436
#   5   4 NA18532 CHB -0.81221780  1.15034938
#   6   5 NA18537 CHB -0.82193989 -0.57296755
#   7   6 NA18540 CHB -1.43953147 -0.14674496
#   8   7 NA18542 CHB -1.73166440 -0.81221780
#
#write.table(Data4b, "beforeCorrection_onlyExpGenes.txt.pedOrder.Transposed.RanSubStrct2", quote=FALSE, row.names=FALSE, col.names=FALSE) #Moved the first two columns of Data1 onto Data4b through multiple steps and later manually added the header column outside of R, hence the use of Data4b out of nowhere here
#> Data1[1:5,1:5]
#         V4         V5         V6         V7         V8
#	  2 -0.5404469  0.5084881  0.9899016 -1.3829941  2.5391848
#	  3 -1.2354403 -0.7554150  1.0013313  1.5689196  1.4395315
#	  4 -0.1467450 -0.2318344 -0.7554150 -1.7316644 -0.5977601
#	  5 -0.8122178  1.1503494 -0.4537622 -0.4079187  1.0013313
#	  6 -0.8219399 -0.5729675 -0.8616341 -0.7461852  0.1256613
#> CovMtrx <- (1/length(Data1[1,])) * (as.matrix(Data1) %*% t(as.matrix(Data1)))           #Creating 1/p * covariance matrix
#>
#> Results1 <- eigen(CovMtrx)              #Getting eigenvectors and eigenvalues
#> write.table(CovMtrx, "beforeCorrection_onlyExpGenes.txt.pedOrder.Transposed.RanSubStrct2.CovMtrx", quote=FALSE, row.names=FALSE, col.names=FALSE)
#> beforeCorrection_onlyExpGenes.txt.pedOrder.Transposed.RanSubStrct2
#> write.table(Results1$vectors, "beforeCorrection_onlyExpGenes.txt.pedOrder.Transposed.RanSubStrct2.eigenVectors.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
#> write.table(Results1$values, "beforeCorrection_onlyExpGenes.txt.pedOrder.Transposed.RanSubStrct2.eigenValues.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)



