
# get command line arguments.
comarg <- commandArgs()
input  <- comarg[6]
n      <- as.numeric(comarg[7])
m      <- as.numeric(comarg[8])
#outdir <- comarg[9]


n <- n*2+1
m<- m*2+1

# read count data per class (do not have to create matrix, as we simply multiply the vectors)
dat <- read.table(input, skip=1, h=F)


## create weights, i.e. (p1*q2)+(p2*q1) values for each count class

# allele frequency vectors for each group
P1 <- seq(0, 1, by=1/(n-1))
P2<- seq (0, 1, by=1/(m-1))
# create all pairwise combinations
# For comparisons in which not the same number of individuals are included the sequence in which these values are generated is crucial to match the SFS/counts contained in d.
# The angsd manual states: For 2dsfs the results is a single line, assuming we have n categories in population1 and m categories in population2, then the first m values will be the SFS for the first category in population1, etc.
weights <- vector(length=n*m)
c <- 1
for (i in 1:n){
  p1=P1[i]
  for (j in 1:m){
    p2=P2[j]
    weights[c] <- p1*(1-p2)+p2*(1-p1)
    c <- c+1
  }
}


## iterate through sfs per window to compute window-wise dxy

out.tab        <- data.frame(matrix(ncol=3,nrow=nrow(dat)))
names(out.tab) <-c( "Boot","dxy", "nsites")


##
for (i in 1:nrow(dat)){
  
  out.tab$Boot[i] <- as.character(i)
  d <- dat[i,1:ncol(dat)]
  out.tab$nsites[i] <- sum(d)
  out.tab$dxy[i] <- sum(d * weights)/out.tab$nsites[i]
  
}
out <- paste(gsub(".sfs","",basename(input)),"_dxy.txt", sep="" )
write.table(out.tab, out, col.names=T, row.names=F, quote=F, sep="\t")