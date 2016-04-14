## This file contains R functions for plotting schematics of genomic regions, similar to a 
## UCSC genome browser shot. Right now it can plot gene bodies, bar plots (for single base
## phyloP scores), and rectangles to indicate genomic regions (like an element from a bed 
## file). I'll add code to plot ChIP peak data (like the layered H3K27ac track later).
## By Rachel Gittelman

##########################################################################################
## set region boundaries
##########################################################################################
chr="chr12"
start=113340750
stop=113464230

##########################################################################################
## plotting functions
##########################################################################################

get_window_average <- function(data, start, size) {

	## This will be a helper function with plot_bar. If there are
	## too many bars to plot, it will instead groups bars into windows
	## and average them.
	
	stop = start + size - 1
	overlap_beginning <- intersect(which(data[,2] >= start), which(data[,2] <= stop))
	overlap_ending <- intersect(which(data[,3] <= stop), which(data[,3] >= start))
	overlap = union(overlap_beginning, overlap_ending)
	#print(overlap)
	data <- data[sort(overlap),]
	return(mean(data[,4]))
}

plot_bar <- function(data, chr, start, stop) {
	
	## This function will plot a bar plot. If it covers a range of more than 1000
	## bases, it will divide the data into 1000 equally spaced regions and average
	## the scores for each region, plotting just those.
	
	data <- data[data[,1] == chr,]
	overlap_beginning <- intersect(which(data[,2] >= start), which(data[,2] <= stop))
	overlap_ending <- intersect(which(data[,3] <= stop), which(data[,3] >= start))
	overlap = union(overlap_beginning, overlap_ending)
	data <- data[sort(overlap),]
	par(mar=c(.5,5,.5,0))

	if(stop - start > 1000)
	{
		window_length = round((stop-start)/1000)
		window_starts <- seq(start,stop,window_length)
		heights <- unlist(lapply(window_starts, get_window_average, data=data, size=window_length))
		plot(c(),xlim=c(start,stop),ylim=c(min(heights),max(heights)),axes=F,xlab="",ylab="", bty="n")
		rect(xleft=window_starts,xright=window_starts+window_length,ybottom=rep(0,length(window_starts)), ytop=heights,col="black")
	}
	else
	{
		plot(c(),xlim=c(start,stop),ylim=c(min(data[,4]),max(data[,4])),axes=F,xlab="",ylab="", bty="n")
		rect(xleft=data[,2]-.5,xright=data[,2]+.5,ybottom=rep(0,length(data[,4])), ytop=data[,4],col="black")
	}
	axis(2,cex.axis=1.5)
}

plot_canon_gene <- function(exon_data, fiveu_data, threeu_data, gene_data, chr, start, stop)
{
	
	## This function will plot gene bodies. It requires 4 bed files, including exons, 5 prime
	## UTRs, 3 prime UTRs, and the overall transcription start and stop coordinates of the genes.
	## Finally, it requires a chromosome (character) and a start and stop for the plotting
	## window (start and stop are numeric)
	
	all_data <- list(exon_data, fiveu_data, threeu_data, gene_data)
	for(i in 1:4)
	{
		data <- all_data[[i]][all_data[[i]][,1] == chr,]
		overlap_beginning <- intersect(which(data[,2] >= start), which(data[,2] <= stop))
		overlap_ending <- intersect(which(data[,3] <= stop), which(data[,3] >= start))
		overlap <- union(overlap_beginning,overlap_ending)
		overlap_whole <- intersect(which(data[,2] <= start), which(data[,3] >= stop))
		overlap <- union(overlap,overlap_whole)
		all_data[[i]] <- data[sort(overlap),]
	}
	##print(all_data)
	par(mar=c(.5,5,.5,0))
	plot(c(),xlim=c(start,stop),ylim=c(0,1), axes=F,xlab="",ylab="")
	if(nrow(all_data[[4]]) > 0){
		par(xpd=F)
		segments(x0=all_data[[4]][,2],x1=all_data[[4]][,3],y0=.5,y1=.5, lwd=2)}
		par(xpd=NA)
	if(nrow(all_data[[1]]) > 0){
		rect(col="blue",xleft=all_data[[1]][,2],xright=all_data[[1]][,3],ybottom=rep(0.05,nrow(all_data[[1]])),ytop=rep(.95,nrow(all_data[[1]])))}
	if(nrow(all_data[[2]]) > 0){
		rect(col="blue",xleft=all_data[[2]][,2],xright=all_data[[2]][,3],ybottom=rep(0.25,nrow(all_data[[2]])),ytop=rep(.75,nrow(all_data[[2]])))}
	if(nrow(all_data[[3]]) > 0){
		rect(col="blue",xleft=all_data[[3]][,2],xright=all_data[[3]][,3],ybottom=rep(0.25,nrow(all_data[[3]])),ytop=rep(.75,nrow(all_data[[3]])))}
}

plot_rect <- function(data, chr, start, stop) {

		## This function will just plot rectangles given a "data" dataframe in the following
		## format:
		
		##"chr1" start stop
		
		## and a chr, start, and stop indicating limits of plotting

        data <- data[data[,1] == chr,]
        overlap_beginning <- intersect(which(data[,2] >= start), which(data[,2] <= stop))
        overlap_ending <- intersect(which(data[,3] <= stop), which(data[,3] >= start))
        overlap = union(overlap_beginning, overlap_ending)
        data <- data[sort(overlap),]
        print(data)
        if(nrow(data) > 0){
                par(mar=c(.5,5,.5,0))
                plot(c(),xlim=c(start,stop),ylim=c(0,1), axes=F,xlab="",ylab="")
                rect(col="black",xleft=data[,2],xright=data[,3],ybottom=rep(0.1,nrow(data)),ytop=rep(.9,nrow(data)))}
}

##########################################################################################
## read in data
##########################################################################################
library(seqminer)

## read in gene coordinate data
exon_data <- read.table("/path/to/coding_exons.canon.bed",sep="\t", header=F)
fiveutr_data <- read.table("/path/to/5primeUTR.canon.bed",sep="\t", header=F)
threeutr_data <- read.table("/path/to/3primeUTR.canon.bed",sep="\t", header=F)
gene_data <- read.table("/path/to/canonical.bed",sep="\t", header=F)

## read in phyloP data
phyloP_data <- tabix.read("/path/to/phyloP.bed.gz", paste(chr,":",start,"-",stop,sep=""))
phyloP_data <- data.frame(do.call(rbind,strsplit(phyloP_data,"\t")), stringsAsFactors=F)
phyloP_data[,2] <- as.numeric(phyloP_data[,2])
phyloP_data[,3] <- as.numeric(phyloP_data[,3])
phyloP_data[,4] <- as.numeric(phyloP_data[,4])

## read in 7 H3K27ac data


##########################################################################################
## plot example locus (OAS on chr 12)
##########################################################################################

layout(matrix(c(1,2,3),ncol=1))
par(oma=c(5,1,1,1))
plot_rect(gene_data,chr,start,stop)
plot_canon_gene(exon_data,fiveutr_data,threeutr_data,gene_data,chr,start,stop)
plot_bar(phyloP_data,chr,start,stop)
axis(side=1,cex=3, lwd=3, cex.axis=2)