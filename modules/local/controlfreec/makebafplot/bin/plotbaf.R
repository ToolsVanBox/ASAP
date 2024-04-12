library("ggplot2")
library("gtools")
library("graphics")
library("pastecs")

args = commandArgs(trailingOnly=TRUE)

baffile <-  type.convert(args[1])
binsize <-  type.convert(args[2])
prefix <-  type.convert(args[3])

f2 <- function(x) sum(unlist(x)[2:length(x)])
f1 <- function(x) sum(unlist(x))

determine_peaks <- function(region) {
  quant <- quantile(region$BAF)
  if(quant[2]==quant[4]) {
    # LOW COMPLEXITY
    df <- data.frame(BAF=quant[2], dens=99.99)
    df$Chromosome <- region$Chromosome[1]
    df$Position <- mean(region$Position)
    df$CALL <- NA
    return(df)
  }

  d <- density(region$BAF, bw="sj")

  ts_y<-ts(d$y)
  tp<-turnpoints(ts_y)
  peaks <- data.frame(BAF=d$x[tp$tppos], dens=d$y[tp$tppos])
  df <- data.frame(subset(peaks, dens>=2.0))
  if (nrow(df)<2) {
    # complex region
    df <- data.frame(subset(peaks, dens>=1.0))
  }
  df$Chromosome <- region$Chromosome[1]
  df$Position <- mean(region$Position)
  df$CALL <- NA

  df <- df[order(df$dens, decreasing=T),]
  return(df)
}

rowdat <- data.frame(read.table(baffile, header=T, sep='\t'))
rowdat$BAF <- as.numeric(as.character(rowdat$BAF))
rowdat$BAF[is.na(rowdat$BAF)] <- 0.0


df <- data.frame()
chromosomes <- mixedsort(as.character(unique((rowdat$Chromosome))))

for (i in c(1:length(chromosomes))) {
  chrom <- chromosomes[i]
  #print(chrom)

  tmp <- subset(rowdat, Chromosome==chrom)

  for (j in seq(1, nrow(tmp), by=binsize/2)) {
    maxy <- j+binsize-1

    if (maxy>nrow(tmp)) {maxy<-nrow(tmp)}

    region <- tmp[j:maxy,]
    if (nrow(region) >= 100) {
      peaks <- determine_peaks(region)
      #print(nrow(peaks))
      df <- rbind(df, peaks[1:2,])
    }

    # check if HETRO / HOMOZYGOUS calls present
    if (peaks$BAF[1]<=0.1 || peaks$BAF[1]>=0.9) {
      # RE-RUN for intermediate BAFs
      region <- subset(region, BAF>0.1 & BAF<0.9)
      if (nrow(region) >= 100) {
        peaks <- determine_peaks(region)
        #print(nrow(peaks))
        df <- rbind(df, peaks)
      }
    }
  }

}

df$BAF <- round(df$BAF,3)
df <- df[!is.na(df$Chromosome),]
df$Chromosome <- factor(df$Chromosome, levels=chromosomes)
#rm(rowdat)

pdf(file=paste0(prefix,"_BAF.pdf"), width=15, height=2, pointsize=6, useDingbats=FALSE)

p <- ggplot(df, aes(Position, BAF, group=Chromosome)) + geom_point(shape=20, size=0.2, alpha=0.5) + facet_wrap(~Chromosome, nrow=1, scales="free_x")

p <- p + theme(
  axis.line.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),

  panel.grid.minor=element_blank(),
  panel.grid.major.x=element_blank(),
  panel.grid.major.y=element_line(colour="grey80", linetype=2, linewidth=0.2),

  legend.position="none",
  panel.background=element_blank(),
  panel.border=element_blank(),
  plot.background=element_blank())

print(p)
dev.off()
rm(p)