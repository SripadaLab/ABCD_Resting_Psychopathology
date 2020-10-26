library(circlize)

ZTHRESH = 2
PercThresh = 7.5
IncludeList = c(1,3,4,5,7,8,11,12)
Filename = "consensus"
mat = as.matrix(read.table("Results/consensus_z_square.txt",header=F,sep=""))
plot_circlegraph(mat,ZTHRESH,PercThresh,Range=180,Filename,IncludeList)

Filename = "component005"
mat = as.matrix(read.table("Results/component005_z_square.txt",header=F,sep=""))
plot_circlegraph(mat,ZTHRESH,PercThresh,400,Filename,IncludeList)

IncludeList = c(1,3,5,7,8,11,12)
Filename = "component139"
mat = as.matrix(read.table("Results/component139_z_square.txt",header=F,sep=""))
plot_circlegraph(mat,ZTHRESH,PercThresh,180,Filename,IncludeList)

build_color_function = function(max) {
  color_function = function(value) {
    percent = abs(value)/max
    r = floor(approx(c(0,1),c(255,0),percent)$y)
    g = floor(approx(c(0,1),c(255,0),percent)$y)
    b = floor(approx(c(0,1),c(255,0),percent)$y)
    a = floor(approx(c(0,1),c(0,255),percent)$y)
    r[value>0] = 255
    b[value<0] = 255
    rgb(r,g,b,a,maxColorValue = 255)
  }
  color_function
}

#if Plot==TRUE, it ignores Filename and just plots the image in the R graphics window
plot_circlegraph = function(mat,ZTHRESH,PercThresh,Range=NULL,Filename,IncludeList=1:16,Plot=TRUE) {
  NetworkNames = c('SM','SMM','CO','AUD','DMN','','VIS','FPN','SAL','SC','VAN','DAN','CER','UNC','CP','RST')
  NetworkNamesFactor = factor(1:16,levels=1:16,labels=NetworkNames)
  
  NetworkColors = c(
    '#00FFFF55',
    '#FF800055',
    '#80008055',
    '#FF00FF55',
    '#FF000055',
    '#80808055',
    '#0000FF55',
    '#FFFF0055',
    '#00000055',
    '#B43C0055',
    '#00808055',
    '#00FF0055',
    '#80FFFF55',
    '#AAAAAA55',
    '#80408055',
    '#80808055')
  names(NetworkColors) = NetworkNames
    
  dat = read.csv("Data/gordon_sub_cere_parcels.csv")
  dat = dat[1:418,]
  
  mat[abs(mat)<ZTHRESH] = 0
  mat = mat + t(mat)
    
  o = order(dat$NetworkNumber)
  onets = dat$NetworkNumber[o]
  onets[onets==2] = 1 #combine somatomotor networks for this vis
  olabels = dat$Community[o]
  
  smallmat = mat[o,o]
  #image(smallmat)
  
  smallmat = smallmat[onets %in% IncludeList,onets %in% IncludeList]
  #image(smallmat)
  
  smallmat[smallmat>0] = 1
  smallmat[smallmat<0] = -1
  rownames(smallmat) = onets[onets %in% IncludeList]
  colnames(smallmat) = onets[onets %in% IncludeList]
  
  smallmat[lower.tri(smallmat)] = 0
  
  df = data.frame(from = rep(rownames(smallmat), times = ncol(smallmat)),
                  to = rep(colnames(smallmat), each = nrow(smallmat)),
                  value = as.vector(smallmat),
                  stringsAsFactors = FALSE)
  df = df[df$value!=0,]
  range(df$value)
  a1 = aggregate(df$value,by=list(from=df$from,to=df$to),function(x){sum(x>0)})
  a2 = aggregate(df$value,by=list(from=df$from,to=df$to),function(x){-1*sum(x<0)})
  a = rbind(a1,a2)
  a = a[a$x!=0,]
  
  #get cell counts
  smallonets = onets[onets %in% IncludeList]
  smallupper = upper.tri(smallmat)
  for (i in 1:nrow(a)) {
    a$total[i] = sum(smallupper[smallonets==a$from[i],smallonets==a$to[i]])
  }
  a$percent = abs(a$x)/a$total
  a$from = NetworkNames[as.integer(a$from)]
  a$to = NetworkNames[as.integer(a$to)]
  
  a_thresh = a[a$percent>(PercThresh/100),]
  DefaultRange = max(abs(range(a_thresh$x)))
  if (is.null(Range)) {
    Range = floor(1.05*DefaultRange)
  }
  col_fun = build_color_function(Range)
  op = par()
  
  if (Plot==TRUE) {
    svg(filename=paste0("Results/",Filename,".svg"), 
        width=8, 
        height=8, 
        pointsize=12)
  }
  circos.clear()
  circos.par(gap.after = c(rep(5,length(IncludeList))))
  chordDiagram(a_thresh[,1:3],
               grid.col=NetworkColors,
               col=col_fun,
               scale=F,
               self.link=1,
               order=NetworkNames[IncludeList],
               link.largest.ontop=T,
               annotationTrack = c("grid"),
               annotationTrackHeight = c(0.15))
  for(si in get.all.sector.index()) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
                facing = "bending.outside", niceFacing = T, col = "black",font=2,cex=2)
  }
  par(op)
  if (Plot==TRUE) {
    dev.off()
  }
}
