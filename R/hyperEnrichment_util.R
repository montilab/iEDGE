#suppressPackageStartupMessages(require(heatmap.plus))

get_enrich_heatmap<-function(tab, outfile, ...){
  if(nrow(tab) == 0) return (NA)
  if(!("set" %in% colnames(tab))) return(NA)
  if(length(unique(as.character(tab[, "set"])))<2) return(NA)
  else {
    pdf(outfile)
    res<-hyper2qmatrix (
        hyper = tab,          
        fdr=c(.25,0.05,.01), 
        do.sort=TRUE,   
        do.heat=TRUE,  
        rm.zero=TRUE,  
        method=c("ward.D"), ...)
      dev.off()
      return(res)
  }

}

#######################################################################
## function: HYPER 2 QMATRIX
##
## Take the output of hyperEnrichment (called on multiple signatures)
## and generate a summary geneset-by-signature matrix (and optionally
## a heatmap) of the genesets' q-values
#######################################################################
hyper2qmatrix <- function
(
    hyper,          # output of hyperEnrichment
    fdr=c(.05,.01), # FDR thresholds (must be in decreasing order)
    do.sort=TRUE,   # sort matrices by HC
    do.heat=FALSE,  # display heatmap
    rm.zero=TRUE,   # remove genesets/rows w/ no hits
                    # hclust methods
    method=c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid"),
    ...             # extra arguments to my.heatmap
)
{
    if ( length(fdr)>1 && any(diff(fdr)>0) )
        stop( "fdr must be in decreasing order" )
    if ( length(unique(hyper[,"set"]))==1 )
        stop( "hyper must contain results for at least two signatures")
    method <- match.arg(method)
    
    SIG <- as.character(unique(hyper[,"set"]))
    gsets <- hyper[hyper[,"set"]==SIG[1],"category"]
    mx <- sapply(SIG,function(z) {
        tmp <- hyper[hyper[,"set"]==z,c("fdr","category")]
        as.numeric(tmp[match(gsets,tmp[,"category"]),"fdr"])
    });
   # print(SIG)
   # mx<-data.frame(mx) #added
    rownames(mx) <- gsets
    colnames(mx)<-SIG #added
    #print(colnames(mx))
    mx01 <- mx
    colnames(mx01)<-SIG #added
    for ( i in 1:length(fdr) ) {
        mx01[mx<=fdr[i]] <- i
    }
    mx01[mx>fdr[1]] <- 0

    if ( rm.zero ) {
        keep.idx <- apply(mx01!=0,1,any)
        if ( sum(keep.idx)==0 ) {
            warning("no significant genesets found")
        }
        else {
            mx01 <- mx01[keep.idx,,drop=FALSE]
            mx <- mx[keep.idx,,drop=FALSE]
        }
    }
    if(dim(mx01)[1]<=2 | dim(mx01)[2]<=2){
      return(NA)
    }
    ## sort rows and columns by HC
    if ( do.sort || do.heat ) {
        hc.col <- hcopt(dist(t(mx01),method="euclidean"),method="ward.D")
        hc.row <- hcopt(dist(mx01,method="euclidean"),method="ward.D")

        if ( do.heat ) {
            ncolors <- length(unique(as.vector(mx01)))

            nr<-nrow(mx01)
            nc<-ncol(mx01)

            cmin<-1
            ccap<-15
            cr<-min(ccap/nr, cmin)
            cc<-min(ccap/nc, cmin)
            margins<-c(13, 23)

            my.heatmap(mx01,Rowv=as.dendrogram(hc.row),Colv=as.dendrogram(hc.col),scale="none",
                       col=col.gradient(c("white","red"),length=ncolors),revC=TRUE,margins = margins, 
                       cexRow = cr, cexCol = cc)
        }
        if ( do.sort ) {
            mx <- mx[hc.row$order,hc.col$order]
            mx01 <- mx01[hc.row$order,hc.col$order]
        }
    }
    return( list(mx01=mx01,mx=mx) )
}


#' @import heatmap.plus
#' @export
my.heatmap <- function
(dat,
 col=NA,
 rng=c(NA,NA),
 Rowv=NA,
 Colv=NA,
 jpg=NULL,
 do.x11=(names(dev.cur())=="null device" || names(dev.cur())=="X11") && do.palette,
 negative=F,
 do.rev=F,
 color.code=FALSE,
 do.palette=F,
 pal.names=NULL,
 ...
 )
{

  if(is.na(col)){
    ramp.br <- colorRamp(c( "blue","white","red"))
    palette.br <- rgb( ramp.br(seq(0, 1, length = 14)), max = 255)
    col<-palette.br
  }
  if (is.null(col)) col <- heat.colors(12)
  rng[is.na(rng)] <- c(min(dat,na.rm=T),max(dat,na.rm=T))[is.na(rng)]
  pal <- seq(rng[1], rng[2], (rng[2]-rng[1])/(length(col)-1) )
  labCol <- round(pal,2); if (negative) labCol <- rev(labCol)

  #heatmap( if (negative) max(dat,na.rm=T)-dat else dat, col=col, Rowv=Rowv, Colv=Colv, ... )
  heatmap.plus( if (negative) max(dat,na.rm=T)-dat else dat, col=col, Rowv=Rowv, Colv=Colv, ... )

  if ( color.code ) {
    CSCn <- apply(CSC01,2,function(z) match(z,unique(z)))
    CSCn[,-1] <- t(t(CSCn[,-1,drop=FALSE])+cumsum(apply(CSCn[,-ncol(CSCn),drop=FALSE],2,max)))
    CSCncol <- unlist(apply(CSC01,2,unique))
  }
  if (do.x11) {
    x11()
  }
  else if ( !is.null(jpg) && do.palette ) {
    dev.off()
    jpeg(jpg)
  }
  if (do.palette)
    heatmap( rbind(pal,pal), labRow=NA, col=col, Rowv=NA, Colv=NA,cexCol=1.00,scale="n",
             labCol=if(is.null(pal.names)) labCol else pal.names, margins=c(40,5) )
  if ( !is.null(jpg) && do.palette ) dev.off()
}
sd.map <- function(dat,max.sd=3,n.sd=1,robust=FALSE)
{
  SD <- apply(dat,1,if(robust) mad else sd)
  MN <- apply(dat,1,if(robust) median else mean)
  TMP <- 
    sapply(1:nrow(dat),function(i)
             cut(dat[i,],
                 breaks=unique(c(-Inf,MN[i]+seq(from=-(max.sd*SD[i]-SD[i]/(2*n.sd)),to=-(SD[i]/(2*n.sd)),by=SD[i]/n.sd),
                   MN[i]+seq(from=SD[i]/(2*n.sd),to=max.sd*SD[i]-SD[i]/(2*n.sd),by=SD[i]/n.sd),+Inf)),
                 include.lowest=T,labels=FALSE))
  dat01 <- dat
  dat01[,] <- t(TMP)
  dat01
}
value.cap <- function(dat,qnt=.95)
{
  cap <- quantile(as.vector(dat),probs=qnt)
  dat[dat>cap] <- cap
  dat
}
clustColors <- function(clust,
                        ncuts,
                        COL=c('magenta','darkgreen','blue','orange','gray','black'),
                        do.matrix=TRUE)
{
    CC <- COL[cutree(clust,k=ncuts)];
    if ( do.matrix ) {
      CC <- cbind(CC,CC)
      colnames(CC) <- c('','cluster')
    }
    return(CC)
}
 