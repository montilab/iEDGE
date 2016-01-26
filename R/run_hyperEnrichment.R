

#' \code{run_hyperEnrichment_wrapper} wrapper for run_hyperEnrichment
#' @param gs.tab DE results table from gistic2ge
#' @param ntot total genes
#' @param gs.pathway list of genesets
#' @return list of data frames of hyperenrichment results
#' @export
#'
run_hyperEnrichment_wrapper<-function(gs.tab, #gistic2ge_sig
  ntot, #total genes for background correction
  gs.pathway #list of genesets
  ){

  res.HE.single<-list()
  for (j in names(gs.tab)){
    de.table<-gs.tab[[j]]
    gs.custom.single<-list(as.character(unique(de.table[, "gene_id"])))
    names(gs.custom.single)<-paste(j, "single", sep = "_")    
    res.HE.single[[j]]<-run_hyperEnrichment(drawn = gs.custom.single, 
      categories = gs.pathway, 
      ntotal = ntot)
  }
  
  #res.HE.single.sig<-lapply(res.HE.single, function(x) 
  #  if (is.null(x)) return(NULL) else return(subset(x, fdr < hyper.fdr)))

  res.HE.groupalt<-list()
 
  for (j in names(gs.tab)){
    de.table<-gs.tab[[j]]
    alt<-unique(de.table$alteration_descriptor)
    gs.custom.groupalt<-sapply(alt,
      function(x) as.character(subset(de.table, alteration_descriptor == x)[, "gene_id"]))
    names(gs.custom.groupalt)<-paste(j, alt, sep = "_")
    res.HE.groupalt[[j]]<-list()
    for (k in names(gs.custom.groupalt)){
      gs.custom<-list(gs.custom.groupalt[[k]])
      names(gs.custom)<-paste(k, sep = "_")
      res.HE.groupalt[[j]][[k]]<-run_hyperEnrichment(drawn = gs.custom, 
        categories = gs.pathway,
        ntotal = ntot)
    }
  }
  res.HE.groupalt.list<-lapply(res.HE.groupalt[!sapply(res.HE.groupalt, is.null)], function(x){
    do.call(rbind, x[sapply(x, function(y) !all(is.na(y)))])
  })
 
#  res.HE.groupalt.list.sig<-lapply(res.HE.groupalt.list, function(x) 
 #   if (is.null(x)) return(NULL) else return(subset(x, fdr < hyper.fdr)))

  cnames <-c("set", "pval","fdr","set annotated","set size","category annotated","total annotated","category","hits")
  nodata<-as.data.frame(setNames(replicate(length(cnames),numeric(0), simplify = F), cnames))

  res<-list(single = res.HE.single, groupalt=res.HE.groupalt.list)
#  res<-list(res.single = res.HE.single, res.groupalt=res.HE.groupalt.list,
 #   res.single.sig = res.HE.single.sig, res.groupalt.sig = res.HE.groupalt.list.sig)
  
  for (i in names(res)){
    for(j in names(res[[i]])){
      if(is.null(res[[i]][[j]])){
        res[[i]][[j]]<-nodata
      }
    }
  }
  return(res)
}


subset_hyperEnrichment<-function(x, #list of data frames
  hyper.fdr){
  res<-list()
  for(i in names(x))
    res[[i]]<-subset(x[[i]], fdr < hyper.fdr)
  return(res)
}

#' \code{run_hyperEnrichment} performs hyperenrichment test from gistic2ge DE results on gene sets
#' @import CBMRtools
#' @param drawn list of genesets, usually from user-defined queries
#' @param categories list of genesets, usually from pathway databases
#' @param ntotal total number of genes 
#' @return limma topTable data frame
#' @import CBMRtools
#' @export
#'
run_hyperEnrichment<-function(drawn, categories, ntotal, min.drawsize = 4, mht = TRUE, verbose = TRUE, order = TRUE){
  require(CBMRtools)
  res.HE<-hyperEnrichment(
      drawn=drawn,         
      categories = categories,   
      ntotal=ntotal,
      min.drawsize=min.drawsize, 
      mht=mht,       
      verbose=verbose)

  res.HE.df<-data.frame(res.HE)
  num_col<-c("pval", "fdr", "set.annotated", "set.size", "category.annotated", "total.annotated")

  for (i in num_col)
    res.HE.df[,i] <-as.numeric(paste(res.HE.df[, i]))

  if (order)
    res.HE.df<-res.HE.df[order(res.HE.df$fdr, decreasing = FALSE),]
    
  return(res.HE.df)
}


run_test<-function(gistic2ge_sig){

  st.res<-read.table("~/Desktop/git_projects/datasets/ricover2010/GISTIC_20110620.4amy/S2X_20110625/HyperEnrichment.MSigDBc2_cp.v2.5.wexp.fdr25.xls",
   sep = "\t", header = TRUE)
  st.res.category<-as.character(st.res$category)
  
  sig.df<-rbind(gistic2ge_sig[["peak_amp_cis_limma"]], gistic2ge_sig[["peak_del_cis_limma"]])
  #thres<-log2(log2(130)/log2(100))
  #sig.df<-subset(sig.df, logFC.lm > thres | logFC.lm < (-thres))

  drawnList<-list(amp_del_region = unique(as.character(sig.df[["gene_id"]])))

  ##hypernerichment

  c2.file<-"/Users/amyli/Desktop/monti_lab/year2/DataIntegration/annot/msigdb_v2.5/c2.cp.v2.5.symbols.gmt"
  c2<-read_gmt(c2.file)$genesets

  c2<-sapply(c2, function(x){
    a<-lapply(x,function(i) strsplit(i, split="[ |_]///[ |_]")[[1]])
    return(unique(do.call(c,a)))
    })

  categoriesList<-c2
  #hyperE <- hyperEnrichment(drawn=hyperSig,categories=getGeneSet(hyperGsets),ntotal=10000)
  hyperE <- run_hyperEnrichment(drawn=drawnList,
    categories=categoriesList,
    ntotal=nrow(gistic_ge_amp$ge$mat),
    min.drawsize = 4, mht = TRUE, verbose = TRUE)
  hyperE.fdr<-as.numeric(hyperE[, "fdr"])
  hyperE.fdr.sig<-hyperE[which(hyperE.fdr < 0.25),]
  hyperE.fdr.sig
  dim(hyperE.fdr.sig)  

  sig.cats<-hyperE.fdr.sig[, "category"]

  length(st.res.category)
  length(sig.cats)
  length(intersect(st.res.category, sig.cats))



}



