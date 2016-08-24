##same output format as GSA.read.gmt except without the annoying text outputs
#' \code{read_gmt} read from gmt file
#' @param f name of gmt file
#' @param split.char delimiter
#' export
read_gmt<-function(f, split.char = '\t'){

  con <- file(f, open = 'r') 

  results.list <- list()
  current.line <- 1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    results.list[[current.line]] <- as.character(unlist(strsplit(line, split=split.char)))
    current.line <- current.line + 1
  } 
  close(con)

  genesets<-lapply(results.list, function(x) x[-(1:2)])
  geneset.names<-unlist(lapply(results.list, function(x) x[1]))
  geneset.descriptions<-unlist(lapply(results.list, function(x) x[2]))
  names(genesets)<-geneset.names
  res<-list(genesets = genesets, geneset.names = geneset.names, geneset.descriptions = geneset.descriptions)
  return(res)
}


run_hyperEnrichment_unpruned<-function(ngenes, gs, gs.file.name, drawnList, f.dir.out, header, header2){
  #run hyperenrichment with all sig cis genes in one geneset

  res<-run_hyperEnrichment(drawn=drawnList,
      categories=gs,
      ntotal=ngenes,
      min.drawsize = min.drawsize, mht = TRUE, verbose = TRUE, order = TRUE)
  f.out<-paste(f.dir.out, "/", header, 
    "_hyperEnrichment_",gs.file.name,"_", header2,".txt", sep = "")  
  cat(paste("Writing table to ", f.out, "\n", sep = ""))
  write.table(res, sep = "\t", col.names = TRUE, row.names = FALSE,
    file = f.out)
  cat("Done!\n")
  return(res)

}

run_hyperEnrichment_pruned<-function(tab, tab.name, gs, ngenes, 
  min.drawsize = 3, 
  hypercol = "fdr", 
  hyperthres = 0.25, 
  f.dir.out){

  drawnList<-list(drawn = unique(tab$trans))
  hyperunion<-run_hyperEnrichment(drawn=drawnList,
    categories=gs,
    ntotal=ngenes,
    min.drawsize = min.drawsize, 
    mht = TRUE, 
    verbose = TRUE, 
    order = TRUE, 
    hypercol = hypercol, hyperthres = hyperthres)
  
  if(nrow(hyperunion)>0)
  write.table(hyperunion, file = paste(f.dir.out, "/", tab.name, ".txt", sep = ""), 
    sep = "\t", row.names = FALSE)
    
  cis<-unique(tab$cis)
  hyperbyalt<-lapply(cis, function(j){
      i.sub<-tab[tab$cis %in% j,]
      drawnList <- list(drawn = unique(i.sub$trans))
      hyper<-run_hyperEnrichment(drawn=drawnList,
          categories=gs,
          ntotal=ngenes,
          min.drawsize = min.drawsize, 
          mht = TRUE, 
          verbose = TRUE, 
          order = TRUE,
          hypercol = hypercol, hyperthres = hyperthres)
      
      if(nrow(hyper)>0)
      write.table(hyper, file = paste(f.dir.out, "/", tab.name, "_", j, ".txt", sep = ""), 
            sep = "\t", row.names = FALSE)
      return(hyper)
      })
  names(hyperbyalt)<-cis

  return(list(hyperunion = hyperunion, hyperbyalt = hyperbyalt))
}


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
  
  res.HE.groupalt<-list()
 
  for (j in names(gs.tab)){
    de.table<-gs.tab[[j]]
    alt<-unique(de.table[, "alteration_descriptor"])
    gs.custom.groupalt<-sapply(alt,
      function(x) 

      as.character(de.table[de.table[, "alteration_descriptor"] == x, "gene_id"])

      #as.character(subset(de.table, alteration_descriptor == x)[, "gene_id"])

      )
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
 
  cnames <-c("set", "pval","fdr","set annotated","set size","category annotated","total annotated","category","hits")
  nodata<-as.data.frame(setNames(replicate(length(cnames),numeric(0), simplify = F), cnames))

  res<-list(single = res.HE.single, groupalt=res.HE.groupalt.list)

  for (i in names(res)){
    for(j in names(res[[i]])){
      if(is.null(res[[i]][[j]])){
        res[[i]][[j]]<-nodata
      }
    }
  }
  return(res)
}

#' \code{run_hyperEnrichment} performs hyperenrichment test from gistic2ge DE results on gene sets

#' @param drawn list of genesets, usually from user-defined queries
#' @param categories list of genesets, usually from pathway databases
#' @param ntotal total number of genes 
#' @return limma topTable data frame
#' @export
#'
run_hyperEnrichment<-function(drawn, categories, ntotal, 
  min.drawsize = 4, mht = TRUE, verbose = FALSE, order = TRUE,
  hypercol = "fdr", hyperthres = NULL){

  if(length(drawn) == 0){
    df.names<-c("set", "pval","fdr","set.annotated","set.size","category.annotated",
    "total.annotated","category", "hits")
    df<-as.data.frame(rep(list(a=numeric(0)), length(df.names)))
    colnames(df)<-df.names
    return(df)
  }
  res.HE<-hyperEnrichment(
      drawn=drawn,         
      categories = categories,   
      ntotal=ntotal,
      min.drawsize=min.drawsize, 
      mht=mht,       
      verbose=verbose)
  if(length(drawn) == 0) return(NULL)
  res.HE.df<-data.frame(res.HE)
  if (nrow(res.HE.df) == 0) return(res.HE.df)
  num_col<-c("pval", "fdr", "set.annotated", "set.size", "category.annotated", "total.annotated")

  for (i in num_col)
    res.HE.df[,i] <-suppressWarnings(as.numeric(paste(res.HE.df[, i])))

  res.HE.df<-res.HE.df[apply(res.HE.df, 1, function(x) !all(is.na(x))),]

  if (order)
    res.HE.df<-res.HE.df[order(res.HE.df$fdr, decreasing = FALSE),]

  if(hypercol %in% colnames(res.HE.df) & !is.null(hyperthres)){
   
    res.HE.df<-res.HE.df[res.HE.df[, hypercol]<=hyperthres,]
  }
    
  return(res.HE.df)
}

hyperEnrichment<-function (drawn, categories, ntotal = length(unique(unlist(categories))), 
    min.drawsize = 4, mht = TRUE, verbose = FALSE) {
    if (!is(categories, "list")) {
        stop("categories expected to be a list of gene sets")
    }
    gene.names <- unique(unlist(categories))
    if (is.list(drawn) && is.null(names(drawn))) {
        stop("drawn must have non-null names when a list")
    }
    if (ntotal < length(unique(unlist(categories)))) {
        warning("background population's size less than unique categories' items: ", 
            ntotal, "<", length(gene.names))
    }
    cnames <- c("pval", "fdr", "set annotated", "set size", "category annotated", 
        "total annotated", "category", "hits")
    if (is.list(drawn)) {
        ncat <- length(categories)
        enrich <- matrix(NA, ncat * length(drawn), length(cnames) + 
            1)
    
        VERBOSE(verbose, "Testing", length(drawn), "drawsets on", 
            ncat, "categories and", length(gene.names), "total items ..\n")
        percent <- 0.1
        base <- 0
        ntst <- 0
        for (i in 1:length(drawn)) {
            VERBOSE(verbose, "*** Testing", names(drawn)[i], 
                ".. ")
            dset <- drawn[[i]]
            tmp <- hyperEnrichment(dset, categories, ntotal = ntotal, 
                verbose = verbose)
            if (is.null(tmp)) {
                VERBOSE(verbose, "not enough items drawn\n")
                next
            }
            ntst <- ntst + 1
            rng <- (base + 1):(base + ncat)
            if (any(!is.na(enrich[rng, ]))) 
                stop("something wrong")
            enrich[rng, ] <- cbind(set = rep(names(drawn)[i], 
                ncat), tmp)
            base <- base + ncat
            if (F && i >= round(length(drawn) * percent)) {
                VERBOSE(verbose, round(100 * percent), "% ", 
                  sep = "")
                percent <- percent + 0.1
            }
            VERBOSE(verbose, " (min fdr: ", signif(min(suppressWarnings(as.numeric(tmp[, 
                "fdr"]))), 2), ")\n", sep = "")
        }
        VERBOSE(verbose, "done.\n")
        colnames(enrich) <- c("set", cnames)
        enrich <- enrich[1:base, , drop = F]
        if (mht) {
            VERBOSE(verbose, "MHT-correction across multiple draws ..")
            enrich[, "fdr"] <- pval2fdr(suppressWarnings(as.numeric(enrich[, "pval"])))
            VERBOSE(verbose, "done.\n")
        }
        VERBOSE(verbose, "Categories tested: ", rjust(length(categories), 
            4), "\n", "Candidate sets:    ", rjust(length(drawn), 
            4), "\n", "Sets tested:       ", rjust(ntst, 4), 
            "\n", "Items tested:      ", rjust(sum(sapply(drawn, 
                length)), 4), " (min,med,max: ", paste(quantile(sapply(drawn, 
                length), probs = c(0, 0.5, 1)), collapse = ","), 
            ")\n", "N(FDR<=0.25):      ", rjust(sum(enrich[, 
                "fdr"] <= 0.25), 4), "\n", "N(FDR<=0.05):      ", 
            rjust(sum(enrich[, "fdr"] <= 0.05), 4), "\n", "N(FDR<=0.01):      ", 
            rjust(sum(enrich[, "fdr"] <= 0.01), 4), "\n", sep = "")
        return(enrich)
    }
    m.idx <- drawn[drawn %in% gene.names]
    if (length(m.idx) < min.drawsize) {
        VERBOSE(verbose, "insufficient annotated genes in the drawn set: ", 
            paste(gene.names[m.idx], collapse = ","), "\n")
        return(NULL)
    }
    VERBOSE(verbose, length(m.idx), "/", length(drawn), " annotated genes found", 
        sep = "")
    nhits <- sapply(categories, function(x, y) length(intersect(x, 
        y)), m.idx)
    ndrawn <- length(drawn)
    ncats <- sapply(categories, length)
    nleft <- ntotal - ncats
    enrich <- phyper(q = nhits - 1, m = ncats, n = nleft, k = ndrawn, 
        lower.tail = F)
    enrich <- cbind(pval = enrich, fdr = pval2fdr(enrich), nhits = nhits, 
        ndrawn = ndrawn, ncats = ncats, ntot = ntotal, category = names(categories))
    enrich <- cbind(enrich, hits = sapply(categories, function(x, 
        y) paste(intersect(x, y), collapse = ","), m.idx))
    ord <- order(suppressWarnings(as.numeric(enrich[, "pval"])))
    enrich <- enrich[ord, , drop = F]
    enrich[, "pval"] <- signif(suppressWarnings(as.numeric(enrich[, "pval"])), 
        2)
    enrich[, "fdr"] <- signif(suppressWarnings(as.numeric(enrich[, "fdr"])), 2)
    colnames(enrich) <- cnames
    rownames(enrich) <- names(categories)[ord]
    return(enrich)
}


