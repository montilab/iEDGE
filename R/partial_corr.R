

wrapper_partial_corr<-function(){
require(ppcor)
n =100
m = 0
s = 1
a<-rnorm(n = n, mean = m, sd = s)
b<-sapply(a, function(x) rnorm(n = 1, mean = x, sd = 0.5))
c<-sapply(a, function(x) rnorm(n = 1, mean = x, sd = 0.5))

pcor.test(a, b, c)


x.res<-lm(formula =Air.Flow ~ Acid.Conc., data = stackloss)$resid
y.res<-lm(formula = Water.Temp ~ Acid.Conc., data = stackloss)$resid
cor(x.res, y.res)

pcor.test(stackloss$Air.Flow, stackloss$Water.Temp, stackloss$Acid.Conc.)


###get data

gistic_ge_amp<-readRDS("~/Desktop/git_projects/gistic2ge_results/lymphoma2010/lymphoma2010_gistic_ge_amp.RDS")
gistic_ge_del<-readRDS("~/Desktop/git_projects/gistic2ge_results/lymphoma2010/lymphoma2010_gistic_ge_del.RDS")

df<-readRDS("~/Desktop/git_projects/gistic2ge_results/lymphoma2010/lymphoma2010_gistic2ge_sig_cis.fdr_0.25_trans.fdr_0.05_sig.FC_NA.RDS")

df.cis<-df$peak_cis



for (i in unique(df$peak_cis$alteration_id)[1]){
	print(i)
	cis.genes<-as.character(subset(df$peak_cis, alteration_id == i, select = gene_id)$gene_id)
	trans.genes<-as.character(subset(df$peak_trans, alteration_id == i, select = gene_id)$gene_id)
	print(cis.genes)
	if( grepl("Amp", i)){
		alt.ind<-which(gistic_ge_amp$gistic$meta$Unique.Name == i)
		alt<-as.numeric(gistic_ge_amp$gistic$mat.cont[alt.ind, ])
		cis<-gistic_ge_amp$ge$mat[cis.genes,]
		trans<-gistic_ge_amp$ge$mat[trans.genes,]
	} else {
		alt.ind<-which(gistic_ge_del$gistic$meta$Unique.Name == i)
		alt<-as.numeric(gistic_ge_del$gistic$mat.cont[alt.ind, ])
		cis<-gistic_ge_del$ge$mat[cis.genes,]
		trans<-gistic_ge_del$ge$mat[trans.genes,]		
	}

}



cis.cor<-sapply(1:nrow(cis), function(i) cor(x=alt, y=as.numeric(log(cis[i,]))))
trans.corr<-sapply(1:nrow(trans), function(i) cor(x=alt, y=as.numeric(log(trans[i,]))))

}
