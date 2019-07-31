require('coloc')
require('data.table')
require('tidyverse')

calc.var<-function(x,ssize){
  x=x^2
  x=x/ssize
  return(x)
}

#leadloc=18652844
window=200000
gwas_size=446237
eqtl_size=117

#read csv

setwd('/Users/pblakele/Documents/ebs_projects/colocalisation/gtex_eqtl/')
sgwas<-fread('/Users/pblakele/Documents/ebs_projects/colocalisation/DBOLT_Potassium_AllleadSNPs_AND_200KB_SPAN.csv')

#unique rows based on position

sgwas<-sgwas %>% distinct(BP, .keep_all = TRUE)

#sgwas<-sgwas[CHR %in% chrom]
#gtexdata<-gtexdata[chr %in% chrom]
#gtexdata<-gtexdata[order(position),]
#sgwas<-sgwas[order(BP),]

ok<-complete.cases(gtexdata)
gtexdata<-gtexdata[ok,]

reslist=reslist[FALSE,]

leads<-sgwas %>% distinct(leadsnp, .keep_all = FALSE)
idrange=length(leads[[1]])

tissues<-c('Brain_frontal_cortex_', 'Artery_Coronary_')

for (i in 1:idrange){
  
  leadid<-leads[i]
  leadid<-as.character(leadid)
  leadloc<-sgwas[sgwas$SNP %like% leadid, ]$BP
  chrom<-sgwas[sgwas$SNP %like% leadid,]$CHR
  ofile=paste("Brain_frontal_cortex_", chrom, ".txt", sep="")
  gtexdata<-fread(ofile)
  headlist<-c('gene_id', 'chr', 'position', 'allele1', 'allele2', 'samp', 'tss_distance', 'ma_samples', 'ma_count', 'maf', 'pval_nominal', 'slope', 'slope_se')
  colnames(gtexdata)<-headlist
  gtexdata<-gtexdata %>% distinct(position, .keep_all = TRUE)
 # chrom<-sgwas[row, CHR]
  #chrom<-sgwas[sgwas$SNP %like% leadid,]$CHR
  sgwas.subs<-sgwas[CHR %in% chrom]
  sgwas.subs<-sgwas.subs[BP<(as.integer(leadloc+window)) & BP>(as.integer(leadloc-window))]
  sgwas.subs<-sgwas.subs[!rowSums(is.na(sgwas.subs)) >= 1, ]
  sgwas.subs<-transform(sgwas.subs, A1FREQ=ifelse(A1FREQ > 0.5, 1-A1FREQ, A1FREQ))
  sgwas.subs<-transform(sgwas.subs, SE=(SE^2)/gwas_size)
  
  gtexdata.subs<-gtexdata[chr %in% chrom]
  gtexdata.subs<-gtexdata.subs[position < (as.integer(leadloc+window)) & position>(as.integer(leadloc-window))]
  gtexdata.subs<-gtexdata.subs[!rowSums(is.na(gtexdata.subs)) >= 1, ]
  gtexdata.subs<-gtexdata.subs[order(position),]
  gtexdata.subs<-gtexdata.subs[position %in% sgwas.subs$BP]
  sgwas.subs<-sgwas.subs[BP %in% gtexdata.subs$position]
  sgwas.subs<-sgwas.subs[order(BP),]
  gtexdata.subs$newid<-with(gtexdata.subs, paste0(chr, "_", position, "_", allele1, "_", allele2))
  gtexdata.subs<-transform(gtexdata.subs, slope_se=(slope_se^2))
  
  my.res <- coloc.abf(dataset1=list(beta=sgwas.subs$BETA, varbeta=sgwas.subs$SE,
                    N=gwas_size,type="quant", snp=sgwas.subs$namesnp, MAF=sgwas.subs$A1FREQ),
                    dataset2=list(beta=gtexdata.subs$slope, varbeta=gtexdata.subs$slope_se,
                    N=eqtl_size,type="quant", snp=gtexdata.subs$newid, MAF=gtexdata.subs$maf)
  )
  
  cat(sprintf("\"%f\" \"%s\" \"%s\" \"%s\" \"%s\"", i, leadid, leadloc, chrom, ofile))
  #colocsum<-cbind(my.res$summary, sgwas.subs$SNP)
  colocsum<-my.res$summary
  print (colocsum)
  if (i==1){
    reslist<-colocsum
  }
  reslist<-rbind(reslist, colocsum)
}

#manhattan plot
plot(gtexdata$position, -log10(gtexdata$pval_nominal), xlab='chromosome position', ylab='-log10 p-value', pch=19)

abf.plot <- function(coloc.obj, Pos=1:nrow(coloc.obj@results),
                     chr=NULL, pos.start=min(Pos), pos.end=max(Pos),
                     trait1="trait 1", trait2="trait 2") {
  
  
  d.pp1 = signif(coloc.obj$summary[3], 3)
  d.pp2 = signif(coloc.obj$summary[4], 3)
  d.pp4 = signif(coloc.obj$summary[6], 3)
  df = coloc.obj$results
  
  df$pp1 <- exp(df$lABF.df1 - logsum(df$lABF.df1))
  df$pp2 <- exp(df$lABF.df2 - logsum(df$lABF.df2))
  df$pos <- Pos
  
  df <- melt(df[,c("snp","pos","pp1","pp2","SNP.PP.H4")], id.vars=c("snp","pos"))
  df$variable <- sub("pp1", paste0("pp1 (", trait1, ") = ", d.pp1) ,df$variable)
  df$variable <- sub("pp2", paste0("pp2 (", trait2, ") = ", d.pp2) ,df$variable)
  df$variable <- sub("SNP.PP.H4", paste0("pp4 (Both) = ", d.pp4) ,df$variable)
  
  
  ## identify and label the top 3 SNPs that have highest pp1, pp2 or pp4
  
  df.ord = df[order(df$value, decreasing=TRUE), ]
  snps = unique(df.ord$snp)[1:3]
  
  label <- NULL # avoid R CMD check NOTE
  df$label <- ifelse(df$snp %in% snps, df$snp,"")
  ttl <- paste0(trait1, ' & ', trait2, ' (chr', chr, ': ', pos.start, '-', pos.end, ')')
  
  ggplot(df, aes_string(x="pos",y="value")) +
    geom_point(data=subset(df,label==""),size=1.5) +
    geom_point(data=subset(df,label!=""),col="red",size=1.5) +
    geom_text(aes_string(label="label"),hjust=-0.1,vjust=0.5,size=2.5,col="red") +
    facet_grid(variable ~ .) +
    theme(legend.position="none") + xlab(paste("Chromosome", chr, sep=' ')) + ylab("Posterior probability") + 
    ggtitle(ttl)
  
}
