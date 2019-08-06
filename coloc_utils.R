require('coloc')
require('data.table')
require('tidyverse')

calc.var<-function(x,ssize){
  x=x^2
  x=x/ssize
  return(x)
}

#leadloc=18652844
#window=200000
#gwas_size=446237
#eqtl_size=117

#gtexdata<-fread('/Users/pblakele/Documents/ebs_projects/colocalisation/Brain_Frontal_Cortex_BA9.allpairs_2.txt')
#sgwas<-fread('/Users/pblakele/Documents/ebs_projects/colocalisation/DBOLT_Potassium_AllleadSNPs_AND_200KB_SPAN.csv')

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=6) {
  stop("At least three argument must be supplied (input file).n", call.=FALSE)
}

working.dir <- args[3]
# change the working dir
setwd(file.path(working.dir))

sgwas = fread(args[1], header=TRUE)
tissues = fread(args[2], header=TRUE)

window <- as.integer(args[4])
gwas_size <- as.integer(args[5])
eqtl_size <- as.integer(args[6])

working.dir2<-paste(working.dir, 'gtex_eqtl', sep='/')
setwd(file.path(working.dir2))

#unique rows based on position
sgwas<-sgwas %>% distinct(BP, .keep_all = TRUE)

#sgwas<-sgwas[CHR %in% chrom]
#gtexdata<-gtexdata[chr %in% chrom]
#gtexdata<-gtexdata[order(position),]
#sgwas<-sgwas[order(BP),]

#ok<-complete.cases(gtexdata)
#gtexdata<-gtexdata[ok,]

#reslist=reslist[FALSE,]

leads<-sgwas %>% distinct(leadsnp, .keep_all = FALSE)
idrange=length(leads[[1]])

tissues<-c('Brain_frontal_cortex', 'Artery_Coronary')
tissue_range<-length(tissues)

for (i in 1:idrange){
  leadid<-leads[i]
  leadid<-as.character(leadid)
  leadloc<-sgwas[sgwas$SNP %like% leadid, ]$BP
  chrom<-sgwas[sgwas$SNP %like% leadid,]$CHR
  sgwas.subs<-sgwas[CHR %in% chrom]
  sgwas.subs<-sgwas.subs[BP<(as.integer(leadloc+window)) & BP>(as.integer(leadloc-window))]
  sgwas.subs<-sgwas.subs[!rowSums(is.na(sgwas.subs)) >= 1, ]
  sgwas.subs<-transform(sgwas.subs, A1FREQ=ifelse(A1FREQ > 0.5, 1-A1FREQ, A1FREQ))
  sgwas.subs<-transform(sgwas.subs, SE=(SE^2)/gwas_size)
    for (z in 1:tissue_range){
      ofile=paste(tissues[z], '_', chrom, ".txt", sep="")
      gtexdata<-fread(ofile)
      print(ofile, z)
        if (ncol(gtexdata)==13){
          headlist<-c('gene_id', 'chr', 'position', 'allele1', 'allele2', 'samp', 'tss_distance', 'ma_samples', 'ma_count', 'maf', 'pval_nominal', 'slope', 'slope_se')
          colnames(gtexdata)<-headlist
          gtexdata$newid<-with(gtexdata, paste0(chrom, "_", position, "_", allele1, "_", allele2))
        } else {
          headlist<-c('gene_id', 'variant_id', 'tss_position', 'ma_samples', 'ma_count', 'maf', 'pval_nominal', 'slope', 'slope_se')
          colnames(gtexdata)<-headlist
          gtexdata$variant_id<-as.character(gtexdata$variant_id)
          #splist<-strsplit(gtexdata$variant_id, "_")
          gtexdata$position<-sapply(strsplit(gtexdata$variant_id,"_"), `[`, 2)
          gtexdata$newid<-gsub("(\\S+_\\S+_\\S+_\\S+)_\\S+$", "\\1", gtexdata$variant_id)
        }
      gtexdata<-gtexdata %>% distinct(position, .keep_all = TRUE)
      gtexdata.subs<-gtexdata[position < (as.integer(leadloc+window)) & position>(as.integer(leadloc-window))]
      gtexdata.subs<-gtexdata.subs[!rowSums(is.na(gtexdata.subs)) >= 1, ]
      gtexdata.subs<-gtexdata.subs[order(position),]
      gtexdata.subs<-gtexdata.subs[position %in% sgwas.subs$BP]
      sgwas.subs<-sgwas.subs[BP %in% gtexdata.subs$position]
      sgwas.subs<-sgwas.subs[order(BP),]
      gtexdata.subs<-transform(gtexdata.subs, slope_se=(slope_se^2))
  
      my.res <- coloc.abf(dataset1=list(beta=sgwas.subs$BETA, varbeta=sgwas.subs$SE,
                    N=gwas_size,type="quant", snp=sgwas.subs$namesnp, MAF=sgwas.subs$A1FREQ),
                    dataset2=list(beta=gtexdata.subs$slope, varbeta=gtexdata.subs$slope_se,
                    N=eqtl_size,type="quant", snp=gtexdata.subs$newid, MAF=gtexdata.subs$maf)
                )
      cat(sprintf("\"%f\" \"%s\" \"%s\" \"%s\" \"%s\"", i, leadid, leadloc, chrom, ofile))
      #colocsum<-cbind(my.res$summary, sgwas.subs$SNP)
      colocsum<-my.res$summary
      tissues.vec<-as.vector(tissues[z])
      leadsnp<-as.vector(c(leadid, leadloc))
      names(tissues.vec)<-'tissue'
      names(leadsnp)<-c('Lead_SNP', 'SNP_LOC')
      colocsum<-append(colocsum, tissues.vec)
      colocsum<-append(colocsum, leadsnp)
      #print (colocsum)
          if (i==1){
            reslist<-colocsum
          }
      reslist<-rbind(reslist, colocsum)
    }
}

setwd(file.path(working.dir))
write.csv(reslist, file='colocABFtest-results.csv')

#manhattan plot 
plot(gtexdata$position, -log10(gtexdata$pval_nominal), xlab='chromosome position', ylab='-log10 p-value', pch=19)
