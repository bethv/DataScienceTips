



# Intro

##### Basic process for RNA seq data

- counts downloaded from recount project
- scale counts
- DeSeq
- limit to genes with counts > 10 in at least 30%
- vst (DeSeq2)
- normalizeQuantiles (limma)


#####Recount website info
main website: <https://jhubiostatistics.shinyapps.io/recount/>
examples: <http://leekgroup.github.io/recount-analyses/>

####WGCNA website

**Can WGCNA be used to analyze RNA-Seq data?(summarized and reproduced from WCGNA website)**

 - remove low counts

 - variance-stabilizing transformation: "For example, package DESeq2 implements the function varianceStabilizingTransformation which we have found useful, but one could also start with normalized counts (or RPKM/FPKM data) and log-transform them using log2(x+1). For highly expressed features, the differences between full variance stabilization and a simple log transformation are small." 

 - normalization: "Finally, we usually check quantile scatterplots to make sure there are no systematic shifts between samples; if sample quantiles show correlations (which they usually do), quantile normalization can be used to remove this effect."

   

# Creation of expression matrices

```R
#source('http://bioconductor.org/biocLite.R')
#biocLite('recount')

library('recount')
library('DESeq2')
library(limma)
```

#### Colon

```R
load("rse_geneSig34.rda")
rse<-rse_geneSig34
rm(rse_geneSig34)
rseSig <- scale_counts(rse_geneSig34)
```

##### scale counts with recount package

```R
rse <- scale_counts(rse)
```

##### DESeq, controlling for RIN and batch

```R
dds <- DESeqDataSet(rse, design = ~as.numeric(smrin)+smgebtch)
dds <- DESeq(dds)
```

##### filter to counts>10 in 30%

```R
n30percent<-.3*dim(rse)[2]
countsGr10in30per<-apply(counts(dds),1,function(x){
  x2<-ifelse(x>10,1,0)
  x3<-ifelse(sum(x2)>n30percent,TRUE,FALSE)
  return(x3)
})

dds <- dds[ countsGr10in30per, ]
nrow(dds)

library(stats)
# smgebtch
rse1 <- assay(rse)
pc1 <- princomp(rse1)
ba1 <- rse$smgebtch
```

##### vsd

```R
vsd <- vst(dds, blind = FALSE)
edata=as.data.frame(assay(vsd))
```

##### scatter plot and normalization


```R
colramp = colorRampPalette(c(3,"white",2))(20)
plot(density(edataSig[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:73){lines(density(edataSig[,i]),lwd=3,col=colramp[i])}

norm_edataSig = normalizeQuantiles(as.matrix(edataSig))
plot(density(norm_edataSig[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:73){lines(density(norm_edataSig [,i]),lwd=3,col=colramp[i])}

datExprSig<-as.data.frame(t(norm_edataSig))
rownames(datExprSig)<-rse_geneSig34$Person

# save(ddsSig,file = "ddsSig.rda")
# save(vsdSig,file = "vsdSig.rda")
# save(edataSig,file = "edataSig.rda")
# save(norm_edataSig,file = "norm_edataSig.rda")
# save(datExprSig,file = "datExprSig.rda")
```

### brain 
```R
cortFilt<-rse_geneCortex$Person%in%samps # limit to samples with sigmoid data
rse_geneCort34<-rse_geneCortex[,cortFilt]
dim(rse_geneCort34) #[1] 58037    34
save(rse_geneCort34,file = "rse_geneCort34.rda")

load("rse_geneCort34.rda")
rseCort34<-scale_counts(rse_geneCort34)
rm(rse_geneCort34)
ddsCort<-DESeqDataSet(rseCort34,~as.numeric(smrin))
ddsCort<-DESeq(ddsCort)

nrow(ddsCort)
ncol(ddsCort)

n30percent<-.3*dim(rseCort34)[2]
countsGr10in30per<-apply(counts(ddsCort),1,function(x){
  x2<-ifelse(x>10,1,0)
  x3<-ifelse(sum(x2)>n30percent,TRUE,FALSE)
  return(x3)
})

ddsCort <- ddsCort[ countsGr10in30per, ]
nrow(ddsCort) #24048
vsdCort <- vst(ddsCort, blind = FALSE)
edataCort=as.data.frame(assay(vsdCort))
colramp = colorRampPalette(c(3,"white",2))(20)

plot(density(edataCort[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:33){lines(density(edataCort[,i]),lwd=3,col=colramp[i])}

norm_edataCort = normalizeQuantiles(as.matrix(edataCort))
plot(density(norm_edataCort[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:33){lines(density(norm_edataCort[,i]),lwd=3,col=colramp[i])}

datExprCort<-as.data.frame(t(norm_edataCort))
rownames(datExprCort)<-rseCort34$Person

# save(ddsCort,file = "ddsCort.rda")
# save(vsdCort,file = "vsdCort.rda")
# save(edataCort,file = "edataCort.rda")
# save(norm_edataCort,file = "norm_edataCort.rda")
# save(datExprCort,file = "datExprCort.rda")
```

### get rid of X and Y
```R
load("genedataAllcolon.rda") # probes

ensemblIDs <- sapply( strsplit( as.character(genedataAll$Probe), split="\\." ), "[", 1 )

library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description","chromosome_name"),
                  filters = "ensembl_gene_id",
                  values = ensemblIDs,
                  mart = ensembl )
idx <- match( ensemblIDs, genemap$ensembl_gene_id )
genedataAll<-data.frame(Probe=genedataAll$Probe, ENS=ensemblIDs,EntrezID=genemap$entrezgene[idx],Symbol=genemap$hgnc_symbol[idx], Name=genemap$description[idx],Chr=genemap$chromosome_name[idx])

# save(genedataAll, file = "genedataAllcolon.rda")

noXY<-as.character(genedataAll$Probe[!genedataAll$Chr%in%c("X","Y")])

load("ConsensusSample.rda")
load("datExprSig.rda")

datExprSigNoXY<-datExprSig[,names(datExprSig)%in%noXY]
datExprSigNoXYconsSamp<-datExprSigNoXY[ConsensusSample,]
# save(datExprSigNoXY,file = "datExprSigNoXY.rda")
# save(datExprSigNoXYconsSamp,file = "datExprSigNoXYconsSamp.rda")
# save(noXY,file = "noXY.rda")

load("datExprCort.rda")
datExprCortNoXY<-datExprCort[,names(datExprCort)%in%noXY]
datExprCortNoXYconsSamp<-datExprCortNoXY[ConsensusSample,]
save(datExprCortNoXY,file = "datExprCortNoXY.rda")
save(datExprCortNoXYconsSamp,file = "datExprCortNoXYconsSamp.rda")

Genes<-intersect(names(datExprCortNoXYconsSamp),names(datExprSigNoXYconsSamp))
datExprBrain<-datExprCortNoXYconsSamp[,Genes]
datExprColon<-datExprSigNoXYconsSamp[,Genes]
multiExpr<-multiData(Colon=datExprColon, Brain=datExprBrain)
# save(multiExpr,file = "multiExpr.rda")
```
