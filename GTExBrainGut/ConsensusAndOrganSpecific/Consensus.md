# Consensus network
Of the 34 samples with sigmoid colon and cerebral cortex there are 29 that are not outliers in either (outliers based on sample clustering in WGCNA)
Also determined that softpower 9 worked for both brain and colon. Switched to signed and doubled soft power to 18 per WGCNA site.
```R
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads() 
```

```R
sp = 18
ms = 80
ds = 1
mch = 0.15

net = blockwiseConsensusModules(
  multiExpr, 
  power = sp, 
  minModuleSize = ms, 
  deepSplit = ds,
  pamRespectsDendro = FALSE, 
  networkType = "signed",
  mergeCutHeight = mch, 
  numericLabels = TRUE,
  minKMEtoStay = 0.2, 
  verbose = 5,
  maxBlockSize = 20000)
setLabels = c("Sigmoid Colon", "Cerebral Cortex")
shortLabels = c("Colon", "Brain")
consMEs = net$multiMEs
consModuleLabels = net$colors;
consModuleColors = labels2colors(consModuleLabels)
consMEsC = multiSetMEs(multiExpr, universalColors = consModuleColors)
mods<-substring(names(consMEsC$Brain$data),3)
mods<-mods[!mods%in%"grey"]
consMEs.ConsOrder<-consensusOrderMEs(consMEsC)
```

#### Gene info and GO

```{r}
load("Sp18_Ms80_Ds1_MCH.15/ConsNet_sp18ms80ds1mch0.15_15mods_midnightblue.rda")
load("genedataAllcolon.rda")
source("gtexFunctions.R")
```

```{r}
# gene info consensus colon-code----
geneModuleMembership = as.data.frame(cor(multiExpr$Colon$data, consMEs.ConsOrder$Colon$data, use = "p"))
names(geneModuleMembership) = substring(names(consMEs.ConsOrder$Colon$data),3)
genedataWGCNA<-genedataAll[genedataAll$Probe%in%names(multiExpr$Colon$data),]
geneInfo.Ccons<-makeGeneInfo(genedataWGCNA,consModuleColors,geneModuleMembership)

write.csv(as.data.frame(table(factor(consModuleColors,levels = names(geneModuleMembership)))),paste(fname.root,"_ModSizes.csv",sep = ""))


# gene info consensus brain-code----
geneModuleMembership = as.data.frame(cor(multiExpr$Brain$data, consMEs.ConsOrder$Brain$data, use = "p"))
names(geneModuleMembership) = substring(names(consMEs.ConsOrder$Brain$data),3)
genedataWGCNA<-genedataAll[genedataAll$Probe%in%names(multiExpr$Brain$data),]
geneInfo.Bcons<-makeGeneInfo(genedataWGCNA,consModuleColors,geneModuleMembership)

# both MM columns
geneInfo.consBC<-geneInfo.Bcons[order(geneInfo.Bcons$ProbeID),]
names(geneInfo.consBC)[6]<-"MM.B"
geneInfo.Ccons<-geneInfo.Ccons[order(geneInfo.Ccons$ProbeID),]
geneInfo.consBC$MM.C<-geneInfo.Ccons$MM
geneInfo.consBC$absMM.B<-abs(geneInfo.consBC$MM.B)
geneInfo.consBC$absMM.C<-abs(geneInfo.consBC$MM.C)
geneInfo.consBC$avgMM<-rowMeans(geneInfo.consBC[,c("absMM.B","absMM.C")])

geneInfo.consBC<-geneInfo.consBC[(geneInfo.consBC$absMM.B>=0.6)|(geneInfo.consBC$absMM.C>=0.6),]

geneOrder = order(factor(geneInfo.consBC$moduleColor,levels = names(geneModuleMembership)), -abs(geneInfo.consBC$MM.C))
geneInfo.consBC<-geneInfo.consBC[geneOrder,]

geneInfo.consBC_all<-geneInfo.consBC

write.csv(geneInfo.consBC, paste(fname.root,"_INFO.csv",sep = ""))
save(geneInfo.consBC,file=paste(fname.root,"_INFO.rda",sep = ""))
```

##### GO - either MM>=0.6

```{r}
library(GOstats)
library(hgu95av2.db)

genedatacons<-genedataAll[genedataAll$Probe%in%names(multiExpr$Brain$data),]
allLLIDscons <- genedatacons$EntrezID
allLLIDscons<allLLIDscons[!allLLIDscons%in%NA]
allLLIDscons<-as.character(allLLIDscons)
geneInfo.consBC<-geneInfo.consBC[geneInfo.consBC$moduleColor!="grey",]
modListBC<-split(geneInfo.consBC,geneInfo.consBC$moduleColor)

allLLIDs<-allLLIDscons
modList<-modListBC

GOlist<-list()

for (module in names(modList)) {
  df1<-view.GOstats(modList[[module]],"BP")
  df1$moduleColor<-module
  df1$moduleSize<-nrow(modList[[module]])
  df1<-df1[order(df1$Pvalue),]
  GOlist[[module]]<-df1
}

save(GOlist,file = paste(fname.root,"_GOlist.rda",sep = ""))

for (module in names(GOlist)) {
  GOlist[[module]]
  GOlist[[module]]$TermPval<-paste(GOlist[[module]]$Term," (",format(GOlist[[module]]$Pvalue,digits = 2,scientific = TRUE),")",sep = "")
}

GOtable<-GOlist$blue[0,]
for (module in names(GOlist)) {
  GOtable<-rbind(GOtable,GOlist[[module]])
}
GOtable<-GOtable[GOtable$Size>=10,]
write.csv(GOtable,paste(fname.root,"_GOtable.csv",sep = ""))
```

##### GO- both 0.6

```{r}
load("Sp18_Ms80_Ds1_MCH.15/ConsNet_sp18ms80ds1mch0.15_15mods_midnightblue_INFO.rda")

library(GOstats)
library(hgu95av2.db)

genedatacons<-genedataAll[genedataAll$Probe%in%names(multiExpr$Brain$data),]
allLLIDscons <- genedatacons$EntrezID
allLLIDscons<-allLLIDscons[!allLLIDscons%in%NA]
allLLIDscons<-as.character(allLLIDscons)
geneInfo.consBC<-geneInfo.consBC[geneInfo.consBC$moduleColor!="grey",]

BandCgeneInfo.conBC<-geneInfo.consBC[geneInfo.consBC$absMM.B>=0.6&geneInfo.consBC$absMM.C>=0.6,]
BgeneInfo.conBC<-geneInfo.consBC[geneInfo.consBC$absMM.B>=0.6,]
CgeneInfo.conBC<-geneInfo.consBC[geneInfo.consBC$absMM.C>=0.6,]

BandCmodListBC<-split(BandCgeneInfo.conBC,BandCgeneInfo.conBC$moduleColor)
BmodListBC<-split(BgeneInfo.conBC,BgeneInfo.conBC$moduleColor)
CmodListBC<-split(CgeneInfo.conBC,CgeneInfo.conBC$moduleColor)

allLLIDs<-allLLIDscons

modList<-BandCmodListBC
fname.root<-"BandCgenes"
source("RunGO.R")

modList<-BmodListBC
fname.root<-"Bgenes"
source("RunGO.R")

modList<-CmodListBC
fname.root<-"Cgenes"
source("RunGO.R")

save(geneInfo.consBC, BandCgeneInfo.conBC, BgeneInfo.conBC, CgeneInfo.conBC, file = "geneInfoDFsConsBC.rda")

```

### Lung neg control


```R
load("rse_gene_lung.Rdata")
rse_gene$person<-gsub("GTEX-","",rse_gene$sampid)
rse_gene$person<-gsub("-.+","",rse_gene$person)

bcfilt<-rse_gene$person%in%subs
rse_gene2<-rse_gene[,bcfilt]
dim(rse_gene2) #19

rseLung <- scale_counts(rse_gene2)
# save(rseLung,file = "rseLung.rda")

library('DESeq2')

dds <- DESeqDataSet(rseLung, ~as.numeric(smrin))
dds <- DESeq(dds)
nrow(dds)
ncol(dds)

n30percent<-.3*dim(rseLung)[2]
countsGr10in30per<-apply(counts(dds),1,function(x){
        x2<-ifelse(x>10,1,0)
        x3<-ifelse(sum(x2)>n30percent,TRUE,FALSE)
        return(x3)
})

dds <- dds[ countsGr10in30per, ]
nrow(dds) #25031
vsd <- vst(dds, blind = FALSE)
edata=as.data.frame(assay(vsd))
colramp = colorRampPalette(c(3,"white",2))(20)
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:73){lines(density(edata[,i]),lwd=3,col=colramp[i])}

library(limma)
norm_edata = normalizeQuantiles(as.matrix(edata))
plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:73){lines(density(norm_edata [,i]),lwd=3,col=colramp[i])}

datExpr<-as.data.frame(t(norm_edata))
rownames(datExpr)<-rseLung$person

# save(dds,file = "ddsLung.rda")
# save(vsd,file = "vsdLung.rda")
# save(edata,file = "edataLung.rda")
# save(norm_edata,file = "norm_edataLung.rda")
# save(datExpr,file = "datExprLung.rda")

ensemblIDs <- sapply( strsplit( rowData(rseLung)$gene_id, split="\\." ), "[", 1 )

library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description","chromosome_name"),
                  filters = "ensembl_gene_id",
                  values = ensemblIDs,
                  mart = ensembl )
idx <- match( ensemblIDs, genemap$ensembl_gene_id )
load("genedataAllcolon.rda")

genedataAll<-data.frame(Probe=genedataAll$Probe, ENS=ensemblIDs,EntrezID=genemap$entrezgene[idx],Symbol=genemap$hgnc_symbol[idx], Name=genemap$description[idx],Chr=genemap$chromosome_name[idx])

# save(genedataAll, file = "genedataAllLung.rda")

noXY<-as.character(genedataAll$Probe[!genedataAll$Chr%in%c("X","Y")])

datExprNoXY<-datExpr[,names(datExpr)%in%noXY]
# save(datExprNoXY,file = "datExprNoXYLung.rda")
# save(noXY,file = "noXYLung.rda")


load("datExprCortNoXYconsSamp.rda")
load("datExprNoXYLung.rda")
lung<-datExprNoXY
brain<-datExprCortNoXYconsSamp[rownames(datExprNoXY),]
rownames(brain);rownames(lung)
Genes<-intersect(names(brain),names(lung))
datExprBrain<-brain[,Genes]
datExprLung<-lung[,Genes]
multiExpr<-multiData(Lung=datExprLung, Brain=datExprBrain)
save(multiExpr,file = "multiExprBrainLung.rda")

source("gtexFunctions.R")
load("multiExprBrainLung.rda")


sp = 18
ms = 80
ds = 1
mch = 0.15
Fname<-paste("BL_sp",sp,"ms",ms,"ds",ds,"mch",mch,"_",sep = "")

net = blockwiseConsensusModules(
  multiExpr, 
  power = sp, 
  minModuleSize = ms, 
  deepSplit = ds,
  pamRespectsDendro = FALSE, 
  networkType = "signed",
  mergeCutHeight = mch, 
  numericLabels = TRUE,
  minKMEtoStay = 0.2, 
  verbose = 5,
  maxBlockSize = 20000)
setLabels = c("Lung", "Cerebral Cortex")
shortLabels = c("Lung", "Brain")
consMEs = net$multiMEs
consModuleLabels = net$colors;
consModuleColors = labels2colors(consModuleLabels)
consMEsC = multiSetMEs(multiExpr, universalColors = consModuleColors)
mods<-substring(names(consMEsC$Brain$data),3)
mods<-mods[!mods%in%"grey"]
consMEs.ConsOrder<-consensusOrderMEs(consMEsC)

# save(multiExpr,consMEs,consModuleLabels,consModuleColors,consMEsC,net,setLabels,shortLabels,consMEs.ConsOrder,mods,file=paste("ConsNet",Fname,".rda",sep = ""))

pdf(paste("EigNet",Fname,".pdf",sep = ""))
    MyPlotEigenNets(consMEs.ConsOrder)
dev.off()


# gene info --------------
fname.root<-"ConsNetLB"
# consensus gene info------------
load("ConsNetBL_sp18ms80ds1mch0.15_.rda")
load("genedataAllcolon.rda")
source("gtexFunctions.R")
# gene info consensus colon-code----
geneModuleMembership = as.data.frame(cor(multiExpr$Lung$data, consMEs.ConsOrder$Lung$data, use = "p"))
names(geneModuleMembership) = substring(names(consMEs.ConsOrder$Lung$data),3)
genedataWGCNA<-genedataAll[genedataAll$Probe%in%names(multiExpr$Lung$data),]
geneInfo.Lcons<-makeGeneInfo(genedataWGCNA,consModuleColors,geneModuleMembership)

write.csv(as.data.frame(table(factor(consModuleColors,levels = names(geneModuleMembership)))),paste(fname.root,"_ModSizes.csv",sep = ""))


# gene info consensus brain-code----
geneModuleMembership = as.data.frame(cor(multiExpr$Brain$data, consMEs.ConsOrder$Brain$data, use = "p"))
names(geneModuleMembership) = substring(names(consMEs.ConsOrder$Brain$data),3)
genedataWGCNA<-genedataAll[genedataAll$Probe%in%names(multiExpr$Brain$data),]
geneInfo.Bcons<-makeGeneInfo(genedataWGCNA,consModuleColors,geneModuleMembership)


# gene info consensus lung-code----
geneModuleMembership = as.data.frame(cor(multiExpr$Lung$data, consMEs.ConsOrder$Lung$data, use = "p"))
names(geneModuleMembership) = substring(names(consMEs.ConsOrder$Lung$data),3)
genedataWGCNA<-genedataAll[genedataAll$Probe%in%names(multiExpr$Lung$data),]
geneInfo.Lcons<-makeGeneInfo(genedataWGCNA,consModuleColors,geneModuleMembership)

# all MM columns
geneInfo.consLB<-geneInfo.Lcons[order(geneInfo.Lcons$ProbeID),]
names(geneInfo.consLB)[6]<-"MM.L"
geneInfo.Bcons<-geneInfo.Bcons[order(geneInfo.Bcons$ProbeID),]
geneInfo.consLB$MM.B<-geneInfo.Bcons$MM

geneInfo.consLB$absMM.L<-abs(geneInfo.consLB$MM.L)
geneInfo.consLB$absMM.B<-abs(geneInfo.consLB$MM.B)

geneInfo.consLB$avgMM<-rowMeans(geneInfo.consLB[,c("absMM.B", "absMM.L")])
geneInfo.consLB<-geneInfo.consLB[(geneInfo.consLB$absMM.B>=0.6)|(geneInfo.consLB$absMM.L>=0.6),]

geneOrder = order(factor(geneInfo.consLB$moduleColor,levels = names(geneModuleMembership)), -abs(geneInfo.consLB$MM.L))
geneInfo.consLB<-geneInfo.consLB[geneOrder,]

# write.csv(geneInfo.consLB, paste(fname.root,"_INFO.csv",sep = ""))
# save(geneInfo.consLB,file=paste(fname.root,"_INFO.rda",sep = ""))

#note: snca s100b sox10 are in grey for LB and LCB


# GO ------------
library(GOstats)
library(hgu95av2.db)
genedatacons<-genedataAll[genedataAll$Probe%in%names(multiExpr$Brain$data),]
allLLIDscons <- genedatacons$EntrezID
allLLIDscons<-allLLIDscons[!allLLIDscons%in%NA]
allLLIDscons<-as.character(allLLIDscons)
save(allLLIDscons, file = "allLLIDscons.rda")
geneInfo.consLB<-geneInfo.consLB[geneInfo.consLB$moduleColor!="grey",]
modListLB<-split(geneInfo.consLB,geneInfo.consLB$moduleColor)

allLLIDs<-allLLIDscons
modList<-modListLB

GOlist<-list()

for (module in names(modList)) {
  df1<-view.GOstats(modList[[module]],"BP")
  df1$moduleColor<-module
  df1$moduleSize<-nrow(modList[[module]])
  df1<-df1[order(df1$Pvalue),]
  GOlist[[module]]<-df1
}

# save(GOlist,file = paste(fname.root,"_GOlist.rda",sep = ""))

for (module in names(GOlist)) {
  GOlist[[module]]
  GOlist[[module]]$TermPval<-paste(GOlist[[module]]$Term," (",format(GOlist[[module]]$Pvalue,digits = 2,scientific = TRUE),")",sep = "")
}

GOtable<-GOlist$blue[0,]
for (module in names(GOlist)) {
  GOtable<-rbind(GOtable,GOlist[[module]])
}
GOtable<-GOtable[GOtable$Size>=10,]
write.csv(GOtable,paste(fname.root,"_GOtable.csv",sep = ""))
```
