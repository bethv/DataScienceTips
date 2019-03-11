MyPlotEigenNets<-function(multiME){
    p<-plotEigengeneNetworks(
        multiME, 
        setLabels, 
        letterSubPlots = FALSE, Letters = NULL, 
        excludeGrey = TRUE, greyLabel = "grey", 
        plotDendrograms = FALSE, plotHeatmaps = TRUE, 
        setMargins = TRUE, marDendro = NULL,
        colorLabels = TRUE, signed = TRUE, 
        heatmapColors = NULL, 
        plotAdjacency = FALSE,
        printAdjacency = FALSE, cex.adjacency = 0.9,
        coloredBarplot = TRUE, barplotMeans = FALSE, barplotErrors = FALSE, 
        plotPreservation = "standard",
        printPreservation = FALSE, cex.preservation = 0.9,
        marHeatmap = c(3,3,2,1),
        zlimPreservation = c(0, 1),
        xLabelsAngle = 90)
    print(p)
}

view.GOstats<-function(df,ontology){
    dat.mod<-df
    selectedEntrezIds<-as.character(dat.mod$EntrezID)
    params <- new("GOHyperGParams", geneIds = selectedEntrezIds, universeGeneIds = allLLIDs, annotation = "hgu95av2.db", ontology = ontology, pvalueCutoff = 0.05, conditional = T, testDirection = "over")
    hgOver <- hyperGTest(params)
    table <- summary(hgOver)
    return(table)
}


makeGeneInfo<-function(geneData,modules,geneModMemb){
    geneInfo0 = data.frame(ProbeID = geneData$Probe, 
                           geneSymbol = geneData$Symbol,
                           EntrezID = geneData$EntrezID,
                           Name=geneData$Name,
                           moduleColor = modules)
    Modmemb<-numeric()
    for(i in 1:nrow(geneInfo0)){
        modColor<-modules[i]
        Modmemb[i]<-geneModMemb[i,modColor]
    }
    geneInfo0$MM<-Modmemb
    geneOrder = order(factor(geneInfo0$moduleColor,levels = names(geneModMemb)), -abs(geneInfo0$MM))
    
    geneInfo<-geneInfo0[geneOrder,]
    geneInfo<-geneInfo[geneInfo$geneSymbol!="",]
    geneInfo<-geneInfo[!is.na(geneInfo$EntrezID),]
    geneInfo$absMM<-abs(geneInfo$MM)
    return(geneInfo)
}
