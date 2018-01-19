context("TF activity matrix")

test_that("the output of 'sepiraRegAct' has the correct dimension", {

    # model gene expression data set
    data("GeneExp")
    # TFs
    data("TFeid")
    # run the function
    net.o <- sepiraInfNet(data=GeneExp, tissue=colnames(GeneExp), toi = "Lung", cft = "Blood", TFs = TFeid, sdth = 0.25, sigth = 0.05, pcorth = 0.2, degth = c(0.3, 0.3), lfcth = c(log2(1.5), 0), minNtgts = 5, ncores = 1)

    # normalized LSCC DNAm data set from TCGA
    data("LUSCmeth")
    # check the activity of lung-specific TFs in inferred network
    TFactLung.v <- sepiraRegAct(GeneExp, type = "exp", regnet = net.o$netTOI, norm = "z", ncores = 1)
    # TF activity in another data set
    TFact <- sepiraRegAct(LUSCmeth, type = "DNAm", regnet = net.o$netTOI, norm = "z", ncores = 1)

    net <- net.o$netTOI
    # check the TFs
    expect_identical(rownames(TFact), colnames(net))
    # check the samples
    expect_identical(colnames(TFact), colnames(LUSCmeth))

})  
