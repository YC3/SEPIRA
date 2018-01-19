context("summary of network")

test_that("the summary information for the network is correct", {

    # model gene expression data set
    data("GeneExp")
    # TFs
    data("TFeid")
    # run the function
    net.o <-  sepiraInfNet(data=GeneExp, tissue=colnames(GeneExp), toi = "Lung", cft = "Blood", TFs = TFeid, sdth = 0.25, sigth = 0.05, pcorth = 0.2, degth = c(0.3, 0.3), lfcth = c(log2(1.5), 0), minNtgts = 5, ncores = 1)

    net.m <- net.o$netTOI
    sum.m <- net.o$sumnet

    # number of TFs
    expect_equal(ncol(net.m), nrow(sum.m))

    # number of targets
    expect_equal(colSums(abs(net.m)), sum.m[, 1])

    # number of avtive interactions
    expect_equal(apply(net.o$netTOI, 2, function(x) {
        i <- which(x == 1)
        return(length(i))
    }), sum.m[, 2])

    # number of negative interactions
    expect_equal(apply(net.o$netTOI, 2, function(x) {
        i <- which(x == -1)
        return(length(i))
    }), sum.m[, 3])

})  

