context("Number of tests")

test_that("the number of tests is correct in LimmaFn", {

    # prepare the data matrix.
    data.m <- matrix(rnorm(120), nrow = 20, ncol = 6)
    # Two different phenotypes
    pheno.2 <- c("C", "C", "C", "T", "T", "T")
    lim.2 <- LimmaFn(pheno.2, data.m)
    # Three different phenotypes
    pheno.3 <- c("C", "C", "N", "T", "T", "T")
    lim.3 <- LimmaFn(pheno.3, data.m)
    # Four different phenotypes
    pheno.4 <- c("C", "C", "N", "T", "T", "Q")
    lim.4 <- LimmaFn(pheno.4, data.m)

    expect_equal(length(lim.2$top), 1)
    expect_equal(length(lim.3$top), 3)
    expect_equal(length(lim.4$top), 6)

})


test_that("the number of comparisons is correct in sepiraInfNet", {
    # model gene expression data set
    data("GeneExp")
    # TFs
    data("TFeid")
    # run the function
    cft <- "Blood"
    net.o <- sepiraInfNet(data=GeneExp, tissue=colnames(GeneExp), toi = "Lung", cft = cft, TFs = TFeid, sdth = 0.25, sigth = 0.05, pcorth = 0.2, degth = c(0.3, 0.3), lfcth = c(log2(1.5), 0), minNtgts = 5, ncores = 1)


    if(is.null(cft)) expect_equal(length(net.o$top), 1)
	if(!is.null(cft)) expect_equal(length(net.o$top), 2)

})

