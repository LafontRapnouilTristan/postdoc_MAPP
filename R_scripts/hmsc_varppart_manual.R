Y <- Y
XFormula = as.formula(XFormula_list[1])
XData = XData
X = NULL
XScale = TRUE
XSelect = NULL
XRRRData = NULL
XRRRFormula = ~. - 1
XRRR = NULL
ncRRR = 2
XRRRScale = TRUE
YScale = FALSE
Loff = NULL
studyDesign = studyDesign
ranLevels = list(Sample = rL.spatial)
ranLevelsUsed = names(ranLevels)
TrFormula = NULL
TrData = NULL
Tr = NULL
TrScale = TRUE
phyloTree = NULL
C = NULL
distr = "lognormal poisson"
truncateNumberOfFactors = TRUE

hM = structure(list(Y = NULL, Loff = NULL, XData = NULL,
        XFormula = NULL, X = NULL, XScaled = NULL, XSelect = NULL,
        XRRRData = NULL, XRRRFormula = NULL, XRRR = NULL, XRRRScaled = NULL,
        YScaled = NULL, XInterceptInd = NULL, studyDesign = NULL,
        ranLevels = NULL, ranLevelsUsed = NULL, dfPi = NULL,
        rL = NULL, Pi = NULL, TrData = NULL, TrFormula = NULL,
        Tr = NULL, TrScaled = NULL, TrInterceptInd = NULL, C = NULL,
        phyloTree = NULL, distr = NULL, ny = NULL, ns = NULL,
        nc = NULL, ncNRRR = NULL, ncRRR = NULL, ncORRR = NULL,
        ncsel = NULL, nr = NULL, nt = NULL, nf = NULL, ncr = NULL,
        ncs = NULL, np = NULL, spNames = NULL, covNames = NULL,
        trNames = NULL, rLNames = NULL, XScalePar = NULL, XRRRScalePar = NULL,
        YScalePar = NULL, TrScalePar = NULL, V0 = NULL, f0 = NULL,
        mGamma = NULL, UGamma = NULL, aSigma = NULL, bSigma = NULL,
        rhopw = NULL, nuRRR = NULL, a1RRR = NULL, b1RRR = NULL,
        a2RRR = NULL, b2RRR = NULL, samples = NULL, transient = NULL,
        thin = NULL, verbose = NULL, adaptNf = NULL, initPar = NULL,
        repN = NULL, randSeed = NULL, postList = NULL), class = "Hmsc")
    Y <- as.matrix(Y)
    if (!(is.numeric(Y) || is.logical(Y)) && is.finite(sum(Y, 
                                                           na.rm = TRUE))) {
        stop("Y must be a numeric matrix of finite values of sampling units times species")
    }
    hM$Y = Y
    hM$ny = nrow(Y)
    hM$ns = ncol(Y)
    if (is.null(colnames(hM$Y))) {
        colnames(hM$Y) = sprintf(sprintf("sp%%.%dd", ceiling(log10(hM$ns))), 
                                 1:hM$ns)
    }
    hM$spNames = colnames(hM$Y)
    if (!is.null(X)) {
        if (!is.null(XData)) 
            stop("only one of XData and X arguments can be specified")
        if (!is.null(XFormula)) 
            XFormula <- NULL
    }
    if (!is.null(XData)) {
        BADvars <- function(DF) {
            cl <- vapply(DF, .MFclass, "")
            ok <- cl %in% c("numeric", "factor", "ordered", 
                            "logical")
            cl <- cl[!ok]
            if (length(cl) > 0) 
                paste(names(cl), cl, sep = ":", collapse = ", ")
            else invisible(NULL)
        }
        if (inherits(XData, "list")) {
            if (length(XData) != hM$ns) {
                stop("the length of XData list argument must be equal to the number of species")
            }
            if (any(!unlist(lapply(XData, is.data.frame)))) {
                stop("each element of X list must be a data.frame")
            }
            if (any(sapply(XData, function(a) class(a)[1L] != 
                           "data.frame"))) {
                for (i in seq_len(length(XData))) XData[[i]] <- as.data.frame(XData, 
                                                                              stringsAsFactors = TRUE)
            }
            if (any(unlist(lapply(XData, function(a) nrow(a) != 
                                  hM$ny)))) {
                stop("for each element of XData list the number of rows must be equal to the number of sampling units")
            }
            if (any(!sapply(lapply(XData, BADvars), is.null))) 
                stop("XData variables had bad types (e.g., matrix, non-numeric)")
            if (any(unlist(lapply(XData, function(a) any(is.na(a)))))) {
                stop("NA values are not allowed in XData")
            }
            hM$XData = XData
            hM$XFormula = XFormula
            hM$X = lapply(XData, function(a) model.matrix(XFormula, 
                                                          a))
            hM$nc = ncol(hM$X[[1]])
        }
        else if (is.data.frame(XData)) {
            if (nrow(XData) != hM$ny) {
                stop("the number of rows in XData must be equal to the number of sampling units")
            }
            if (any(is.na(XData))) {
                stop("NA values are not allowed in XData")
            }
            if (!is.null(baddie <- BADvars(XData))) 
                stop("XData variables had bad types: ", baddie)
            if (class(XData)[1L] != "data.frame") 
                XData <- as.data.frame(XData, stringsAsFactors = TRUE)
            hM$XData = XData
            hM$XFormula = XFormula
            hM$X = model.matrix(XFormula, XData)
            hM$nc = ncol(hM$X)
        }
        else {
            stop("XData must be either a data.frame or a list of data.frame objects")
        }
    }
    if (!is.null(X)) {
        switch(class(X)[1L], list = {
            if (length(X) != hM$ns) {
                stop("the length of X list argument must be equal to the number of species")
            }
            if (any(!unlist(lapply(X, is.matrix)))) {
                stop("each element of X list must be a matrix")
            }
            if (any(unlist(lapply(X, function(a) nrow(a) != 
                                  hM$ny)))) {
                stop("for each element of X list the number of rows must be equal to the number of sampling units")
            }
            if (any(unlist(lapply(X, function(a) any(is.na(a)))))) {
                stop("NA values are not allowed in X")
            }
            hM$XData = NULL
            hM$XFormula = NULL
            hM$X = X
            hM$nc = ncol(hM$X[[1]])
        }, matrix = {
            if (nrow(X) != hM$ny) {
                stop("the number of rows in X must be equal to the number of sampling units")
            }
            if (any(is.na(X))) {
                stop("NA values are not allowed in X")
            }
            hM$XData = NULL
            hM$XFormula = NULL
            hM$X = X
            hM$nc = ncol(hM$X)
        }, {
            stop("X must be a matrix or a list of matrix objects")
        })
    }
    if (is.null(XData) && is.null(X)) {
        X = matrix(NA, hM$ny, 0)
        hM$nc = 0
    }
    switch(class(hM$X)[1L], matrix = {
        if (is.null(colnames(hM$X))) {
            colnames(hM$X) = sprintf(sprintf("cov%%.%dd", ceiling(log10(hM$nc))), 
                                     1:hM$nc)
        }
    }, if (is.null(colnames(hM$X[[1]]))) {
        list = {
            for (j in 1:hM$ns) colnames(hM$X[[j]]) = sprintf(sprintf("cov%%.%dd", 
                                                                     ceiling(log10(hM$nc))), 1:hM$nc)
        }
    })
    switch(class(hM$X)[1L], matrix = {
        hM$covNames = colnames(hM$X)
    }, list = {
        hM$covNames = colnames(hM$X[[1]])
    })
    if (!XScale) {
        hM$XScalePar = rbind(rep(0, hM$nc), rep(1, hM$nc))
        hM$XScaled = hM$X
        hM$XInterceptInd = NULL
    }
    else {
        switch(class(hM$X)[1L], matrix = {
            XStack = hM$X
        }, list = {
            XStack = Reduce(rbind, hM$X)
        })
        XInterceptInd = which(colnames(XStack) %in% c("Intercept", 
                                                      "(Intercept)"))
        if (length(XInterceptInd) > 1) {
            stop("only one column of X matrix can be named Intercept or (Intercept)")
        }
        if (!all(XStack[, XInterceptInd] == 1)) {
            stop("intercept column in X matrix must be a column of ones")
        }
        if (length(XInterceptInd) == 1) {
            hM$XInterceptInd = XInterceptInd
        }
        else hM$XInterceptInd = NULL
        XScalePar = rbind(rep(0, hM$nc), rep(1, hM$nc))
        XScaled = XStack
        if (XScale) {
            scaleInd = apply(XStack, 2, function(a) !all(a %in% 
                                                             c(0, 1)))
        }
        else {
            scaleInd = XScale
        }
        scaleInd[XInterceptInd] = FALSE
        if (length(XInterceptInd) > 0) {
            sc = scale(XStack)
            XScalePar[, scaleInd] = rbind(attr(sc, "scaled:center"), 
                                          attr(sc, "scaled:scale"))[, scaleInd]
        }
        else {
            sc = scale(XStack, center = FALSE)
            XScalePar[2, scaleInd] = attr(sc, "scaled:scale")[scaleInd]
        }
        XScaled[, scaleInd] = sc[, scaleInd]
        hM$XScalePar = XScalePar
        switch(class(hM$X)[1L], matrix = {
            hM$XScaled = XScaled
        }, list = {
            hM$XScaled = lapply(split(XScaled, rep(1:hM$ns, 
                                                   each = hM$ny)), function(a) matrix(a, hM$ny, 
                                                                                      hM$nc))
        })
    }
    hM$ncsel = length(XSelect)
    hM$XSelect = XSelect
    for (i in seq_len(hM$ncsel)) {
        XSel = hM$XSelect[[i]]
        if (max(XSel$covGroup) > hM$nc) {
            stop("covGroup for XSelect cannot have values greater than number of columns in X")
        }
    }
    hM$ncNRRR = hM$nc
    if (!is.null(XRRRData)) {
        if (!is.data.frame(XRRRData)) {
            stop("XRRRData must be a data.frame")
        }
        if (nrow(XRRRData) != hM$ny) {
            stop("the number of rows in XRRRData must be equal to the number of sampling units")
        }
        if (any(is.na(XRRRData))) {
            stop("XRRRData must contain no NA values")
        }
        hM$XRRRData = XRRRData
        hM$XRRRFormula = XRRRFormula
        hM$XRRR = model.matrix(XRRRFormula, XRRRData)
        hM$ncORRR = ncol(hM$XRRR)
        hM$ncRRR = ncRRR
    }
    if (!is.null(XRRR)) {
        if (!is.matrix(XRRR)) {
            stop("XRRR must be a matrix")
        }
        if (nrow(XRRR) != hM$ny) {
            stop("the number of rows in XRRR must be equal to the number of sampling units")
        }
        if (any(is.na(XRRR))) {
            stop("XRRR must contain no NA values")
        }
        hM$XRRRData = NULL
        hM$XRRRFormula = NULL
        hM$XRRR = XRRR
        hM$ncORRR = ncol(hM$XRRR)
        hM$ncRRR = ncRRR
    }
    if (is.null(XRRRData) && is.null(XRRR)) {
        X = matrix(NA, hM$ny, 0)
        hM$XRRR = NULL
        hM$ncORRR = 0
        hM$ncRRR = 0
    }
    if (hM$ncRRR > 0) {
        if (is.null(colnames(hM$XRRR))) {
            colnames(hM$XRRR) = sprintf(sprintf("covRRR%%.%dd", 
                                                ceiling(log10(hM$ncORRR))), 1:hM$ncORRR)
        }
        for (k in seq_len(hM$ncRRR)) {
            hM$covNames = c(hM$covNames, paste0("XRRR_", as.character(k)))
        }
        hM$nc = hM$ncNRRR + hM$ncRRR
        if (!XRRRScale) {
            hM$XRRRScalePar = rbind(rep(0, hM$ncORRR), rep(1, 
                                                           hM$ncORRR))
            hM$XRRRScaled = hM$XRRR
        }
        else {
            if (!XScale) {
                stop("XRRR cannot be scaled if X is not scaled")
            }
            XRRRStack = hM$XRRR
            if (XRRRScale) {
                XRRRscaleInd = apply(XRRRStack, 2, function(a) !all(a %in% 
                                                                        c(0, 1)))
            }
            else {
                XRRRscaleInd = XRRRScale
            }
            XRRRScalePar = rbind(rep(0, hM$ncORRR), rep(1, hM$ncORRR))
            XRRRScaled = XRRRStack
            if (length(XInterceptInd) > 0) {
                XRRRsc = scale(XRRRStack)
                XRRRScalePar[, XRRRscaleInd] = rbind(attr(XRRRsc, 
                                                          "scaled:center"), attr(XRRRsc, "scaled:scale"))[, 
                                                                                                          XRRRscaleInd]
            }
            else {
                XRRRsc = scale(XRRRStack, center = FALSE)
                XRRRScalePar[2, XRRRscaleInd] = attr(sc, "scaled:scale")[XRRRscaleInd]
            }
            XRRRScaled[, XRRRscaleInd] = XRRRsc[, XRRRscaleInd]
            hM$XRRRScalePar = XRRRScalePar
            hM$XRRRScaled = XRRRScaled
        }
    }
    if (!is.null(TrData)) {
        if (!is.null(Tr)) {
            stop("only one of TrData and Tr arguments can be specified")
        }
        if (is.null(TrFormula)) {
            stop("TrFormula argument must be specified if TrData is provided")
        }
    }
    if (!is.null(TrData)) {
        if (nrow(TrData) != hM$ns) {
            stop("the number of rows in TrData should be equal to number of columns in Y")
        }
        if (!all(rownames(TrData) == colnames(Y))) {
            stop("rownames of TrData must match species names in Y")
        }
        if (any(is.na(TrData))) {
            stop("TrData parameter must not contain any NA values")
        }
        hM$TrData = TrData
        hM$TrFormula = TrFormula
        hM$Tr = model.matrix(TrFormula, TrData)
    }
    if (!is.null(Tr)) {
        if (!is.matrix(Tr)) {
            stop("Tr must be a matrix")
        }
        if (nrow(Tr) != hM$ns) {
            stop("the number of rows in Tr should be equal to number of columns in Y")
        }
        if (!all(rownames(Tr) == colnames(Y))) {
            stop("rownames of Tr must match species names in Y")
        }
        if (any(is.na(Tr))) {
            stop("Tr parameter must not contain any NA values")
        }
        hM$TrData = NULL
        hM$Tr = Tr
    }
    if (is.null(hM$Tr)) {
        hM$Tr = matrix(1, hM$ns, 1)
    }
    hM$nt = ncol(hM$Tr)
    if (is.null(colnames(hM$Tr))) {
        colnames(hM$Tr) = sprintf(sprintf("tr%%.%dd", ceiling(log10(hM$nt))), 
                                  1:hM$nt)
    }
    hM$trNames = colnames(hM$Tr)
    if (!TrScale) {
        hM$TrScalePar = rbind(rep(0, hM$nt), rep(1, hM$nt))
        hM$TrScaled = hM$Tr
        hM$TrInterceptInd = NULL
    }
    else {
        TrInterceptInd = which(colnames(hM$Tr) %in% c("Intercept", 
                                                      "(Intercept)"))
        if (length(TrInterceptInd) > 1) {
            stop("only one column of Tr matrix can be named Intercept or (Intercept)")
        }
        if (!all(hM$Tr[, TrInterceptInd] == 1)) {
            stop("intercept column in Tr matrix must be a column of ones")
        }
        if (length(TrInterceptInd) == 1) {
            hM$TrInterceptInd = TrInterceptInd
        }
        else hM$TrInterceptInd = NULL
        TrScalePar = rbind(rep(0, hM$nt), rep(1, hM$nt))
        TrScaled = hM$Tr
        if (TrScale) {
            scaleInd = apply(hM$Tr, 2, function(a) !all(a %in% 
                                                            c(0, 1)))
        }
        else {
            scaleInd = TrScale
        }
        scaleInd[TrInterceptInd] = FALSE
        if (length(TrInterceptInd) > 0) {
            sc = scale(hM$Tr)
            TrScalePar[, scaleInd] = rbind(attr(sc, "scaled:center"), 
                                           attr(sc, "scaled:scale"))[, scaleInd]
        }
        else {
            sc = scale(hM$Tr, center = FALSE)
            TrScalePar[2, scaleInd] = attr(sc, "scaled:scale")[scaleInd]
        }
        TrScaled[, scaleInd] = sc[, scaleInd]
        hM$TrScalePar = TrScalePar
        hM$TrScaled = TrScaled
    }
    if (!is.null(C) && !is.null(phyloTree)) {
        stop("only one of phyloTree and C arguments can be specified")
    }
    if (!is.null(phyloTree)) {
        corM = vcv.phylo(phyloTree, model = "Brownian", corr = TRUE)
        corM = corM[hM$spNames, hM$spNames]
        hM$phyloTree = phyloTree
        hM$C = corM
    }
    if (!is.null(C)) {
        if (any(dim(C) != hM$ns)) {
            stop("the size of square matrix C must be equal to number of species")
        }
        hM$C = C
    }
    if (is.null(studyDesign)) {
        hM$dfPi = NULL
        hM$Pi = matrix(NA, hM$ny, 0)
        hM$np = integer(0)
        hM$nr = 0
        hM$rLNames = character(0)
        if (!is.null(ranLevels)) {
            if (length(ranLevels) > 0) {
                stop("studyDesign is empty, but ranLevels is not")
            }
        }
    }
    else {
        if (nrow(studyDesign) != hM$ny) {
            stop("the number of rows in studyDesign must be equal to number of rows in Y")
        }
        if (!all(sapply(studyDesign, is.factor))) 
            stop("studyDesign columns must be factors")
        if (!all(ranLevelsUsed %in% names(ranLevels))) {
            stop("ranLevels must contain named elements corresponding to all levels listed in ranLevelsUsed")
        }
        if (!all(ranLevelsUsed %in% colnames(studyDesign))) {
            stop("studyDesign must contain named columns corresponding to all levels listed in ranLevelsUsed")
        }
        hM$studyDesign = studyDesign
        if (length(ranLevels) && !is.list(ranLevels)) 
            stop("'ranLevels' must be a list of 'HmscRandomLevel' objects")
        if (!all(sapply(ranLevels, inherits, what = "HmscRandomLevel"))) 
            stop("'ranLevels' must be 'HmscRandomLevel' objects")
        hM$ranLevels = ranLevels
        hM$ranLevelsUsed = ranLevelsUsed
        hM$dfPi = studyDesign[, ranLevelsUsed, drop = FALSE]
        hM$rL = ranLevels[ranLevelsUsed]
        hM$rLNames = colnames(hM$dfPi)
        hM$Pi = matrix(NA, hM$ny, ncol(hM$dfPi), dimnames = list(NULL, 
                                                                 hM$rLNames))
        for (r in seq_len(ncol(hM$dfPi))) hM$Pi[, r] = as.numeric(as.factor(hM$dfPi[, 
                                                                                    r]))
        hM$np = apply(hM$Pi, 2, function(a) return(length(unique(a))))
        hM$nr = ncol(hM$Pi)
        if (truncateNumberOfFactors) {
            for (r in seq_len(hM$nr)) {
                hM$rL[[r]]$nfMax = min(hM$rL[[r]]$nfMax, hM$ns)
                hM$rL[[r]]$nfMin = min(hM$rL[[r]]$nfMin, hM$rL[[r]]$nfMax)
            }
        }
    }
    if (is.matrix(distr)) {
        if (NROW(distr) != hM$ns) 
            stop("no. of rows in distr matrix must be equal to the no. of species")
        if (NCOL(distr) < 2) 
            stop("distr matrix should have 2 columns")
    }
    knownDistributions <- c("normal", "probit", "poisson", "lognormal poisson")
    if (length(distr) == 1) {
        switch(match.arg(distr, knownDistributions), normal = {
            distr = matrix(0, hM$ns, 2)
            distr[, 1] = 1
            distr[, 2] = 1
        }, probit = {
            distr = matrix(0, hM$ns, 2)
            distr[, 1] = 2
            distr[, 2] = 0
        }, poisson = {
            distr = matrix(0, hM$ns, 2)
            distr[, 1] = 3
            distr[, 2] = 0
        }, `lognormal poisson` = {
            distr = matrix(0, hM$ns, 2)
            distr[, 1] = 3
            distr[, 2] = 1
        })
    }
    if (length(distr) > 1 && !is.matrix(distr)) {
        if (length(distr) != hM$ns) 
            stop("length of distr should be 1 or equal to the number of species")
        distr2 = matrix(0, hM$ns, 2)
        for (i in 1:hM$ns) {
            switch(match.arg(distr[i], knownDistributions), 
                   normal = {
                       distr2[i, 1] = 1
                       distr2[i, 2] = 1
                   }, probit = {
                       distr2[i, 1] = 2
                       distr2[i, 2] = 0
                   }, poisson = {
                       distr2[i, 1] = 3
                       distr2[i, 2] = 0
                   }, `lognormal poisson` = {
                       distr2[i, 1] = 3
                       distr2[i, 2] = 1
                   })
        }
        distr = distr2
    }
    if (NCOL(distr) > 2) {
        warning("keeping only two first columns of 'distr' matrix")
        distr <- distr[, 1:2, drop = FALSE]
    }
    colnames(distr) = c("family", "variance")
    if (any(distr[, 1] == 0)) {
        stop("some of the distributions ill defined")
    }
    hM$distr = distr
    if (!YScale) {
        hM$YScalePar = rbind(rep(0, hM$ns), rep(1, hM$ns))
        hM$YScaled = hM$Y
    }
    else {
        scaleInd = which(hM$distr[, 1] == 1)
        YScalePar = rbind(rep(0, hM$ns), rep(1, hM$ns))
        YScaled = hM$Y
        if (length(scaleInd) > 0) {
            sc = scale(hM$Y)
            YScalePar[, scaleInd] = rbind(attr(sc, "scaled:center"), 
                                          attr(sc, "scaled:scale"))[, scaleInd]
            YScaled[, scaleInd] = sc[, scaleInd]
        }
        hM$YScalePar = YScalePar
        hM$YScaled = YScaled
    }
    if (!is.null(Loff)) {
        if (!is.matrix(Loff)) 
            stop("Loff argument must be NULL or a numeric matrix")
        if (nrow(Loff) != hM$ny) 
            stop("number of rows in Loff argument must be equal to ny")
        if (ncol(Loff) != hM$ns) 
            stop("number of columns in Loff argument must be equal to ns")
        hM$Loff = Loff
    }
    else {
        hM$Loff = NULL
    }
    hM = setPriors(hM, setDefault = TRUE)
    hM$call <- match.call()
    hM$HmscVersion <- packageVersion("Hmsc")
    hM
}
