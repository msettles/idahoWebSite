require(Biobase)
require(affy)
require(simpleaffy)
require(xtable)
require(affyPLM)
require(RColorBrewer)
require(genefilter)
require(lattice)
### Fix for long tables, or in general for many arrays
"affyQAReport2" <- 
function (affyB, output = "pdf", outdir = file.path(getwd(), 
    "affyQA"), overwrite = FALSE, repName) 
{
## debug
# affyB <- RawData
# output = "pdf"
# outdir = outdir
# overwrite = overwrite
# repName="affyQAReport")
## end debug
    if (missing(repName)) 
        repName = deparse(substitute(affyB))
    if (!inherits(affyB, "AffyBatch")) 
        stop("QA reports can only be generated for AffyBatch instances")
    sN = sampleNames(affyB)
    lcS = lcSuffix(sN)
    if (nchar(lcS) > 0) 
        sN = gsub(paste(lcS, "$", sep = ""), "", sN)
    sN = gsub("_", ".", sN)
    sampleNames(affyB) = sN
    if (!file.exists(outdir)) 
        if (!dir.create(outdir)) 
            stop("could not create output directory")
    outdir = file.path(outdir, repName)
    if (file.exists(outdir)) 
        if (overwrite) 
            unlink(outdir, recursive = TRUE)
        else stop("report already exists")
    if (!dir.create(outdir)) 
        stop("could not create report directory")
    cd = setwd(outdir)
    on.exit(setwd(cd))
    outfiles = list(sA = "simpleaffy", RNAdeg = "rnadeg", hist = "affyhist", 
        bxp = "affybxp", NUSE = "NUSE", RLE = "RLE")
    outf = lapply(outfiles, function(x) paste(x, output, sep = "."))
    numArrays = nrow(pData(affyB))
    qcStats = qc(affyB)
    dfout = data.frame(AvBg = avbg(qcStats), ScaleF = sfs(qcStats), 
        PerCPres = percent.present(qcStats))
    tab1 = xtable(dfout, label = "table1", caption = "Average background, scale factor and percent present calls.")
    tcon = textConnection("TAB1", "w", local = TRUE)
    print(tab1, file = tcon)
    close(tcon)
    TAB1 = paste(TAB1, collapse = "\n")
    CompStats = sapply(dfout, function(x) max(x)/min(x))
    CS = paste("$", round(CompStats, 3), "$")
    txtstrs = c("Since this ratio is less than 3 there is unlikely to be a problem.  ", 
        "Since this ratio is larger than 3 there is a potential problem.  ")
    CSoutstrs = ifelse(CompStats <= 3, txtstrs[1], txtstrs[2])
    names(CSoutstrs) = NULL
    qcratios = ratios(qcStats)
    rnrat = colnames(qcratios)
    rn1 = gsub("^AFFX-", "", rnrat)
    rn2 = strsplit(rn1, "\\.")
    grumpf = if (any(listLen(rn2) != 2)) {
        rn1
    } else sapply(rn2, paste, collapse = "\n")
    colnames(qcratios) = letters[seq_len(ncol(qcratios))]
    tab2 = xtable(qcratios, label = "table2", caption = paste("3'/5' ratios. ", 
        paste(colnames(qcratios), ") ", grumpf, sep = "", collapse = " "), 
        ".", sep = ""))
    tcon = textConnection("TAB2", "w", local = TRUE)
    print(tab2, file = tcon)
    close(tcon)
    TAB2 = paste(TAB2, collapse = "\n")
    biobOut = qcStats@spikes
    cn = colnames(biobOut)
    cn = gsub("AFFX-", "", cn, fixed = TRUE)
    cn = gsub("-3_at", "", cn, fixed = TRUE)
    colnames(biobOut) = cn
    bbc = qcStats@bioBCalls
    names(bbc) = gsub(".present", "", names(bbc))
    bbO = cbind(BioBCall = bbc, round(biobOut, digits = 3))
    tab3 = xtable(bbO, label = "table3", caption = "BioB and friends")
    tcon = textConnection("TAB3", "w", local = TRUE)
    print(tab3, file = tcon)
    close(tcon)
    TAB3 = paste(TAB3, collapse = "\n")
    acol = sample(brewer.pal(8, "Dark2"), numArrays, replace = (8 < 
        numArrays))
    argb = sapply(acol, substring, first = c(2, 4, 6), last = c(3, 
        5, 7))
    argb = apply(argb, 2, function(h) as.integer(paste("0x", 
        h, sep = ""))/255)
    definecolor = paste(paste("\\definecolor{farbe", seq_along(acol), 
        "}{rgb}{", argb[1, ], ",", argb[2, ], ",", argb[3, ], 
        "}", sep = ""), collapse = "\n")
    arrayNamesInColors = paste(paste("\\textcolor{farbe", seq_along(acol), 
        "}{", sN, "}", sep = ""), collapse = ", ")
    pdf(file = outf$sA)
    plot(qcStats)
    dev.off()
    pdf(file = outf$hist)
    hist(affyB, lty = 1, col = acol)
    dev.off()
    pdf(file = outf$bxp)
    boxplot(affyB, col = acol, las = 3)
    dev.off()
    rnaDeg = AffyRNAdeg(affyB)
    pdf(file = outf$RNAdeg)
    plotAffyRNAdeg(rnaDeg, cols = acol, lwd = 2)
    dev.off()
    pp1 = preprocess(affyB)
    epp = log2(exprs(pp1))
    colnames(epp) = sN
    medArray = rowMedians(epp)
    M = epp - medArray
    A = (epp + medArray)/2
    app = 4 + 2 * (sum(numArrays > c(4, 6)))
    nfig = ceiling(numArrays/8)
    plotNames = paste("MA", 1:nfig, sep = "")
    fNames = paste(plotNames, "pdf", sep = ".")
    xlim = quantile(A, probs = 1e-04 * c(1, -1) + c(0, 1))
    ylim = quantile(M, probs = 1e-04 * c(1, -1) + c(0, 1))
    dummy.df <- data.frame(sN = factor(sN, levels = sN), x = seq_along(sN), 
        y = seq_along(sN))
    trobj <- xyplot(y ~ x | sN, dummy.df, xlim = xlim, ylim = ylim, 
        xlab = "A", ylab = "M", panel = function(x, y, ...) {
            x <- A[, x]
            y <- M[, y]
            panel.smoothScatter(x, y, ...)
        }, layout = c(app/2, 2, 1))
    id.firstpage <- seq_len(app)
    for (i in seq_len(nfig)) {
        pdf(file = fNames[i])
        id.thispage <- (i - 1) * app + id.firstpage
        id.thispage <- id.thispage[id.thispage <= numArrays]
        print(update(trobj, index.cond = list(id.thispage)))
        dev.off()
    }
    MALatex = paste("\\begin{figure}[tp]", "\\begin{center}", 
        paste("\\includegraphics{", plotNames, "}", sep = ""), 
        "\\caption{\\label{fig:ma}MA plots.", "A \\textit{reference array} array is calculated from", 
        "the median across arrays, and for each array $M$ and $A$ values are", 
        "calculated for the comparison to that reference.}", 
        "\\end{center}\\end{figure}", sep = " \n", collapse = "\n\n")
    outM = dist2(epp)
    pdf(file = "MADimage.pdf", height = 6, width = (6 - 0.7) * 
        1.25 + 0.7)
    par(mai = c(0.7, 0.7, 0.01, 0.01))
    layout(cbind(1, 2), widths = c(4, 1))
    imcol = colorRampPalette(brewer.pal(9, "RdPu"))(256)
    rg = range(outM, na.rm = TRUE)
    image(1:numArrays, 1:numArrays, outM, xlab = "", ylab = "", 
        zlim = rg, main = "", col = imcol, axes = FALSE)
    axis(1, at = 1:numArrays, labels = sN, las = 3)
    axis(2, at = 1:numArrays, labels = sN, las = 2)
    image(1, seq(rg[1], rg[2], length = length(imcol)), rbind(seq_along(imcol)), 
        xaxt = "n", ylab = "", col = imcol)
    text(1, 0, "MAD", xpd = NA)
    dev.off()
    dataPLM = fitPLM(pp1, background = FALSE, normalize = FALSE)
    pdf(file = outf$RLE)
    Mbox(dataPLM, ylim = c(-1, 1), names = sN, col = "lightblue", 
        whisklty = 0, staplelty = 0, main = "RLE", las = 3)
    dev.off()
    pdf(file = outf$NUSE)
    boxplot(dataPLM, ylim = c(0.95, 1.5), names = sN, outline = FALSE, 
        col = "lightblue", main = "NUSE", las = 3)
    dev.off()
    texTemp = system.file("Templates/affyQAtemplate.tex", package = "affyQCReport")
    pkVers = packageDescription("affyQCReport")$Version
    sessInfo = paste(toLatex(sessionInfo()), collapse = "\n")
    symVals = c(repName = repName, outfiles, TABLE1 = TAB1, TABLE2 = TAB2, 
        TABLE3 = TAB3, MAPLOTS = MALatex, MADimage = "MADimage", 
        affyQCVersNO = pkVers, sessionInfo = sessInfo, definecolor = definecolor, 
        arrayNamesInColors = arrayNamesInColors, numArrays = as.character(numArrays), 
        chipName = affyB@cdfName, BGRATIO = CS[1], SFRATIO = CS[2], 
        PPRATIO = CS[3], BGRATIOTEXT = CSoutstrs[1], SFRATIOTEXT = CSoutstrs[2], 
        PPRATIOTEXT = CSoutstrs[3])
        
##***********##
##    outFile = file.path(outdir, paste(repName, ".tex", sep = ""))
    outFile = paste(repName, ".tex", sep = "")
##***********##

    copySubstitute(texTemp, outFile, symbolValues = symVals)
    if (.Platform$OS.type == "windows") 
        outFile = shortPathName(outFile)
    syscall = paste("pdflatex", outFile)
    system(syscall)
    system(syscall)
    return(list(qcStats = qcStats, affyPLM = dataPLM, MADS = outM, 
        loc = outdir, name = repName))
}

