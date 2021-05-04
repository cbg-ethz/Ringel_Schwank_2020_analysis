
Sgenes <- c("ARID1A", "ARID2", "SMARCA4", "SMARCB1", "CNIH4", "KEAP1", "FAM122A", "NBAS")

core <- c("TGFBR1", "TGFBR2", "SMAD4", "SMAD3")

all <- c(Sgenes, core, "APC")

all <- sort(Sgenes)

if (set %in% c("COADREAD", "COAD")) {

    type <- "TCGA-COAD"

    path <- "mutclust/"

    load(paste0(path, type, "_final.rda"))

    colnames(DF)[1:ncol(D)] <- colnames(D)
    rownames(D) <- rownames(DF)

    if (!file.exists(paste0(path, type, "_biomart.rda"))) {
        library(biomaRt)

        ensembl=useMart("ensembl")
        ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
        filters = listFilters(ensembl)
        attributes = listAttributes(ensembl)
        geneinfo <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'),
                          filters = 'ensembl_gene_id',
                          values = rownames(D),
                          mart = ensembl)
        save(geneinfo, file = paste0(path, type, "_biomart.rda"))
    } else {
        load(paste0(path, type, "_biomart.rda"))
    }

    D <- D[order(rownames(D)), ]
    geneinfo <- geneinfo[which(duplicated(geneinfo[, 2]) == FALSE), ]
    geneinfo <- geneinfo[order(geneinfo[, 2]), ]

    print(all(rownames(D) == geneinfo[, 2]))

    DF <- DF[order(rownames(DF)), ]

    rownames(DF) <- geneinfo[, 1]

    DM <- DF[, which(colnames(D) %in% colnames(M))]

    M <- M[, colnames(DM)]

    Mb <- M
    DMb <- DM
    Db <- D

}

if (set %in% c("COADREAD", "READ")) {

    type <- "TCGA-READ"

    load(paste0(path, type, "_final.rda"))

    colnames(DF)[1:ncol(D)] <- colnames(D)
    rownames(D) <- rownames(DF)

    if (!file.exists(paste0(path, type, "_biomart.rda"))) {
        library(biomaRt)

        ensembl=useMart("ensembl")
        ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
        filters = listFilters(ensembl)
        attributes = listAttributes(ensembl)
        geneinfo <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'),
                          filters = 'ensembl_gene_id',
                          values = rownames(D),
                          mart = ensembl)
        save(geneinfo, file = paste0(path, type, "_biomart.rda"))
    } else {
        load(paste0(path, type, "_biomart.rda"))
    }

    D <- D[order(rownames(D)), ]
    geneinfo <- geneinfo[which(duplicated(geneinfo[, 2]) == FALSE), ]
    geneinfo <- geneinfo[order(geneinfo[, 2]), ]

    print(all(rownames(D) == geneinfo[, 2]))

    DF <- DF[order(rownames(DF)), ]

    rownames(DF) <- geneinfo[, 1]

    DM <- DF[, which(colnames(D) %in% colnames(M))]

    M <- M[, colnames(DM)]

    ## sum(!(c(core, Sgenes) %in% rownames(M)))

    M <- rbind(M, TGFBR2 = 0)

}


if (set %in% c("COADREAD")) {

    ## COADREAD:

    M <- M[order(rownames(M)), ]

    Mb <- Mb[order(rownames(Mb)), ]

    Mgenes <- intersect(rownames(M), rownames(Mb))

    M <- cbind(M[which(rownames(M) %in% Mgenes), ], Mb[which(rownames(Mb) %in% Mgenes), ])

    DMgenes <- intersect(rownames(D), rownames(Db))

    DM <- cbind(DM[which(rownames(D) %in% DMgenes), ], DMb[which(rownames(Db) %in% DMgenes), ])

    D <- cbind(D[which(rownames(D) %in% DMgenes), ], Db[which(rownames(Db) %in% DMgenes), ])

}

## ## my plan:

## DMlog <- log2(DM[core, ] + 1)

## hist(DMlog)

## tgfpos <- which(apply(DMlog, 2, mean) > mean(DMlog) + 1*sd(DMlog))

## tgfneg <- which(apply(DMlog, 2, mean) < mean(DMlog) - (-1)*sd(DMlog))

## design <- matrix(0, ncol(DM), 4)

## design[, 1] <- 0.5

## colnames(design) <- c("TGF", "APCKO", "ARID1AKO", "SMARCA4KO")

## design[tgfpos, "TGF"] <- 1; design[tgfneg, "TGF"] <- 0 # expression tgf

## ## design[which(M["SMAD4", ] == 0), "TGF"] <- 1; design[, "noTGF"] <- 1 - design[, "TGF"] # mutation

## design[which(M["APC", ] > 0), "APCKO"] <- 1

## design[which(M["ARID1A", ] > 0), "ARID1AKO"] <- 1

## design[which(M["SMARCA4", ] > 0), "SMARCA4KO"] <- 1

## DSMURF2 <- list(A = DM["SMURF2", which(design[, 2] == 1 & design[, 1] == 0)],
##                  B = DM["SMURF2", which(design[, 2] == 1 & design[, 1] == 1)],
##                  C = DM["SMURF2", which(design[, 2] == 1 & design[, 1] == 0 & design[, 3] == 1)],
##                  D = DM["SMURF2", which(design[, 2] == 1 & design[, 1] == 1 & design[, 3] == 1)],
##                  E = DM["SMURF2", which(design[, 2] == 1 & design[, 1] == 0 & design[, 4] == 1)],
##                  F = DM["SMURF2", which(design[, 2] == 1 & design[, 1] == 1 & design[, 4] == 1)])

## DSKIL <- list(A = DM["SKIL", which(design[, 2] == 1 & design[, 1] == 0)],
##                  B = DM["SKIL", which(design[, 2] == 1 & design[, 1] == 1)],
##                  C = DM["SKIL", which(design[, 2] == 1 & design[, 1] == 0 & design[, 3] == 1)],
##                  D = DM["SKIL", which(design[, 2] == 1 & design[, 1] == 1 & design[, 3] == 1)],
##                  E = DM["SKIL", which(design[, 2] == 1 & design[, 1] == 0 & design[, 4] == 1)],
##                  F = DM["SKIL", which(design[, 2] == 1 & design[, 1] == 1 & design[, 4] == 1)])

## DITGA5 <- list(A = DM["ITGA5", which(design[, 2] == 1 & design[, 1] == 0)],
##                  B = DM["ITGA5", which(design[, 2] == 1 & design[, 1] == 1)],
##                  C = DM["ITGA5", which(design[, 2] == 1 & design[, 1] == 0 & design[, 3] == 1)],
##                  D = DM["ITGA5", which(design[, 2] == 1 & design[, 1] == 1 & design[, 3] == 1)],
##                  E = DM["ITGA5", which(design[, 2] == 1 & design[, 1] == 0 & design[, 4] == 1)],
##                  F = DM["ITGA5", which(design[, 2] == 1 & design[, 1] == 1 & design[, 4] == 1)])

## par(mfrow=c(1,1))
## boxplot(c(DSMURF2, DSKIL, DITGA5), col = rep(rgb(c(0,0,1), c(0,1,0), c(1,0,0), 0.5), each = 2), xaxt = "n")
## axis(1, 1:18, rep(c("-", "+"), 9))
## par(las=3)
## axis(1, seq(1.5, 17.5, 2), rep(c("APCKO", "APCKO/ARID1AKO", "APCKO/SMARCA4KO"), 3), cex.axis = 0.6, tick = 0)
## abline(v = c(6.5, 12.5), lty = 2)

## Till's suggestion:

library(edgeR)

checkSingles <- function(x, disc = 1, verbose = FALSE) {
    if (disc == 1) {
        x[which(x > 0)] <- 1
    } else if (disc == 2) {
        x[which(x < max(x))] <- 0
        x[which(x > 0)] <- 1
    }
    y <- apply(x, 2, sum)
    z <- apply(x, 1, function(u) {
        v <- sum(u == 1 & y == 1)
        return(v)
    })
    if (verbose) {
        print(z)
    }
    return(list(x = x, z = z))
}

MC <- checkSingles(M, #[c("TGFBR1", "TGFBR2", "SMAD3", "SMAD4", "ARID1A", "SMARCA4", "ARID2", "SMARCB1"), ],
                   disc = 2)

M1 <- MC$x

design <- matrix(0, ncol(DM), 2)

colnames(design) <- c("A", "B")

if (comp %in% "AB") {

    ## A vs B:

    design[which(apply(M1[c("TGFBR1", "TGFBR2", "SMAD3", "SMAD4"), ], 2, sum) == 0 & apply(M1[c("ARID1A", "SMARCA4", "ARID2", "SMARCB1"), ], 2, sum) > 0), 1] <- 1

    design[which(apply(M1[c("TGFBR1", "TGFBR2", "SMAD3", "SMAD4"), ], 2, sum) == 0 & apply(M1[c("ARID1A", "SMARCA4", "ARID2", "SMARCB1"), ], 2, sum) == 0), 2] <- 1

}

if (comp %in% "A1B") {

    ## A1 vs B:

    design[which(apply(M1[c("TGFBR1", "TGFBR2", "SMAD3", "SMAD4"), ], 2, sum) == 0 & apply(M1["ARID1A", , drop = FALSE], 2, sum) > 0), 1] <- 1

    design[which(apply(M1[c("TGFBR1", "TGFBR2", "SMAD3", "SMAD4"), ], 2, sum) == 0 & apply(M1[c("ARID1A", "SMARCA4", "ARID2", "SMARCB1"), ], 2, sum) == 0), 2] <- 1
}

if (comp %in% "A2B") {

    ## A2 vs B:

    design[which(apply(M1[c("TGFBR1", "TGFBR2", "SMAD3", "SMAD4"), ], 2, sum) == 0 & apply(M1["SMARCA4", , drop = FALSE], 2, sum) > 0), 1] <- 1

    design[which(apply(M1[c("TGFBR1", "TGFBR2", "SMAD3", "SMAD4"), ], 2, sum) == 0 & apply(M1[c("ARID1A", "SMARCA4", "ARID2", "SMARCB1"), ], 2, sum) == 0), 2] <- 1

}

if (comp %in% "full") {

    ## A1 vs A2 vs B:

    design <- matrix(0, ncol(DM), 3)

    design[which(apply(M1[c("TGFBR1", "TGFBR2", "SMAD3", "SMAD4"), ], 2, sum) == 0 & apply(M1["ARID1A", , drop = FALSE], 2, sum) > 0), 1] <- 1

    design[which(apply(M1[c("TGFBR1", "TGFBR2", "SMAD3", "SMAD4"), ], 2, sum) == 0 & apply(M1["SMARCA4", , drop = FALSE], 2, sum) > 0), 2] <- 1

    design[which(apply(M1[c("TGFBR1", "TGFBR2", "SMAD3", "SMAD4"), ], 2, sum) == 0 & apply(M1[c("ARID1A", "SMARCA4", "ARID2", "SMARCB1"), ], 2, sum) == 0), 3] <- 1

    colnames(design) <- c("A1", "A2", "B")

}

## design <- M1

epiNEM::HeatmapOP(design, col = "RdBu")

print(apply(design, 2, sum))

entrez <- rownames(D)

high <- which(apply(DM, 1, median) >= 10)

DM1 <- DM[high, ]

entrez <- entrez[high]

y <- DGEList(counts=DM1)

y <- calcNormFactors(y)

y <- estimateDisp(y, design)

fit <- glmQLFit(y,design)

if (comp %in% "full") {

    ## A1/A2 vs B:

    for (con in 1:2) {

        if (con == 1) {

            qlf <- glmQLFTest(fit,contrast=c(1,0,-1))

        } else {

            qlf <- glmQLFTest(fit,contrast=c(0,1,-1))

        }

        ## show results:

        topTags(qlf)

        tmp <- qlf$table

        topsort <- order(tmp$PValue)

        tmp <- cbind(ENTREZ = entrez[topsort], tmp[topsort, ], FDR = p.adjust(tmp$PValue[topsort], method = "BH")) # 233, 4, 5

        rownames(tmp)[grep("^\\.[0-9]", rownames(tmp))] <- paste0("no_symbol_", seq_len(length(grep("^\\.[0-9]", rownames(tmp)))))

        write.csv(tmp, file = paste0("DE_analysis_", set, "_", comp, "_A", i, ".csv"))

        hist(qlf$table$PValue)

        myVolcano <- function(x, sigfc = 1, FDR = FALSE, ...) {
            pvs <- x$table$PValue
            if (FDR) { pvs <- p.adjust(pvs, method = "BH") }
            logfcs <- qlf$table$logFC
            plot(logfcs, -log(pvs), ...)
            abline(h=-log(0.05), v = c(-sigfc, sigfc), col = "red")
        }

        pdf(paste0("DE_analysis_", set, "_", comp, "_A", i, ".pdf"), width = 10, height = 10)
        myVolcano(qlf, col = rgb(0,0,0,0.5), xlab = "log foldchanges", ylab = "-log p-values")
        dev.off()

    }

} else {

    ## A/A1/A2 vs B:

    qlf <- glmQLFTest(fit,contrast=c(1,-1))

    ## show results:

    topTags(qlf)

    tmp <- qlf$table

    topsort <- order(tmp$PValue)

    tmp <- cbind(ENTREZ = entrez[topsort], tmp[topsort, ], FDR = p.adjust(tmp$PValue[topsort], method = "BH")) # 233, 4, 5

    rownames(tmp)[grep("^\\.[0-9]", rownames(tmp))] <- paste0("no_symbol_", seq_len(length(grep("^\\.[0-9]", rownames(tmp)))))

    write.csv(tmp, file = paste0("DE_analysis_", set, "_", comp, ".csv"))

    hist(qlf$table$PValue)

    myVolcano <- function(x, sigfc = 1, FDR = FALSE, ...) {
        pvs <- x$table$PValue
        if (FDR) { pvs <- p.adjust(pvs, method = "BH") }
        logfcs <- qlf$table$logFC
        plot(logfcs, -log(pvs), ...)
        abline(h=-log(0.05), v = c(-sigfc, sigfc), col = "red")
    }

    pdf(paste0("DE_analysis_", set, "_", comp, ".pdf"), width = 10, height = 10)
    myVolcano(qlf, col = rgb(0,0,0,0.5), xlab = "log foldchanges", ylab = "-log p-values")
    dev.off()

}

stop()

## run script:

set <- "COAD"; comp <- "AB"; source("~/Documents/testing/vignettes/schwankIII.r")

set <- "COADREAD"; comp <- "AB"; source("~/Documents/testing/vignettes/schwankIII.r")

set <- "COAD"; comp <- "A1B"; source("~/Documents/testing/vignettes/schwankIII.r")

set <- "COADREAD"; comp <- "A1B"; source("~/Documents/testing/vignettes/schwankIII.r")

set <- "COAD"; comp <- "A2B"; source("~/Documents/testing/vignettes/schwankIII.r")

set <- "COADREAD"; comp <- "A2B"; source("~/Documents/testing/vignettes/schwankIII.r")

set <- "COAD"; comp <- "full"; source("~/Documents/testing/vignettes/schwankIII.r")

set <- "COADREAD"; comp <- "full"; source("~/Documents/testing/vignettes/schwankIII.r")

## overlap:

iv <- 1

plot(1, length(intersect(rownames(COAD)[1:iv], rownames(COADREAD)[1:iv])), xlim = c(0,110/iv), ylim = c(0,110))

for (i in (2:(100/iv))) {
    lines(i, length(intersect(rownames(COAD)[1:(i*iv)], rownames(COADREAD)[1:(i*iv)])), type = "p")
}
abline(0,iv)

## MSS:

if (!file.exists(paste0(path, type, "_msi_mss.rda"))) {
    query <- GDCquery(project = type,
                      data.category = "Other",
                      legacy = TRUE,
                      access = "open",
                      data.type = "Auxiliary test",
                      barcode = gsub("-01$", "", colnames((D))))
    GDCdownload(query)
    msires <- GDCprepare_clinic(query, "msi")
    save(msires, file = paste0(path, type, "_msi_mss.rda"))
} else {
    load(paste0(path, type, "_msi_mss.rda"))
}

## check mutation amount with mss/msi status:

colnames(M) <- gsub("-01$", "", colnames(M))

MSS <- M[, which(colnames(M) %in% msires[which(msires[, 3] %in% "MSS"), 1])]

MSI <- M[, which(colnames(M) %in% msires[which(msires[, 3] %in% c("MSI-L", "MSI-H")), 1])]

MSI[which(MSI < 4)] <- 0

MSS[which(MSS < 4)] <- 0

if (type %in% "TCGA-COAD") {
    breaks <- breaks2 <- 50
    ymax <- 200
} else {
    breaks <- 100
    breaks2 <- 1000
    ymax <- 25
}

pdf("temp.pdf", height = 10, width = 10)
MSIsum <- apply(MSI, 2, sum)
MSSsum <- apply(MSS, 2, sum)
hist(MSIsum, col = rgb(1,0,0,0.5), breaks = breaks, ylim = c(0, ymax), main = "distribution of mutation rate of patients\n(MSI: red, MSS: green)", xlab = "number of mutations", ylab = "number of patients")
hist(MSSsum, col = rgb(0,1,0,0.5), add = TRUE, breaks = breaks2)
dev.off()

## new analysis:

library(survival)

survdat <- clinical[which(clinical[, 1] %in% msires[, 1] & clinical[, 1] %in% gsub("-01$", "", colnames(M))), ]

survdat <- survdat[order(survdat[, 1]), ]

tmp <- t(M[, which(gsub("-01$", "", colnames(M)) %in% survdat[, 1])])

rownames(tmp) <- gsub("-01$", "", rownames(tmp))

tmp <- tmp[survdat[, 1], ]

tmp[tmp > 0] <- 1

tmp2 <- msires

rownames(tmp2) <- tmp2[, 1]

tmp2 <- as.numeric(tmp2[survdat[, 1], 3])

tmp2[which(tmp2 %in% c(2,3))] <- 2
tmp2[which(tmp2 %in% 4)] <- 3

survdat <- cbind(survdat, tmp[, all], msimss = tmp2)

sfits <- list()
coxs <- list()
pvals <- numeric(length(all))
for (j in 1:length(all)) {
    i <- all[j]
    coxs[[i]] <- coxph(Surv(days_to_death, vital_status) ~ get(i)+tumor_stage+age_at_diagnosis+msimss, survdat)
    pvals[j] <- summary(coxs[[i]])$coefficients[1, 5]
}
qvals <- p.adjust(pvals, method = "BH")

## pdf(paste0("Schwank/", type, "_mut.pdf"), width = 5, height = 5)
pdf(paste0(type, "_mut.pdf"), width = 5, height = 5)
for (j in 1:length(all)) {
    i <- all[j]
    sfits[[i]] <- survfit(Surv(days_to_death, vital_status) ~ get(i), survdat)
    plot(sfits[[i]], col = 1:2, ylab = "Survival", xlab = "Time in days", lwd = 2, mark.time = TRUE,
         main = paste0(round(qvals[j], 3), " (", round(pvals[j], 3), ")"))
    legend(max(sfits[[i]]$time), 0.9, c(paste0("wild type (",
                                               sum(survdat[, grep(i, colnames(survdat))] == 0), ")"),
                                        paste0(i, " mutated\n(",
                                               sum(survdat[, grep(i, colnames(survdat))] == 1)
                                              ,")")),
           col=1:2, lwd=2, bty='n', xjust = 1, cex = 1)
}
dev.off()

##

survtmp <- survdat[which(survdat$msimss == 3), ]

sfits <- list()
coxs <- list()
pvals <- numeric(length(all))
for (j in 1:length(all)) {
    i <- all[j]
    coxs[[i]] <- coxph(Surv(days_to_death, vital_status) ~ get(i)+tumor_stage+age_at_diagnosis+msimss, survtmp)
    pvals[j] <- summary(coxs[[i]])$coefficients[1, 5]
}
qvals <- p.adjust(pvals, method = "BH")

## pdf(paste0("Schwank/", type, "_mut.pdf"), width = 5, height = 5)
pdf(paste0(type, "_mut.pdf"), width = 5, height = 5)
for (j in 1:length(all)) {
    i <- all[j]
    sfits[[i]] <- survfit(Surv(days_to_death, vital_status) ~ get(i), survtmp)
    plot(sfits[[i]], col = 1:2, ylab = "Survival", xlab = "Time in days", lwd = 2, mark.time = TRUE,
         main = paste0(round(qvals[j], 3), " (", round(pvals[j], 3), ")"))
    legend(max(sfits[[i]]$time), 0.9, c(paste0("wild type (",
                                               sum(survtmp[, grep(i, colnames(survtmp))] == 0), ")"),
                                        paste0(i, " mutated\n(",
                                               sum(survtmp[, grep(i, colnames(survtmp))] == 1)
                                              ,")")),
           col=1:2, lwd=2, bty='n', xjust = 1, cex = 1)
}
dev.off()

## permutation test for overlap of genelists:

A2full <- read.delim("result--ARID1A--over--APC-0h.txt")
S2full <- read.delim("result--SMARCA4--over--APC-0h.txt")

type <- "COADREAD"

S1 <- read.csv(paste0("DE_analysis_", type, "_A1B_smarca4.csv"))
A1 <- read.csv(paste0("DE_analysis_", type, "_A2B_arid1a.csv"))

S2 <- read.csv("result--SMARCA4--over--APC-0h_small.csv")
A2 <- read.csv("result--ARID1A--over--APC-0h_small.csv")

print(all(A2full$Identifier == A2$Identifier))
print(all(S2full$Identifier == S2$Identifier))

S2 <- S2[which(S2full$isPresent == TRUE), ]
A2 <- A2[which(A2full$isPresent == TRUE), ]

## order for same columns

A1 <- A1[, c(2,1,6,7,3)]
S1 <- S1[, c(2,1,6,7,3)]

## only use genes which are in both sets?

# A1 <- A1[which(A1[, 1] %in% A2[, 1]), ]
# S1 <- S1[which(S1[, 1] %in% S2[, 1]), ]
# A2 <- A2[which(A2[, 1] %in% A1[, 1]), ]
# S2 <- S2[which(S2[, 1] %in% S1[, 1]), ]

## order for pvale or others

orderfor <- 3 # pval 3, fdr 4, logfc 5

if (orderfor == 3) {
S1 <- S1[order(S1[, orderfor]), ]
A1 <- A1[order(A1[, orderfor]), ]
S2 <- S2[order(S2[, orderfor]), ]
A2 <- A2[order(A2[, orderfor]), ]
} else {
S1 <- S1[order(abs(S1[, orderfor]), decreasing = TRUE), ]
A1 <- A1[order(abs(A1[, orderfor]), decreasing = TRUE), ]
S2 <- S2[order(abs(S2[, orderfor]), decreasing = TRUE), ]
A2 <- A2[order(abs(A2[, orderfor]), decreasing = TRUE), ]
}

## get up and down regulated lists separate

S1up <- S1[which(S1[, 5] > 0), ]
S2up <- S2[which(S2[, 5] > 0), ]
A1up <- A1[which(A1[, 5] > 0), ]
A2up <- A2[which(A2[, 5] > 0), ]

S1dn <- S1[which(S1[, 5] < 0), ]
S2dn <- S2[which(S2[, 5] < 0), ]
A1dn <- A1[which(A1[, 5] < 0), ]
A2dn <- A2[which(A2[, 5] < 0), ]

## compute overlap for certain cutoffs

cutoff <- (1:10)*100
cutoff <- 100*(1:10)
#cutoff <- 25*(1:8)
overlap1 <- matrix(0, length(cutoff), 4)
colnames(overlap1) <- c("SMARCA4 up", "SMARCA4 dn", "ARID1A up", "ARID1A dn")
rownames(overlap1) <- cutoff
for (n in seq_len(length(cutoff))) {
  overlap1[n, 1] <- sum(S1up[1:cutoff[n], 1] %in% S2up[1:cutoff[n], 1])
  overlap1[n, 2] <- sum(S1dn[1:cutoff[n], 1] %in% S2dn[1:cutoff[n], 1])
  overlap1[n, 3] <- sum(A1up[1:cutoff[n], 1] %in% A2up[1:cutoff[n], 1])
  overlap1[n, 4] <- sum(A1dn[1:cutoff[n], 1] %in% A2dn[1:cutoff[n], 1])
}

print(overlap1)

## permutation test:

# perms <- 10000 # too low?
# overlap <- array(0, c(length(cutoff), 4, perms),
#                  list(cutoff = paste0("cutoff", cutoff),
#                       comps = paste0("comparison", 1:4),
#                       perms = paste0("permutations", 1:perms)))
# for (m in seq_len(perms)) {
#   for (n in seq_len(length(cutoff))) {
#     ## prob only necessary to permute one list?
#     S1up <- S1up[sample(1:nrow(S1up), nrow(S1up)), ]
#     S1dn <- S1dn[sample(1:nrow(S1dn), nrow(S1dn)), ]
#     A1up <- A1up[sample(1:nrow(A1up), nrow(A1up)), ]
#     A1dn <- A1dn[sample(1:nrow(A1dn), nrow(A1dn)), ]
#     S2up <- S2up[sample(1:nrow(S2up), nrow(S2up)), ]
#     S2dn <- S2dn[sample(1:nrow(S2dn), nrow(S2dn)), ]
#     A2up <- A2up[sample(1:nrow(A2up), nrow(A2up)), ]
#     A2dn <- A2dn[sample(1:nrow(A2dn), nrow(A2dn)), ]
#     overlap[n, 1, m] <- sum(S1up[1:cutoff[n], 1] %in% S2up[1:cutoff[n], 1])
#     overlap[n, 2, m] <- sum(S1dn[1:cutoff[n], 1] %in% S2dn[1:cutoff[n], 1])
#     overlap[n, 3, m] <- sum(A1up[1:cutoff[n], 1] %in% A2up[1:cutoff[n], 1])
#     overlap[n, 4, m] <- sum(A1dn[1:cutoff[n], 1] %in% A2dn[1:cutoff[n], 1])
#   }
# }

doPerm <- function(m, cutoff) {
  overlap <- matrix(0, length(cutoff), 4)
  for (n in seq_len(length(cutoff))) {
    ## prob only necessary to permute one list?
    S1up <- S1up[sample(1:nrow(S1up), nrow(S1up)), ]
    S1dn <- S1dn[sample(1:nrow(S1dn), nrow(S1dn)), ]
    A1up <- A1up[sample(1:nrow(A1up), nrow(A1up)), ]
    A1dn <- A1dn[sample(1:nrow(A1dn), nrow(A1dn)), ]
    S2up <- S2up[sample(1:nrow(S2up), nrow(S2up)), ]
    S2dn <- S2dn[sample(1:nrow(S2dn), nrow(S2dn)), ]
    A2up <- A2up[sample(1:nrow(A2up), nrow(A2up)), ]
    A2dn <- A2dn[sample(1:nrow(A2dn), nrow(A2dn)), ]
    overlap[n, 1] <- sum(S1up[1:cutoff[n], 1] %in% S2up[1:cutoff[n], 1])
    overlap[n, 2] <- sum(S1dn[1:cutoff[n], 1] %in% S2dn[1:cutoff[n], 1])
    overlap[n, 3] <- sum(A1up[1:cutoff[n], 1] %in% A2up[1:cutoff[n], 1])
    overlap[n, 4] <- sum(A1dn[1:cutoff[n], 1] %in% A2dn[1:cutoff[n], 1])
  }
  return(overlap)
}

perms <- 10000
library(snowfall)
sfInit(parallel = TRUE, cpus = 1)
sfExport("cutoff", "doPerm", "S1up", "S2up", "S1dn", "S2dn", "A1up", "A2up", "A1dn", "A2dn")
overlap <- array(unlist(sfLapply(seq_len(perms), doPerm, cutoff)), dim = c(length(cutoff), 4, perms))
sfStop()

## save(overlap, file = paste0(type, "_perm_test_overlap_", orderfor, "_", max(cutoff), "_", perms, ".rda"))

load(paste0(type, "_perm_test_overlap_", orderfor, "_", max(cutoff), "_", perms, ".rda"))

## compute statistics

permpval <- overlap1*0 + 1
fishpval <- overlap1*0 + 1
overmean <- overlap1*0
over5 <- overX <- overlap1*0
for (i in 1:length(cutoff)) {
  for (j in 1:4) {
    permpval[i, j] <- sum(overlap[i, j, ] >= overlap1[i, j])/perms
    overmean[i, j] <- median(overlap[i, j, ])
    over5[i, j] <- quantile(overlap[i, j, ], 0.95)
    overX[i, j] <- quantile(overlap[i, j, ], 0.5)
    fishpval[i, j] <- fisher.test(matrix(c(overlap1[i, j],
                                             cutoff[i]-overlap1[i, j],
                                             cutoff[i]-overlap1[i, j],
                                             nrow(S1)-cutoff[i]*2), 2), alternative = "greater")$p.value
  }
}

print(overmean)
print(over5)
print(fishpval)
print(permpval)

## plot results

plot(overlap1[, 1], type = "b", col = 2, pch = 1, xaxt = "n", ylab = "# overlap", xlab = "top n genes")
axis(1, 1:length(cutoff), cutoff)
lines(overlap1[, 2], type = "b", col = 2, pch = 2)
lines(overlap1[, 3], type = "b", col = 3, pch = 3)
lines(overlap1[, 4], type = "b", col = 3, pch = 4)

legendy <- max(cbind(overlap1, apply(over5, 1, max)))

legend(1, legendy, c("SMARCA4-up", "SMARCA4-down", "ARID1A-up", "ARID1A-down", "~5% significance level"),
       col = c(2,2,3,3,4),
       pch = c(1,2,3,4,5), cex = 1.75)

lines(apply(over5, 1, max), type = "b", col = 4, pch = 5)

# lines(apply(overX, 1, max), type = "b", col = 5, pch = 6)
