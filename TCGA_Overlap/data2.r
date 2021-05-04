
## bsub -R "rusage[mem=1000]" -q normal.4h -n 1 -e schwank_err.txt -o schwank_out.txt < testing/vignettes/schwankII.r

## colorectal cancer

## CRISPR screen from organoids identifies several interesting genes.

## TCGA analysis: (to use the gdc query you have to set your internal machine clock to before 15th of june 2018!)

## possible focus on subtypes: stable (low mut. freq., APC, annotated?); (early stage tumors?)

Sgenes <- c("ARID1A", "ARID2", "SMARCA4", "SMARCB1", "CNIH4", "KEAP1", "FAM122A", "NBAS")

core <- c("TGFBR1", "TGFBR2", "SMAD4", "SMAD3")

all <- c(Sgenes, core, "APC")

all <- sort(Sgenes)

if (apc) { all <- c(all, "APC") }

chrom <- c(1,12,19,22,1,19,9,2,9,3,18,15)

strand <- c(1,1,1,1,1,0,1,0,1,1,1,1)

start <- c(26693236,45729665,10960922,23786963,224356811,10486120,68780034,99104038,30606493,51028394,67063763)

end <- c(26782110,4590800,11065395,23834518,10503741,68784608,99154192,30694142,51085045,67195195)

library(TCGAbiolinks)

## ## variables:
## takeall <- FALSE ## take all mutations from all methods
## checksize <- TRUE ## for best kmeans k (FALSE is k=3)
## giant <- TRUE ## make big plots with all curves on one page

## takeall <- TRUE; checksize <- TRUE; giant <- TRUE; source("testing/vignettes/schwankII.r")
## takeall <- TRUE; checksize <- FALSE; giant <- TRUE; source("testing/vignettes/schwankII.r")
## takeall <- FALSE; checksize <- TRUE; giant <- TRUE; source("testing/vignettes/schwankII.r")
## takeall <- FALSE; checksize <- FALSE; giant <- TRUE; source("testing/vignettes/schwankII.r")
## takeall <- TRUE; checksize <- TRUE; giant <- FALSE; source("testing/vignettes/schwankII.r")
## takeall <- TRUE; checksize <- FALSE; giant <- FALSE; source("testing/vignettes/schwankII.r")
## takeall <- FALSE; checksize <- TRUE; giant <- FALSE; source("testing/vignettes/schwankII.r")
## takeall <- FALSE; checksize <- FALSE; giant <- FALSE; source("testing/vignettes/schwankII.r")

## use that
## apc <- FALSE; takeall <- FALSE; checksize <- TRUE; giant <- TRUE; source("testing/vignettes/schwankII.r")
## apc <- FALSE; takeall <- FALSE; checksize <- FALSE; giant <- TRUE; source("testing/vignettes/schwankII.r")
## apc <- TRUE; takeall <- FALSE; checksize <- TRUE; giant <- TRUE; source("testing/vignettes/schwankII.r")
## apc <- TRUE; takeall <- FALSE; checksize <- FALSE; giant <- TRUE; source("testing/vignettes/schwankII.r")

for (type in c("COAD", "READ", "both")) {

    ## mutations:
    ## type <- "COAD" # type <- "READ"
    ## mutation1 <- GDCquery_Maf(tumor = type, save.csv = TRUE, pipeline = "varscan2") # mutation file
    ## mutation2 <- GDCquery_Maf(tumor = type, save.csv = TRUE, pipeline = "muse") # mutation file
    ## mutation3 <- GDCquery_Maf(tumor = type, save.csv = TRUE, pipeline = "somaticsniper") # mutation file
    ## mutation4 <- GDCquery_Maf(tumor = type, save.csv = TRUE, pipeline = "mutect2") # mutation file

    if (file.exists(paste0(type, "_mut.rda"))) {
        load(paste0(type, "_mut.rda"))
    } else {
        library(data.table)
        if (type %in% "COAD") {
            mutation1 <- fread("GDCdata/TCGA.COAD.varscan.8177ce4f-02d8-4d75-a0d6-1c5450ee08b0.DR-10.0.somatic.maf.csv")
            mutation2 <- fread("GDCdata/TCGA.COAD.muse.70cb1255-ec99-4c08-b482-415f8375be3f.DR-10.0.somatic.maf.csv")
            mutation3 <- fread("GDCdata/TCGA.COAD.somaticsniper.70835251-ddd5-4c0d-968e-1791bf6379f6.DR-10.0.somatic.maf.csv")
            mutation4 <- fread("GDCdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv")
        } else {
            mutation1 <- fread("GDCdata/TCGA.READ.varscan.b2689e8f-3b64-4214-8a87-dc7e7cf6fe5e.DR-10.0.somatic.maf.csv")
            mutation2 <- fread("GDCdata/TCGA.READ.muse.ec8ec3ad-f08d-46eb-9571-42806e304b37.DR-10.0.somatic.maf.csv")
            mutation3 <- fread("GDCdata/TCGA.READ.somaticsniper.e48ffb82-9208-4be3-8a47-0a1168a07054.DR-10.0.somatic.maf.csv")
            mutation4 <- fread("GDCdata/TCGA.READ.mutect.faa5f62a-2731-4867-a264-0e85b7074e87.DR-10.0.somatic.maf.csv")
        }
        library(snowfall)
        countSamples <- function(i, genes, mut.mat, coln) {
            i <- which(rownames(mut.mat) %in% genes[i])
            tmp <- mutation[which(mutation$Hugo_Symbol %in% rownames(mut.mat)[i]), coln]
            tmp2 <- mut.mat[i, ]
            tmp2[which(colnames(mut.mat) %in% tmp)] <- 1
            return(tmp2)
        }

        for (meth in 1:4) {
            mutation <- get(paste0("mutation", meth))
            allsub <- sort(unique(mutation$Hugo_Symbol)) # all[which(all %in% mutation$Hugo_Symbol)]
            mut.mat <- matrix(0, length(unique(mutation$Hugo_Symbol)), length(unique(mutation$Tumor_Sample_Barcode)))
            colnames(mut.mat) <- sort(unique(mutation$Tumor_Sample_Barcode))
            rownames(mut.mat) <- sort(unique(mutation$Hugo_Symbol))
            coln <- which(colnames(mutation) %in% "Tumor_Sample_Barcode")
            ## mutation <- mutation[which(mutation$Hugo_Symbol %in% all), ]
            sfInit(parallel = T, cpus = 4)
            sfExport("mutation", "mut.mat", "coln")
            M <- sfLapply(as.list(1:length(unique(mutation$Hugo_Symbol))), countSamples, allsub, mut.mat, coln)
            sfStop()
            tmp <- do.call("rbind", M)
            rownames(tmp) <- allsub
            colnames(tmp) <- colnames(mut.mat)
            assign(paste0("M", meth), tmp)
        }
        save(M1, M2, M3, M4, file = paste0(type, "_mut.rda"))
    }
    for (meth in 1:4) {
        tmp <- get(paste0("M", meth))
        if (meth == 1) {
            mut.mat <- tmp
        } else {
            newrow <- length(unique(c(rownames(mut.mat), rownames(tmp))))
            newcol <- length(unique(c(colnames(mut.mat), colnames(tmp))))
            new <- matrix(0, newrow, newcol)
            rownames(new) <- sort(unique(c(rownames(mut.mat), rownames(tmp))))
            colnames(new) <- sort(unique(c(colnames(mut.mat), colnames(tmp))))
            new[which(rownames(new) %in% rownames(mut.mat)), which(colnames(new) %in% colnames(mut.mat))] <- mut.mat
            new[which(rownames(new) %in% rownames(tmp)), which(colnames(new) %in% colnames(tmp))] <- new[which(rownames(new) %in% rownames(tmp)), which(colnames(new) %in% colnames(tmp))] + tmp
            mut.mat <- new
        }
    }
    if (takeall) {
        mut.mat[which(mut.mat > 0)] <- 1
    } else {
        mut.mat[which(mut.mat < 4)] <- 0
        mut.mat[which(mut.mat > 0)] <- 1
    }
    colnames(mut.mat) <- unlist(lapply(colnames(mut.mat), function(x) {
        y <- unlist(strsplit(x, "-"))
        y <- y[1:3]
        y <- paste(y, collapse = "-")
        return(y)
        }))
    mut.matr <- mut.mat[which(rownames(mut.mat) %in% all), ]
    assign(paste("mut.", type, sep = ""), mut.matr)
    colnames(mut.matr) <- unlist(lapply(colnames(mut.matr), function(x) {
        y <- unlist(strsplit(x, "-"))
        y <- y[1:3]
        y <- paste(y, collapse = "-")
        return(y)
        }))

    ## expression:

    library(SummarizedExperiment)
    library(biomaRt)
    if (file.exists(paste0(type, "_data.rda"))) {
        load(paste0(type, "_data.rda"))
    } else {
        assign(paste("query.", type, ".exprs", sep = ""), GDCquery(project = paste("TCGA-", type, sep = ""), sample.type = "Primary solid Tumor",
                                                                   data.category = "Transcriptome Profiling"
                                                                 , data.type = "Gene Expression Quantification"
                                                                 , workflow.type = "HTSeq - Counts"
                                                                   )
               )
        GDCdownload(get(paste("query.", type, ".exprs", sep = "")))
        assign(paste("exprs.", type, sep = ""), GDCprepare(get(paste("query.", type, ".exprs", sep = ""))))
        data <- assay(get(paste("exprs.", type, sep = "")))
        save(data, file = paste0(type, "_data.rda"))
    }
    if (file.exists(paste0(type, "_biomart.rda"))) {
        load(paste0(type, "_biomart.rda"))
    } else {
        mart = useMart("ensembl")
        mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
        tmp <- listAttributes(mart)
        gene.symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position", "chromosome_name"), mart = mart, values = rownames(data), filters = "ensembl_gene_id")
        save(gene.symbols, file = paste0(type, "_biomart.rda"))
    }
    gs <- gene.symbols[-which(duplicated(gene.symbols[, 1]) == TRUE), ]
    rownames(data) <- gs[, 2]
    ginfor <- gene.symbols[which(gene.symbols[, 2] %in% all), ]
    library(edgeR)
    nf <- calcNormFactors(data)
    datan <- t(t(data)*nf)
    datar <- datan[which(rownames(data) %in% all), ]
    colnames(datar) <- unlist(lapply(colnames(datar), function(x) {
        y <- unlist(strsplit(x, "-"))
        y <- y[1:3]
        y <- paste(y, collapse = "-")
        return(y)
    }))
    datar <- datar[, order(colnames(datar))]
    dups <- which(duplicated(colnames(datar)) == TRUE)
    for (i in unique(colnames(datar)[dups])) {
        tmp <- datar[, which(colnames(datar) %in% i)]
        tmp2 <- apply(tmp, 1, function(x) {
            if (any(is.na(x) == TRUE)) {
                if (all(is.na(x) == TRUE)) {
                    y <- NA
                } else {
                    y <- median(x[-which(is.na(x) == TRUE)])
                }
            } else {
                y <- median(x)
            }
            return(y)
        })
        datar[, which(colnames(datar) %in% i)[1]] <- tmp2
        datar <- datar[, -which(colnames(datar) %in% i)[-1]]
    }
    datarsave <- datar
    datar <- log(datar+1)
    library(cluster)
    for (i in all) {
        ## k-means:
        d <- dist(datar[i, ])
        bestk <- 3
        scs <- numeric(9)
        if (checksize) {
            bestsc <- 0
            for (k in 2:10) {
                tmp <- kmeans(d, k, nstart = 15)
                sc <- silhouette(tmp$cluster, d)
                sc <- sum(sc[, 3])
                if (sc > bestsc) {
                    bestsc <- sc
                    bestk <- k
                }
                scs[k-1] <- sc
            }
        }
        print(bestk)
        ##print(scs)
        #kcenters <- quantile(datar[i, ], probs = seq(0, 1, length.out = bestk))
        tmp <- kmeans(datar[i, ], bestk, nstart = 15) # centers = kcenters)
        low <- which.min(tmp$centers)
        low <- which(tmp$cluster == low)
        high <- which.max(tmp$centers)
        high <- which(tmp$cluster == high)
        tmp <- datar[i, ]*0
        tmp[low] <- -1
        tmp[high] <- 1
        datar[i, ] <- tmp
    }

    if (type %in% "READ") {
        readsamples <- colnames(datar)
    }
    if (type %in% "COAD") {
        coadsamples <- colnames(datar)
    }
    if (type %in% "both") {
        print("check if clusters are just read and coad")
        print("read")
        print(table(datar[, which(colnames(datar) %in% readsamples)]))
        print(apply(datar, 1, function(x) {
            y <- table(x[which(colnames(datar) %in% readsamples)])
            return(y)
            }))
        print("coad")
        print(table(datar[, which(colnames(datar) %in% coadsamples)]))
        print(apply(datar, 1, function(x) {
            y <- table(x[which(colnames(datar) %in% coadsamples)])
            return(y)
            }))
    }

    ## try zscores:

    dataz <- t(apply(log(datarsave+0.1), 1, function(x) {
        x <- scale(x)
        low <- which(x < -2)
        high <- which(x > 2)
        zeros <- which(abs(x) <= 2)
        x[low] <- -1
        x[zeros] <- 0
        x[high] <- 1
        return(x)
    }))
    colnames(dataz) <- colnames(datar)
    rownames(dataz) <- paste0(rownames(datar), "_expz")

    ## DawnRank:
    ## get string network:
    if (!file.exists(paste0(type, "_dawnrank.rda")) & !(type %in% "both")) {
        library(STRINGdb)
        data.normal <- GDCquery(project = paste("TCGA-", type, sep = ""), sample.type = "Solid Tissue Normal",
                                data.category = "Transcriptome Profiling"
                              , data.type = "Gene Expression Quantification"
                              , workflow.type = "HTSeq - Counts"
                                )
        GDCdownload(data.normal)
        data.normal <- GDCprepare(data.normal)
        data.normal <- assay(data.normal)
        rownames(data.normal) <- gs[, 2]
        string_db <- STRINGdb$new( version="10", species=9606, score_threshold=0.9, input_directory="" )
        graph <- string_db$get_graph()
        library(igraph)
        adj <- as.matrix(as_adjacency_matrix(graph))
        adj <- adj[order(rownames(adj)), order(colnames(adj))]
        gsII <- string_db$get_aliases()
        gsII <- gsII[order(gsII[, 1]), ]
        gsIII <- gsII[-which(duplicated(gsII[, 1]) == TRUE), ]
        adj2 <- adj
        colnames(adj2) <- rownames(adj2) <- gsIII[which(rownames(adj2) %in% gsIII[, 1]), 2]
        save(data, data.normal, adj2, mut.mat, file = paste0(type, "_dawnrank.rda"))
    }

    ## copy number (use gistic):
    if (file.exists(paste0(type, "_cnv.rda"))) {
        load(paste0(type, "_cnv.rda"))
    } else {
        cnv <- getGistic(type, type = "thresholded")
        save(cnv, file = paste0(type, "_cnv.rda"))
    }
    cnvr <- cnv[which(cnv[, 1] %in% all), ]
    rownames(cnvr) <- cnvr[, 1]
    cnvr <- cnvr[, -(1:3)]
    cnvr <- apply(cnvr, c(1,2), as.numeric)
    colnames(cnvr) <- unlist(lapply(colnames(cnvr), function(x) {
        y <- unlist(strsplit(x, "-"))
        y <- y[1:3]
        y <- paste(y, collapse = "-")
        return(y)
    }))

    ## clinical:
    if (file.exists(paste0(type, "_clin.rda"))) {
        load(paste0(type, "_clin.rda"))
    } else {
        clinical <- GDCquery_clinic(project = paste0("TCGA-", type), type = "clinical")
        save(clinical, file = paste0(type, "_clin.rda"))
    }

    ## add mut and exprs to clinical and process:
    clinical <- clinical[order(clinical$bcr_patient_barcode), ]
    mutf <- mut.matr[order(rownames(mut.matr)), order(colnames(mut.matr))]
    dataf <- datar[order(rownames(datar)), order(colnames(datar))]
    datazf <- dataz[order(rownames(dataz)), order(colnames(dataz))]
    cnvf <- cnvr[order(rownames(cnvr)), order(colnames(cnvr))]
    mutf <- mutf[, which(colnames(mutf) %in% clinical$bcr_patient_barcode)]
    dataf <- dataf[, which(colnames(dataf) %in% clinical$bcr_patient_barcode)]
    datazf <- datazf[, which(colnames(datazf) %in% clinical$bcr_patient_barcode)]
    cnvf <- cnvf[, which(colnames(cnvf) %in% clinical$bcr_patient_barcode)]
    rownames(dataf) <- paste0(rownames(dataf), "_exp")
    rownames(cnvf) <- paste0(rownames(cnvf), "_cnv")
    rownames(mutf) <- paste0(rownames(mutf), "_mut")

    ## prep further:
    clinf <- clinical
    clinf$days_to_last_follow_up[which(is.na(clinf$days_to_death) == FALSE)] <- clinf$days_to_death[which(is.na(clinf$days_to_death) == FALSE)]
    clinf$days_to_last_follow_up <- as.numeric(clinf$days_to_last_follow_up)
    clinf$vital_status[which(clinf$vital_status %in% "dead")] <- 1
    clinf$vital_status[which(clinf$vital_status %in% "alive")] <- 0
    clinf$vital_status <- as.numeric(clinf$vital_status)
    stages <- sort(unique(clinf$tumor_stage))
    for (i in 1:length(unique(clinf$tumor_stage))) {
        clinf$tumor_stage[which(clinf$tumor_stage %in% stages[i])] <- i
    }
    clinf$tumor_stage <- as.numeric(clinf$tumor_stage)

    ## start survival analysis:
    library(survival)

    ## do survival for mutations:
    print(all(clinf[which(clinf$bcr_patient_barcode %in% colnames(mutf)), 1] == colnames(mutf)))
    clinmut <- cbind(clinf[which(clinf$bcr_patient_barcode %in% colnames(mutf)), ], t(mutf))
    if (apc) {
        clinmutsave <- clinmut
        clinmut <- clinmut[which(clinmut[, 1] %in% clinmutsave[which(clinmutsave[, "APC_mut"] == 1), 1]), ]
    }
    sfits <- list()
    coxs <- list()
    pvals <- numeric(length(all))
    for (j in 1:length(all)) {
        i <- all[j]
        coxs[[i]] <- coxph(Surv(days_to_last_follow_up, vital_status) ~ get(paste0(i, "_mut"))+tumor_stage+age_at_diagnosis, clinmut)
        pvals[j] <- summary(coxs[[i]])$coefficients[1, 5]
    }
    qvals <- p.adjust(pvals, method = "BH")
    if (giant) {
        pdf(paste0("Schwank/", type, "_mut.pdf"), width = 14, height = 7)
        par(mfrow=c(2,4))
    } else {
        pdf(paste0("Schwank/", type, "_mut.pdf"), width = 5, height = 5)
    }
    for (j in 1:length(all)) {
        i <- all[j]
        sfits[[i]] <- survfit(Surv(days_to_last_follow_up, vital_status) ~ get(paste0(i, "_mut")), clinmut)
        plot(sfits[[i]], col = 1:2, ylab = "Survival", xlab = "Time in days", lwd = 2, mark.time = TRUE,
             main = paste0(round(qvals[j], 3), " (", round(pvals[j], 3), ")"))
        legend(max(sfits[[i]]$time), 0.9, c(paste0("wild type (",
                                                   sum(clinmut[, grep(paste0(i, "_mut"), colnames(clinmut))] == 0), ")"),
                                            paste0(i, " mutated\n(",
                                                   sum(clinmut[, grep(paste0(i, "_mut"), colnames(clinmut))] == 1)
                                                  ,")")),
               col=1:2, lwd=2, bty='n', xjust = 1, cex = 1)
    }
    dev.off()

    ## do survival for low expression:
    print(all(clinf[which(clinf$bcr_patient_barcode %in% colnames(dataf)), 1] == colnames(dataf)))
    dataf2 <- dataf
    dataf2[which(dataf2 != -1)] <- 0
    dataf2 <- abs(dataf2)
    clinexp <- cbind(clinf[which(clinf$bcr_patient_barcode %in% colnames(dataf2)), ], t(dataf2))
    if (apc) {
        clinexp <- clinexp[which(clinexp[, 1] %in% clinmutsave[which(clinmutsave[, "APC_mut"] == 1), 1]), ]
    }
    sfits <- list()
    coxs <- list()
    pvals <- numeric(length(all))
    for (j in 1:length(all)) {
        i <- all[j]
        coxs[[i]] <- coxph(Surv(days_to_last_follow_up, vital_status) ~ get(paste0(i, "_exp"))+tumor_stage+age_at_diagnosis, clinexp)
        pvals[j] <- summary(coxs[[i]])$coefficients[1, 5]
    }
    qvals <- p.adjust(pvals, method = "BH")

    if (giant) {
        pdf(paste0("Schwank/", type, "_exp.pdf"), width = 14, height = 7)
        par(mfrow=c(2,4))
    } else {
        pdf(paste0("Schwank/", type, "_exp.pdf"), width = 5, height = 5)
    }
    for (j in 1:length(all)) {
        i <- all[j]
        sfits[[i]] <- survfit(Surv(days_to_last_follow_up, vital_status) ~ get(paste0(i, "_exp")), clinexp)
        plot(sfits[[i]], col = 1:2, ylab = "Survival", xlab = "Time in days", lwd = 2, mark.time = TRUE,
             main = paste0(round(qvals[j], 3), " (", round(pvals[j], 3), ")"))
        legend(max(sfits[[i]]$time), 0.9, c(paste0("wild type (",
                                                   sum(clinexp[, grep(paste0(i, "_exp"), colnames(clinexp))] == 0), ")"),
                                            paste0(i, " low expression\n(",
                                                   sum(clinexp[, grep(paste0(i, "_exp"), colnames(clinexp))] == 1)
                                                  ,")")),
               col=1:2, lwd=2, bty='n', xjust = 1, cex = 1)
    }
    dev.off()

    ## ## do survival for low expression zscored:
    print(all(clinf[which(clinf$bcr_patient_barcode %in% colnames(datazf)), 1] == colnames(datazf)))
    datazf2 <- datazf
    datazf2[which(datazf2 != -1)] <- 0
    datazf2 <- abs(datazf2)
    clinexp <- cbind(clinf[which(clinf$bcr_patient_barcode %in% colnames(datazf2)), ], t(datazf2))
    if (apc) {
        clinexp <- clinexp[which(clinexp[, 1] %in% clinmutsave[which(clinmutsave[, "APC_mut"] == 1), 1]), ]
    }
    sfits <- list()
    coxs <- list()
    pvals <- numeric(length(all))
    for (j in 1:length(all)) {
        i <- all[j]
        coxs[[i]] <- coxph(Surv(days_to_last_follow_up, vital_status) ~ get(paste0(i, "_expz"))+tumor_stage+age_at_diagnosis, clinexp)
        pvals[j] <- summary(coxs[[i]])$coefficients[1, 5]
    }
    qvals <- p.adjust(pvals, method = "BH")
    if (giant) {
        pdf(paste0("Schwank/", type, "_exp_z.pdf"), width = 14, height = 7)
        par(mfrow=c(2,4))
    } else {
        pdf(paste0("Schwank/", type, "_exp_z.pdf"), width = 5, height = 5)
    }
    for (j in 1:length(all)) {
        i <- all[j]
        sfits[[i]] <- survfit(Surv(days_to_last_follow_up, vital_status) ~ get(paste0(i, "_expz")), clinexp)
        plot(sfits[[i]], col = 1:2, ylab = "Survival", xlab = "Time in days", lwd = 2, mark.time = TRUE,
             main = paste0(round(qvals[j], 3), " (", round(pvals[j], 3), ")"))
        legend(max(sfits[[i]]$time), 0.9, c(paste0("wild type (",
                                                   sum(clinexp[, grep(paste0(i, "_expz"), colnames(clinexp))] == 0), ")"),
                                            paste0(i, " low expression\n(",
                                                   sum(clinexp[, grep(paste0(i, "_expz"), colnames(clinexp))] == 1)
                                                  ,")")),
               col=1:2, lwd=2, bty='n', xjust = 1, cex = 1)
    }
    dev.off()

    ## ## do survival for copynumber var:
    print(all(clinf[which(clinf$bcr_patient_barcode %in% colnames(cnvf)), 1] == colnames(cnvf)))
    cnvf2 <- cnvf
    cnvf2[which(abs(cnvf2) <= 1)] <- 0
    cnvf2[which(abs(cnvf2) > 1)] <- 1
    clincnv <- cbind(clinf[which(clinf$bcr_patient_barcode %in% colnames(cnvf)), ], t(cnvf2))
    ## sfits <- list()
    ## coxs <- list()
    ## pvals <- numeric(length(all))
    ## for (j in 1:length(all)) {
    ##     i <- all[j]
    ##     coxs[[i]] <- coxph(Surv(days_to_last_follow_up, vital_status) ~ get(paste0(i, "_cnv"))+tumor_stage+age_at_diagnosis, clincnv)
    ##     pvals[j] <- summary(coxs[[i]])$coefficients[1, 5]
    ## }
    ## qvals <- p.adjust(pvals, method = "BH")
    ## if (giant) {
    ##     pdf(paste0("Schwank/", type, "_cnv.pdf"), width = 14, height = 7)
    ##     par(mfrow=c(2,4))
    ## } else {
    ##     pdf(paste0("Schwank/", type, "_cnv.pdf"), width = 5, height = 5)
    ## }
    ## for (j in 1:length(all)) {
    ##     i <- all[j]
    ##     sfits[[i]] <- survfit(Surv(days_to_last_follow_up, vital_status) ~ get(paste0(i, "_cnv")), clincnv)
    ##     plot(sfits[[i]], col = 1:2, ylab = "Survival", xlab = "Time in days", lwd = 2, mark.time = TRUE,
    ##          main = paste0(round(qvals[j], 3), " (", round(pvals[j], 3), ")"))
    ##     legend(max(sfits[[i]]$time), 0.9, c(paste0("wild type (",
    ##                                                sum(clincnv[, grep(paste0(i, "_cnv"), colnames(clincnv))] == 0), ")"),
    ##                                         paste0(i, " low expression\n(",
    ##                                                sum(clincnv[, grep(paste0(i, "_cnv"), colnames(clincnv))] == 1)
    ##                                               ,")")),
    ##            col=1:2, lwd=2, bty='n', xjust = 1, cex = 1)
    ## }
    ## dev.off()

    ## ## do survival for copynumber var plus mutated:
    ## ## combine (bad practice?):
    ## clinmutcnv <- clinmut
    ## colnames(clinmutcnv)[grep("_mut", colnames(clinmutcnv))] <- gsub("_mut", "_cnv", colnames(clinmutcnv)[grep("_mut", colnames(clinmutcnv))])
    ## clinmutcnv <- rbind(clinmutcnv,
    ##                     clincnv[which(!(clincnv[, 1] %in% clinmut[, 1])), ])
    ## clinmutcnv <- clinmutcnv[order(clinmutcnv[, 1]), ]
    ## clinmutcnv[, grep("_cnv", colnames(clinmutcnv))] <- 0
    ## clinmutcnv[which(clinmutcnv[, 1] %in% clinmut[, 1]), grep("_cnv", colnames(clinmutcnv))] <-
    ##     clinmut[which(clinmut[, 1] %in% clinmutcnv[, 1]), grep("_mut", colnames(clinmut))]
    ## clinmutcnv[which(clinmutcnv[, 1] %in% clincnv[, 1]), grep("_cnv", colnames(clinmutcnv))] <-
    ##     clincnv[which(clincnv[, 1] %in% clinmutcnv[, 1]), grep("_cnv", colnames(clincnv))] + clinmutcnv[which(clinmutcnv[, 1] %in% clincnv[, 1]), grep("_cnv", colnames(clinmutcnv))]
    ## clinmutcnv[, grep("_cnv", colnames(clinmutcnv))] <- apply(clinmutcnv[, grep("_cnv", colnames(clinmutcnv))],
    ##                                                           c(1,2), function(x) {
    ##                                                               y <- as.numeric(x)
    ##                                                               if (y >= 1) {
    ##                                                                   y <- 1
    ##                                                               } else {
    ##                                                                   y <- 0
    ##                                                               }
    ##                                                               return(y)
    ##                                                           })
    ## clinmutcnvcomb <- clinmutcnv
    ## intersect:
    clinmutcnv <- clinmut[which(clinmut[, 1] %in% clincnv[, 1]), ]
    colnames(clinmutcnv)[grep("_mut", colnames(clinmutcnv))] <- gsub("_mut", "_cnv", colnames(clinmutcnv)[grep("_mut", colnames(clinmutcnv))])
    clinmutcnv <- clinmutcnv[order(clinmutcnv[, 1]), ]
    clinmutcnv[, grep("_cnv", colnames(clinmutcnv))] <- 0
    clinmutcnv[which(clinmutcnv[, 1] %in% clinmut[, 1]), grep("_cnv", colnames(clinmutcnv))] <-
        clinmut[which(clinmut[, 1] %in% clinmutcnv[, 1]), grep("_mut", colnames(clinmut))]
    clinmutcnv[which(clinmutcnv[, 1] %in% clincnv[, 1]), grep("_cnv", colnames(clinmutcnv))] <-
        clincnv[which(clincnv[, 1] %in% clinmutcnv[, 1]), grep("_cnv", colnames(clincnv))] + clinmutcnv[which(clinmutcnv[, 1] %in% clincnv[, 1]), grep("_cnv", colnames(clinmutcnv))]
    clinmutcnv[, grep("_cnv", colnames(clinmutcnv))] <- apply(clinmutcnv[, grep("_cnv", colnames(clinmutcnv))],
                                                              c(1,2), function(x) {
                                                                  y <- as.numeric(x)
                                                                  if (y >= 1) {
                                                                      y <- 1
                                                                  } else {
                                                                      y <- 0
                                                                  }
                                                                  return(y)
                                                              })
    clinmutcnvint <- clinmutcnv
    ## choose:
    ##clinmutcnv <- clinmutcnvcomb
    if (apc) {
        clinmutcnv <- clinmutcnv[which(clinmutcnv[, 1] %in% clinmutsave[which(clinmutsave[, "APC_mut"] == 1), 1]), ]
    }
    sfits <- list()
    coxs <- list()
    pvals <- numeric(length(all))
    for (j in 1:length(all)) {
        i <- all[j]
        coxs[[i]] <- coxph(Surv(days_to_last_follow_up, vital_status) ~ get(paste0(i, "_cnv"))+tumor_stage+age_at_diagnosis, clinmutcnv)
        pvals[j] <- summary(coxs[[i]])$coefficients[1, 5]
    }
    qvals <- p.adjust(pvals, method = "BH")
    if (giant) {
        pdf(paste0("Schwank/", type, "_cnv_mut.pdf"), width = 14, height = 7)
        par(mfrow=c(2,4))
    } else {
        pdf(paste0("Schwank/", type, "_cnv_mut.pdf"), width = 5, height = 5)
    }
    for (j in 1:length(all)) {
        i <- all[j]
        sfits[[i]] <- survfit(Surv(days_to_last_follow_up, vital_status) ~ get(paste0(i, "_cnv")), clinmutcnv)
        plot(sfits[[i]], col = 1:2, ylab = "Survival", xlab = "Time in days", lwd = 2, mark.time = TRUE,
             main = paste0(round(qvals[j], 3), " (", round(pvals[j], 3), ")"))
        legend(max(sfits[[i]]$time), 0.9, c(paste0("wild type (",
                                                   sum(clinmutcnv[, grep(paste0(i, "_cnv"), colnames(clinmutcnv))] == 0), ")"),
                                            paste0(i, " mutated\n(",
                                                   sum(clinmutcnv[, grep(paste0(i, "_cnv"), colnames(clinmutcnv))] == 1)
                                                  ,")")),
               col=1:2, lwd=2, bty='n', xjust = 1, cex = 1)
    }
    dev.off()

}

stop("done")

source("testing/vignettes/schwankII.r")

### other:

source("~/Documents/mnem/R/mnems.r")
source("~/Documents/mnem/R/mnems_low.r")
library(naturalsort)

for (type in c("COAD", "READ")) {

    R <- D <- t(get(paste0("mut.", type)))

    R[which(R == 0)] <- -1

    best <- nemEst(R)

    nemEst.bootstrap <- function(R, num, pct, ...) {
        best <- nemEst(R, ...)
        phisum <- best$phi*0
        thetasum <- best$theta
        for (i in 1:num) {
            subset <- sample(1:nrow(R), ceiling(nrow(R)*pct), replace = TRUE)
            Rsub <- R[subset, ]
            tmp <- nemEst(Rsub, ...)
            phisum <- phisum + tmp$phi # transitive.reduction(tmp$phi)
            thetasum[, subset] <- thetasum[, subset] + tmp$theta
        }
        phisum <- phisum/num
        thetasum <- thetasum/num
        return(list(phi = phisum, theta = thetasum, best = best))
    }

    boot <- nemEst.bootstrap(R, 1000, 0.8)

    ## compute exact me matrix:

    MUE <- t(D)%*%D
    MUE[which(MUE > 0)] <- 1
    MUE <- MUE[order(rownames(MUE)), order(colnames(MUE))]

    cutoff <- 0.5
    cutoff2 <- 0.05

    sup <- boot$phi
    sup[which(boot$phi > cutoff)] <- 1
    sup[which(boot$phi <= cutoff)] <- 0
    diag(sup) <- 0

    sup3 <- boot$phi*0
    sup3[which(boot$phi/t(boot$phi) > 0.5 & boot$phi/t(boot$phi) < 2 & boot$phi+t(boot$phi) > cutoff)] <- 1
    diag(sup3) <- 0

    sup2 <- boot$phi
    sup2[which(boot$phi <= cutoff2)] <- 1
    sup2[which(boot$phi > cutoff2)] <- 0
    diag(sup2) <- 0
    sup2 <- sup2*t(sup2)
    sup2 <- sup2*(MUE)

    dnf4 <- adj2dnf(1-MUE)

    ## apply(boot$theta, 1, sum) # E-gene attachement: how to use it?

    trsup <- sup#transitive.reduction(sup)
    dnf <- adj2dnf(trsup)
    freq <- trsup*boot$phi
    freq <- c(rep(0, length(dnf) - length(grep("=", dnf))), t(freq)[which(t(freq) != 0)])

    trsup2 <- sup2
    dnf2 <- adj2dnf(trsup2)
    freq2 <- trsup2*(1-boot$phi)
    diag(freq2) <- 0
    freq2 <- c(rep(0, length(dnf2) - length(grep("=", dnf2))), t(freq2)[which(t(freq2) != 0)])

    trsup3 <- sup3
    dnf3 <- adj2dnf(trsup3)
    freq3 <- trsup3*(1-boot$phi)
    diag(freq3) <- 0
    freq3 <- c(rep(0, length(dnf3) - length(grep("=", dnf3))), t(freq2)[which(t(freq3) != 0)])

    pdf(paste0("Schwank/", type, "_ME_CO.pdf"), width = 16, height = 16)

    par(mfrow=c(3,2))
    if (length(grep("=", dnf)) != 0) {
        plotDnf(dnf, freq = freq, edgelabel = freq, main = paste0("bootstrap support > ", round(cutoff*100), "% for one sided co-occurence"))
    }
    if (length(grep("=", dnf2)) != 0) {
        plotDnf(dnf2, freq = freq2, edgelabel = freq2, main = paste0("bootstrap support > ", round((1-cutoff2)*100), "% for mutual exclusivity"))
    }
    plotDnf(adj2dnf(transitive.reduction(boot$best$phi)), main = "MLE")
    if (length(grep("=", dnf4)) != 0) {
        plotDnf(dnf4, main = "Exact mutually exclusive")
    }
    if (length(grep("=", dnf3)) != 0) {
        plotDnf(dnf3, main = "Co-occurence")
    }

    dev.off()

}

###### DEBUGGING:

runs <- 1000

n <- 100

dat <- rnorm(n)

cs <- matrix(0, n, runs)

for (i in 1:runs) {

    tmp <- kmeans(dat, c(min(dat), median(data), max(dat)), nstart = 100)
    low <- which.min(tmp$centers)
    low <- which(tmp$cluster == low)
    high <- which.max(tmp$centers)
    high <- which(tmp$cluster == high)
    tmp2 <- tmp$cluster*0 +2
    tmp2[low] <- 1
    tmp2[high] <- 3

    cs[, i] <- tmp2

}

HeatmapOP(cor(cs), breaks = 100)

HeatmapOP(cs, colv = F, Rowv = F, breaks = 3)

csun <- apply(cs, 2, function(x) return(paste0(x, collapse = "")))

cs2 <- cs[, -which(duplicated(csun) == TRUE)]

HeatmapOP(cor(cs2), breaks = 100)

HeatmapOP(cs2, colv = F, Rowv = F, breaks = 3)

## combine read and coad:

load("COAD_clin.rda")
load("COAD_cnv.rda")
load("COAD_data.rda")
load("COAD_mut.rda")
load("COAD_biomart.rda")

clinical0 <- clinical
cnv0 <- cnv
data0 <- data
M10 <- M1
M20 <- M2
M30 <- M3
M40 <- M4

load("READ_clin.rda")
load("READ_cnv.rda")
load("READ_data.rda")
load("READ_mut.rda")
load("READ_biomart.rda")

clinical1 <- clinical
cnv1 <- cnv
data1 <- data
M11 <- M1
M21 <- M2
M31 <- M3
M41 <- M4

clinical <- rbind(clinical0, clinical1)
cnv <- cbind(cnv0, cnv1)
data <- cbind(data0, data1)
M1 <- cbind(M10[which(rownames(M10) %in% rownames(M11)), ], M11[which(rownames(M11) %in% rownames(M10)), ])
M2 <- cbind(M20[which(rownames(M20) %in% rownames(M21)), ], M21[which(rownames(M21) %in% rownames(M20)), ])
M3 <- cbind(M30[which(rownames(M30) %in% rownames(M31)), ], M31[which(rownames(M31) %in% rownames(M30)), ])
M4 <- cbind(M40[which(rownames(M40) %in% rownames(M41)), ], M41[which(rownames(M41) %in% rownames(M40)), ])

save(clinical, file = "both_clin.rda")
save(cnv, file = "both_cnv.rda")
save(data, file = "both_data.rda")
save(M1, M2, M3, M4, file = "both_mut.rda")
save(gene.symbols, file = "both_biomart.rda")

## do dawnrank:

load("COAD_dawnrank.rda")

data.full <- cbind(data, data.normal)

nf <- calcNormFactors(data.full)

dataDR <- t(t(data.full)*nf)

library(DawnRank)

data(goldStandard)

dataDR <- DawnNormalize(tumorMat = dataDR[, 1:ncol(data)], normalMat = data.full[, (ncol(data)+1):ncol(data.full)])

colnames(dataDR) <- unlist(lapply(colnames(dataDR), function(x) {
    y <- unlist(strsplit(x, "-"))
    y <- y[1:3]
    y <- paste(y, collapse = "-")
    return(y)
}))

dataDR2 <- dataDR[which(rownames(dataDR) %in% rownames(mut.mat)), ]
mut2 <- mut.mat[which(rownames(mut.mat) %in% rownames(dataDR)), ]

dataDR2 <- dataDR2[which(rownames(dataDR2) %in% rownames(adj2)), ]
mut3 <- mut2[which(rownames(mut2) %in% rownames(adj2)), ]
adj3 <- adj2[which(rownames(adj2) %in% rownames(mut3)), which(colnames(adj2) %in% rownames(mut2))]

dataDR2 <- dataDR2[, which(colnames(dataDR2) %in% colnames(mut3))]
mut3 <- mut3[, which(colnames(mut3) %in% colnames(dataDR2))]

dataDR3 <- matrix(0, nrow(mut3), ncol(dataDR2))
rownames(dataDR3) <- unique(rownames(dataDR2))
colnames(dataDR3) <- colnames(dataDR2)

for (i in 1:length(unique(rownames(dataDR2)))) {
    gene <- unique(rownames(dataDR2))[i]
    tmp <- dataDR2[which(rownames(dataDR2) %in% gene), ]
    if (sum(rownames(dataDR2) %in% gene) > 1) {
        dataDR3[i, ] <- apply(tmp, 2, median)
    } else {
        dataDR3[i, ] <- tmp
    }
}

dataDR2 <- dataDR3

dataDR3 <- matrix(0, nrow(dataDR2), ncol(mut3))
rownames(dataDR3) <- rownames(dataDR2)
colnames(dataDR3) <- unique(colnames(dataDR2))

for (i in 1:length(unique(colnames(dataDR2)))) {
    sample <- unique(colnames(dataDR2))[i]
    tmp <- dataDR2[, which(colnames(dataDR2) %in% sample)]
    if (sum(colnames(dataDR2) %in% sample) > 1) {
        dataDR3[, i] <- apply(tmp, 1, median)
    } else {
        dataDR3[, i] <- tmp
    }
}

dataDR3 <- dataDR3[order(rownames(dataDR3)), order(colnames(dataDR3))]
mut3 <- mut3[order(rownames(mut3)), order(colnames(mut3))]
adj3 <- adj3[order(rownames(adj3)), order(colnames(adj3))]

dawnRankScore <- DawnRank(adjMatrix = adj3, mutationMatrix = mut3, expressionMatrix = dataDR3, mu = 3, goldStandard = goldStandard, parallel = 4)

dawnRankFrame <- dawnRankScore[[3]]
head(dawnRankFrame)

aggregateDawnRankScore <- condorcetRanking(scoreMatrix = dawnRankScore[[2]], mutationMatrix = mut.mat, parallel = 4)

top10 <- aggregateDawnRankScore[[2]][1:10]
top10

dawnRankFrame$isCGC <- dawnRankFrame$isGoldStandard

significance <- patspeccutoff(patient = "TCGA-A2-A04P", ms = dawnRankFrame, default = 95)

significance[[1]]

