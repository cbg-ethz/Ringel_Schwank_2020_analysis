
library(DawnRank)
library(STRINGdb)

for (type in c("COAD", "READ")) {

    ## mutations:

    ## type <- "COAD"

    ## mutation1 <- GDCquery_Maf(tumor = type, save.csv = TRUE, pipeline = "varscan2") # mutation file
    
    ## mutation2 <- GDCquery_Maf(tumor = type, save.csv = TRUE, pipeline = "muse") # mutation file
    
    ## mutation3 <- GDCquery_Maf(tumor = type, save.csv = TRUE, pipeline = "somaticsniper") # mutation file
    
    ## mutation4 <- GDCquery_Maf(tumor = type, save.csv = TRUE, pipeline = "mutect2") # mutation file

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
    
    mutation <- rbind(mutation1, mutation2, mutation3, mutation4)

    #mutation <- mutation[which(mutation$Hugo_Symbol %in% all), ]

    mut.mat <- matrix(0, length(unique(mutation$Hugo_Symbol)), length(unique(mutation$Tumor_Sample_Barcode)))

    colnames(mut.mat) <- sort(unique(mutation$Tumor_Sample_Barcode))

    rownames(mut.mat) <- sort(unique(mutation$Hugo_Symbol))

    coln <- which(colnames(mutation) %in% "Tumor_Sample_Barcode")

    library(snowfall)

    sfInit(parallel = T, cpus = 4)

    sfExport("mutation", "mut.mat", "coln")

    countSamples <- function(i, genes, mut.mat) {
        
        i <- which(rownames(mut.mat) %in% genes[i])
        
        tmp <- mutation[which(mutation$Hugo_Symbol %in% rownames(mut.mat)[i]), coln]
        
        tmp2 <- mut.mat[i, ]
        
        tmp2[which(colnames(mut.mat) %in% tmp)] <- 1
        
        return(tmp2)
        
    }
    
    M <- sfLapply(as.list(1:nrow(mut.mat)), countSamples, colnames(mut.mat), mut.mat)

    mut.mat <- do.call("rbind", M)

    colnames(mut.mat) <- sort(unique(mutation$Tumor_Sample_Barcode))

    rownames(mut.mat) <- all

    sfStop()

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

    assign(paste("query.", type, ".exprs", sep = ""), GDCquery(project = paste("TCGA-", type, sep = ""), sample.type = "Primary solid Tumor",
                                                                     data.category = "Transcriptome Profiling"
                                                                   , data.type = "Gene Expression Quantification"
                                                                   , workflow.type = "HTSeq - Counts"
                                                                     )
           )
    
    GDCdownload(get(paste("query.", type, ".exprs", sep = "")))

    assign(paste("exprs.", type, sep = ""), GDCprepare(get(paste("query.", type, ".exprs", sep = ""))))

    data <- assay(get(paste("exprs.", type, sep = "")))

    mart = useMart("ensembl")

    mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

    tmp <- listAttributes(mart)

    gene.symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position", "chromosome_name"), mart = mart, values = rownames(data), filters = "ensembl_gene_id")

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

    datar <- log(datar+0.1)

    for (i in all) {
        ## k-means:
        d <- dist(datar[i, ])
        bestk <- 3
        ## bestsc <- 0
        ## for (k in 2:10) {
        ##     tmp <- kmeans(d, k)
        ##     sc <- silhouette(tmp$cluster, d)
        ##     sc <- sum(sc[, 3])
        ##     if (sc > bestsc) {
        ##         bestsc <- sc
        ##         bestk <- k
        ##     }
        ##     bestk <- 3 # force 3 clusters
        ## }
        ## print(bestk)
        tmp <- kmeans(datar[i, ], bestk)
        low <- which.min(tmp$centers)
        low <- which(tmp$cluster == low)
        high <- which.max(tmp$centers)
        high <- which(tmp$cluster == high)
        tmp <- datar[i, ]*0
        tmp[low] <- -1
        tmp[high] <- 1
        datar[i, ] <- tmp
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

}

stop("dawnrank")

## colorectal cancer

## CRISPR screen from organoids identifies several interesting genes.

## TCGA analysis:

## possible focus on subtypes: stable (low mut. freq., APC, annotated?); (early stage tumors?)

Sgenes <- c("ARID1A", "ARID2", "SMARCA4", "SMARCB1", "CNIH4", "KEAP1", "FAM122A", "NBAS")

core <- c("TGFBR1", "TGFBR2", "SMAD4", "SMAD3")

all <- c(Sgenes, core, "APC")

chrom <- c(1,12,19,22,1,19,9,2,9,3,18,15)

strand <- c(1,1,1,1,1,0,1,0,1,1,1,1)

start <- c(26693236,45729665,10960922,23786963,224356811,10486120,68780034,99104038,30606493,51028394,67063763)

end <- c(26782110,4590800,11065395,23834518,10503741,68784608,99154192,30694142,51085045,67195195)

library("cgdsr")

mycgds <- CGDS("http://www.cbioportal.org/public-portal/")

tmp <- getCancerStudies(mycgds)

study <- tmp[grep(paste0("^", tolower(type), "_"), tmp[, 1]), 1]

getGeneticProfiles(mycgds,study)[,c(1:2)]

getCaseLists(mycgds,study)[,c(1:2)]

cnv <- getProfileData(mycgds, all, paste0(tolower(type), "_tcga_gistic"), paste0(tolower(type), "_tcga_all"))

rownames(cnv) <- gsub(".01", "", tolower(rownames(cnv)))

# Get clinical data for the case list
clindata <- getClinicalData(mycgds,paste0(tolower(type), "_tcga_all"))

clindata$VITAL_STATUS[which(clindata$VITAL_STATUS %in% "Alive")] <- 0

clindata$VITAL_STATUS[which(clindata$VITAL_STATUS %in% "Dead")] <- 1

clindata$VITAL_STATUS[which(clindata$VITAL_STATUS %in% "")] <- NA

clindata$VITAL_STATUS <- as.numeric(clindata$VITAL_STATUS)

cnv <- apply(cnv, c(1,2), function(x) {
    if (is.nan(x) == TRUE) { x <- NA
    } else if (x %in% c(1,2,-1,-2)) { x <- 1
    } else { x <- x }
})

gene <- "CNIH4"

sfit <- survfit(Surv(DAYS_TO_LAST_FOLLOWUP, VITAL_STATUS) ~ get(gene), cbind(clindata, cnv)); plot(sfit); tmp <- coxph(Surv(DAYS_TO_LAST_FOLLOWUP, VITAL_STATUS) ~ get(gene), cbind(clindata, cnv))

## 1) look for mutation rates and downregulation:

## 2) check survival for mutation and downregulation:

## 3) mutual exclusivity with tgfbeta (targets?); possibly pool genes to increase power; additive effect with "soft" mutual exclusivity:


## 1) Check for sig. driver gene: done, calculate dawnrank pvalue

types <- c("COAD", "READ")

for (i in 1:length(types)) {

    ## i <- 2

    tmp <- read.csv(paste0("DawnRank_", types[i], "_genelist.csv"))

    if (i == 1) {
        ## plot(density(tmp[, 2]), pch = as.character(i), type = "b", ylim = c(0,3.5))
        hist(tmp[, 2], col = rgb(1,0,0,0.5), freq=FALSE)
    } else {
        ## lines(density(tmp[, 2]), col = i, pch = as.character(i), type = "b")
        hist(tmp[, 2], add = TRUE, col = rgb(0,1,0,0.5), freq = FALSE)
    }
    
    a <- which(tmp[, 1] %in% Sgenes)

    print(wilcox.test(a, 1:nrow(tmp)))
    
    print(tmp[a, ])

    b <- which(tmp[, 1] %in% core)

    print(wilcox.test(b, 1:nrow(tmp)))
    
    print(tmp[b, ])
    
}

abline(v=0.5)

## 2) check survival:

library(RTCGAToolbox)

stddata <- getFirehoseRunningDates()
stddata

for (type in c("COAD", "READ", "COADREAD", "STAD")) {

    ## type <- "COAD"

    data <- getFirehoseData(dataset=type, runDate="20160128", RNASeqGene = TRUE, RNASeq2GeneNorm = TRUE,
                            forceDownload=FALSE, clinical=TRUE, Mutation=TRUE, CNASNP = TRUE,
                            CNASeq = TRUE, Methylation = TRUE)

    mut <- table(data@Mutation$Hugo_Symbol)

    mut[Sgenes]

    ## ARID1A   ARID2 SMARCA4    <NA>    <NA>    <NA>    <NA>    NBAS 
    ##      5       8       3                                       4

    mut[core]

    ## TGFBR1   <NA>  SMAD4  SMAD3 
    ##      3             8      3 

    library(survival)

    clin <- as.data.frame(data@clinical)

    cin <- clin[order(rownames(clin)), ]

    mut <- data@Mutation

    tmp <- mut$Tumor_Sample_Barcode

    tmp <- unlist(lapply(tmp, function(x) {
        x <- unlist(strsplit(x, "-"))
        x <- paste(x[1:3], collapse = ".")
        x <- tolower(x)
        return(x)
    }))

    mut <- cbind(mut, id = tmp)

    cnasnp <- data@CNASNP

    cnasnp[, 1] <- unlist(lapply(cnasnp[, 1], function(x) {
        x <- unlist(strsplit(x, "-"))
        x <- paste(x[1:3], collapse = ".")
        x <- tolower(x)
        return(x)
    }))

    cnasnp$Segment_Mean <- scale(cnasnp$Segment_Mean) 

    for (i in all) {
        ## add mutations:
        clin <- cbind(clin, 0)
        colnames(clin)[ncol(clin)] <- paste0(i, "_mut")
        ids <- mut$id[which(mut$Hugo_Symbol %in% i)]
        clin[which(rownames(clin) %in% ids), paste0(i, "_mut")] <- 1
        ## add cnasnp:
        clin <- cbind(clin, 0)
        colnames(clin)[ncol(clin)] <- paste0(i, "_cnasnp")
        
        ids <- cnasnp$Sample[which(cnasnp$Start <= start[which(all %in% i)] & cnasnp$End >= end[which(all %in% i)] & cnasnp$Chromosome == chrom[which(all %in% i)])]
        tmp <- cnasnp$Segment_Mean[which(cnasnp$Start <= start[which(all %in% i)] & cnasnp$End >= end[which(all %in% i)] & cnasnp$Chromosome == chrom[which(all %in% i)])]
        tmp <- tmp[which(duplicated(ids) == FALSE)]
        ids <- ids[which(duplicated(ids) == FALSE)]
        tmp <- tmp[order(ids)]
        ids <- ids[order(ids)]
        tmp <- tmp[which(ids %in% rownames(clin))]
        ids <- ids[which(ids %in% rownames(clin))]
        clin[which(rownames(clin) %in% ids), paste0(i, "_cnasnp")] <- tmp
    }
    
    clin$days_to_last_followup <- as.numeric(clin$days_to_last_followup)
    
    clin$days_to_death <- as.numeric(clin$days_to_death)

    clin <- cbind(clin, event = clin$days_to_last_followup)

    clin$event[which(is.na(clin$days_to_death) == FALSE)] <- clin$days_to_death[which(is.na(clin$days_to_death) == FALSE)]
    
    clin$vital_status <- as.numeric(clin$vital_status)

    stages <- sort(unique(clin$pathologic_stage))
    
    for (i in 1:length(stages)) {
        
        clin$pathologic_stage <- gsub(paste0("^", stages[i], "$"), i, clin$pathologic_stage)

    }
    
    clin$pathologic_stage <- as.numeric(clin$pathologic_stage)

    ## sfit <- survfit(Surv(event, vital_status) ~ ARID1A_cnasnp, clin); plot(sfit);  tmp <- coxph(Surv(event, vital_status) ~ ARID1A_cnasnp, clin)

    sfits <- list()
    coxs <- list()

    pdf(paste0("Schwank_Firehose_", type, ".pdf"), width = 8, height = 8)

    for (i in c(Sgenes, core)) {
        coxs[[i]] <- coxph(Surv(event, vital_status) ~ get(paste0(i, "_mut"))+pathologic_stage, clin)
        sfits[[i]] <- survfit(Surv(event, vital_status) ~ get(paste0(i, "_mut")), clin)
        plot(sfits[[i]], col = 1:2, ylab = "Survival", xlab = "Time in days", lwd = 2,
             main = round(summary(coxs[[i]])$coefficients[1, 5], 3))
        legend(max(sfits[[i]]$time), 0.9, c("wild type", paste0(i, " mutation")),
               col=1:2, lwd=2, bty='n', xjust = 1, cex = 1)
    }

    dev.off()

    pvals <- numeric(length(Sgenes))
    for (i in 1:length(Sgenes)) {
        if (length(coxs[[i]]$wald.test) == 0) {
            pvals[i] <- NA
        } else {
            pvals[i] <- summary(coxs[[i]])$logtest[3]
        }
    }

    pvals <- pvals[!is.na(pvals)]

    x <- -2*sum(log(pvals))

    1 - pchisq(x, 2*length(pvals))

### gene expression:

    norm <- data@RNASeq2GeneNorm

    lnorm <- log2(norm + 0.1)

    znorm <- t(scale(t(lnorm)))

    colnames(znorm) <- unlist(lapply(colnames(znorm), function(x) {
        x <- unlist(strsplit(x, "-"))
        x <- paste(x[1:3], collapse = ".")
        x <- tolower(x)
        return(x)
    }))

    znorm <- znorm[, order(colnames(znorm))]

    clin <- clin[order(rownames(clin)), ]

    clin2 <- clin[which(rownames(clin) %in% colnames(znorm)), ]

    znorm <- znorm[, which(colnames(znorm) %in% rownames(clin2))]

    dups <- which(duplicated(colnames(znorm)) == TRUE)

    for (i in unique(colnames(znorm)[dups])) {
        tmp <- znorm[, which(colnames(znorm) %in% i)]
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
        znorm[, which(colnames(znorm) %in% i)[1]] <- tmp2
        znorm <- znorm[, -which(colnames(znorm) %in% i)[-1]]
    }
    
    library(cluster)
    
    for (i in c(Sgenes, core)) {
        ## ## hard cut at 2 sd:
        ## tmp <- znorm[i, ]
        ## low <- which(tmp < -2)
        ## tmp[low] <- 1
        ## tmp[-low] <- 0
        ## k-means:
        d <- dist(znorm[i, ])
        bestk <- 2
        bestsc <- 0
        for (k in 2:10) {
            tmp <- kmeans(d, k)
            sc <- silhouette(tmp$cluster, d)
            sc <- sum(sc[, 3])
            if (sc > bestsc) {
                bestsc <- sc
                bestk <- k
            }
            bestk <- 3 # force 3 clusters
        }
        print(bestk)
        tmp <- kmeans(znorm[i, ], bestk)
        low <- which.min(tmp$centers)
        low <- which(tmp$cluster == low)
        tmp <- znorm[i, ]
        tmp[low] <- 1
        tmp[-low] <- 0
        clin2 <- cbind(clin2, tmp)
        colnames(clin2)[ncol(clin2)] <- paste0(i, "_low")
        clin2 <- cbind(clin2, znorm[i, ])
        colnames(clin2)[ncol(clin2)] <- paste0(i, "_cont")
    }

    sfits <- list()
    coxs <- list()

    pdf(paste0("Schwank_Firehose_", type, "_rnaseq.pdf"), width = 8, height = 8)

    for (i in c(Sgenes, core)) {
        
        coxs[[i]] <- coxph(Surv(event, vital_status) ~ get(paste0(i, "_low"))+pathologic_stage, clin2)
        sfits[[i]] <- survfit(Surv(event, vital_status) ~ get(paste0(i, "_low")), clin2)
        plot(sfits[[i]], col = 1:2, ylab = "Survival", xlab = "Time in days", lwd = 2,
             main = round(summary(coxs[[i]])$coefficients[1, 5], 3))
        legend(max(sfits[[i]]$time), 0.9, c("wild type", paste0(i, " low expression")),
               col=1:2, lwd=2, bty='n', xjust = 1, cex = 1)
    }

    dev.off()

    sfits2 <- list()
    coxs2 <- list()

    pdf(paste0("Schwank_Firehose_", type, "_rnaseq_cont.pdf"), width = 8, height = 8)

    for (i in c(Sgenes, core)) {
        
        coxs2[[i]] <- coxph(Surv(event, vital_status) ~ get(paste0(i, "_cont"))+pathologic_stage, clin2)
        sfits2[[i]] <- survfit(Surv(event, vital_status) ~ get(paste0(i, "_cont")), clin2)
        plot(sfits2[[i]], col = 1:2, ylab = "Survival", xlab = "Time in days", lwd = 2,
             main = round(summary(coxs2[[i]])$coefficients[1, 5], 3))
        legend(max(sfits[[i]]$time), 0.9, c("wild type", paste0(i, " low expression")),
               col=1:2, lwd=2, bty='n', xjust = 1, cex = 1)
    }

    dev.off()

    all <- c(Sgenes, core)

    mut <- mut2 <- joint <- list()

    for (i in 1:length(all)) {
        for (j in 1:length(all)) {
            if (i <= j) { next() }
            A <- paste0(all[i], "_mut")
            B <- paste0(all[j], "_mut")
            a <- sum(clin[, A] == 1 & clin[, B] == 1)
            b <- sum(clin[, A] == 1 & clin[, B] == 0)
            c <- sum(clin[, A] == 0 & clin[, B] == 1)
            d <- sum(clin[, A] == 0 & clin[, B] == 0)
            mut[[paste0(all[i], "_", all[j])]] <- phyper(a, a+c, b+d, a+b, lower.tail = TRUE)
            mut2[[paste0(all[i], "_", all[j])]] <- phyper(a, a+c, b, a+b, lower.tail = TRUE)
            joint[[paste0(all[i], "_", all[j])]] <- a
        }
    }

    mut <- c(unlist(mut), unlist(mut2))

    if (any(is.na(mut)) | (min(mut) == max(mut) & min(mut) == 1)) { next() }

    joint <- c(unlist(joint), unlist(joint))
    
    sig <- which(mut < 0.5)

    ## pdf(paste0("Schwank_Firehose_", type, "_mut_ex.pdf"), width = 8, height = 8)
    
    ## plot(mut, pch = 4, col = "red")
    ## if (length(sig) > 0) {
    ##     text(sig, mut[sig], paste0(names(mut)[sig], " ", joint[sig]), col = rgb(0,0,1,0.5), cex = 0.75)
    ## }
    ## abline(v=length(joint)/2 + 0.5, col = "black", lty = 2)
    
    ## dev.off()
    
    mm <- clin2[, grep("_mut", colnames(clin2))]
    
    pdf(paste0("Schwank_Firehose_", type, "_heatmap.pdf"), width = 8, height = 8)

    print(HeatmapOP(mm[order(apply(mm, 1, sum)), order(apply(mm, 2, sum))], Colv = FALSE, Rowv = FALSE, bordercol = "transparent"))

    dev.off()

}

stop("done")

## try cbioportal:

library("cgdsr")

mycgds = CGDS("http://www.cbioportal.org/public-portal/")

test(mycgds)

# Get list of cancer studies at server
studies <- getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[grep("stad_", studies[, 1]),1]
mycaselist = getCaseLists(mycgds,mycancerstudy[4])

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy[4])

# Get data slices for a specified list of genes, genetic profile and case list
profiledata <- getProfileData(mycgds,all,mycaselist[5, 1], mycaselist[3, 1])

cnv <- getProfileData(mycgds, all, c("stad_tcga_gistic"), "stad_tcga_all")

rownames(cnv) <- gsub(".01", "", tolower(rownames(cnv)))

# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,mycaselist$case_list_id[1],unlist(strsplit(mycaselist$case_ids, " ")))

# documentation
help('cgdsr')
help('CGDS')

####

source("testing/vignettes/schwank.r")

## get "logodds":

lods <- znorm # t(scale(t(abs(znorm))))

mutdata <- clin[, grep("mut", colnames(clin))]

index <- which(apply(mutdata, 1, sum) >= 1)

lods <- lods[, which(colnames(lods) %in% rownames(clin)[index])]

mutdata <- mutdata[which(rownames(mutdata) %in% colnames(lods)), ]

colnames(lods) <- unlist(lapply(colnames(lods), function(x) {
    y <- paste(gsub("_mut", "", colnames(mutdata)[which(mutdata[which(rownames(mutdata) %in% x), ] == 1)]), collapse = "_")
    return(y)
    }))

source("~/Documents/mnem/R/mnems_low.r")
library(naturalsort)
library(nem)
library(Rgraphviz)

lods[is.na(lods)] <- 0

res <- nemEst(lods)

plotDnf(adj2dnf(transitive.reduction(res$phi)))

