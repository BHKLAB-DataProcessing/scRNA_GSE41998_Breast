library(MultiAssayExperiment)
library(SummarizedExperiment)
library(dplyr)
library(qs)
library(biomaRt)
library(purrr)
library(GEOquery)

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[[1]]
filename <- args[[2]]

print('get study')

data <- getGEO("GSE41998", GSEMatrix = TRUE)[[1]]
assay_data <- makeSummarizedExperimentFromExpressionSet(data)
coldata <- data@phenoData@data
drugs <- sapply(strsplit(coldata$characteristics_ch1.3, ": "), "[[", 2)
coldata$characteristics_ch1.3 <- drugs
coldata <- coldata[coldata$characteristics_ch1.3 == "Paclitaxel", ]
assays <- as.data.frame(assay(assay_data))
new_coldata <- data.frame(patientid = strsplit(coldata$characteristics_ch1, ": ") %>% map_chr(`[`, 2))
rownames(coldata) <- new_coldata$patientid
colnames(assays) <- new_coldata$patientid
new_coldata$treatmentid <- "Paclitaxel"
new_coldata$tissueid <- "Breast"
new_coldata$response <- strsplit(coldata$characteristics_ch1.4, ": ") %>% map_chr(`[`, 2)
new_coldata <- new_coldata[(new_coldata$response != "0") & (new_coldata$response != "unable to determine"), ]
new_coldata[(new_coldata$response == "progressive disease") | (new_coldata$response == "stable disease"), ]$response <- "NR"
new_coldata[(new_coldata$response == "partial response") | (new_coldata$response == "complete response"), ]$response <- "R"
new_coldata$`survival_time_pfs/survival_time_os` <- NA
new_coldata$survival_unit <- NA
new_coldata$`event_occurred_pfs/event_occurred_os` <- NA
coldata <- coldata[new_coldata$patientid, ]
new_coldata$sex <- coldata$`gender:ch1`
new_coldata$age <- coldata$`age:ch1`
new_experiment <- MultiAssayExperiment(colData = new_coldata)

print('getBM')

assays <- assays[, new_coldata$patientid]
genes <- data@featureData@data
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
affyids <- genes$ID
id_version <- "affy_hg_u133a_2"

mapping <- getBM(
  attributes = c(
    id_version,
    "hgnc_symbol", "start_position", "end_position", "ensembl_gene_id_version"
  ), filters = id_version,
  values = affyids, mart = ensembl
)

print('get MultiAssay')

mapping <- mapping[mapping$ensembl_gene_id_version != "", ]
assays <- assays[as.character(mapping[[id_version]]), ]
mapping <- mapping[!duplicated(mapping$ensembl_gene_id_version), ]
assays <- assays[as.character(mapping[[id_version]]), ]
rownames(assays) <- mapping$ensembl_gene_id_version
genes <- genes[!duplicated(genes$ID), ]
rownames(genes) <- genes$ID
genes <- genes[mapping$affy_hg_u133a_2, ]
rownames(genes) <- mapping$ensembl_gene_id_version
genes$ID <- mapping$ensembl_gene_id_version
colnames(assays) <- strsplit(coldata$characteristics_ch1, ": ") %>% map_chr(`[`, 2)
counts <- SummarizedExperiment(assays = list(assays), colData = new_coldata, rowData = genes)

# assays <- assays[mapping$affy_hg_u133a_2,]
# rownames(assays) <- mapping$ensembl_gene_id_version
length <- mapping$end_position - mapping$start_position
x <- assays / length
tpm <- t(t(x) * 1e6 / colSums(x))
tpm <- as.data.frame(tpm)
tpm_counts <- SummarizedExperiment(assays = list(tpm), colData = new_coldata, rowData = genes)
experiment_list <- list(expr_gene_counts = counts, expr_gene_tpm = tpm)
experiment_list <- ExperimentList(experiment_list)
new_experiment@ExperimentList <- experiment_list

saveRDS(new_experiment, paste0(work_dir, filename))
