BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))


db <- read.csv("GSE144552_gene_count_matrix.csv", header = T)

db$gene_id <- gsub("\\..", "", db$gene_id)

#replacing Ensembl gene ID with conventional gene IDs

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- db$gene_id

#here, constructed the mart functions with the geneset data from ensembl.
#then, created the genes object and passed the ensembl gene IDs from the


# Jurkat database


list <- getBM(filters = "ensembl_gene_id", 
               attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
               values = genes, mart = mart)

db<-merge(db, list, by.x = "gene_id", by.y = "ensembl_gene_id")

#move gene names next to the Ensembl gene IDs
db <- db %>% relocate(external_gene_name, .after = gene_id)
db <- subset(db, select = -entrezgene_id)

#the dataframe is now annotated and genes are ready to be analyzed using
# differential gene analysis using DEseq

BiocManager::install("DESeq2")
library(DESeq2)

#generate sample information table called coldata

coldata <- matrix(c("E0A_1_S1", "E0B_2_S2", "E0C_3_S3", "E1A_4_S4", "E1B_5_S5", "E1C_6_S6",
                    "E8A_7_S7",	"E8B_8_S8",	"E8C_9_S9","E16A_10_S10",	"E16B_11_S11",
                    "E16C_12_S12", "S0A_1_S1",	"S0B_2_S2",	"S0C_3_S3","S2A_4_S4", "S2B_5_S5",	"S2C_6_S6",
                    "S4A_7_S7",	"S4B_8_S8",	"S4C_9_S9", "S8A_10_S10",	"S8B_11_S11",	"S8C_12_S12",
                    rep("untreated_EPH334", 3), rep("1h_EPH334",3),rep("8h_EPH334",3),
                    rep("16h_EPH334",3), rep("untreated_SAHA",3), rep("2h_SAHA", 3), 
                    rep("4h_SAHA",3), rep("8h_SAHA",3)), ncol = 2)
colnames(coldata) <- c("batch", "condition")

coldata <- as.table(coldata)

#convert the input gene dataframe into a matrix

geneID <- db$gene_id
sampleIndex <- grepl("")