#Packages usados:

library("TCGAbiolinks")
library("Biobase")
library("DESeq2")
library("ggbeeswarm")
library("genefilter")
library("pheatmap")
library("org.Hs.eg.db")
library("fgsea")
library("ggplot2")

#Obtenção dos dados

query_LGG <- GDCquery(project = "TCGA-LGG", 
                      data.category = "Transcriptome Profiling", 
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "STAR - Counts")

#GDCdownload(query_LGG)

data_rna_LGG <- GDCprepare(query_LGG, summarizedExperiment = TRUE, save = TRUE, save.filename = "TCGA_LGG.rda")
load("TCGA_LGG.rda")

class(data_rna_LGG)
dim(data_rna_LGG)
names(data_rna_LGG)
colnames(data_rna_LGG)


meta_LGG = colData(data_rna_LGG)
dim(meta_LGG)


#Metadados e sua análise descritiva

cols_with_not_reported <- which(sapply(meta_LGG,function(x) sum(x == "not reported", na.rm = TRUE)) > 50)

cols_with_Not_Reported <- which(sapply(meta_LGG,function(x) sum(x == "Not Reported", na.rm = TRUE)) > 50)

cols_with_NA <- which(sapply(meta_LGG, function(x) sum(is.na(x))) > 50)

metadata_matriz_clean <- meta_LGG[, -c(cols_with_not_reported, cols_with_Not_Reported, cols_with_NA)]

dim(metadata_matriz_clean)

metadados_teste=as.data.frame(metadata_matriz_clean)

summary(metadados_teste)

freq_table <- table(metadata_matriz_clean$primary_diagnosis)

pie(freq_table,
    main = "Distribuição de diagnósticos primários",
    col = rainbow(length(freq_table)),
    border = "white",
    labels = paste(names(freq_table), ": ", freq_table, sep = ""))

counts <- table(metadata_matriz_clean$gender)
colors <- c("pink", "lightblue")
pie(counts, col = colors, main = "Distribuição por género")
legend("topright", legend = names(counts), fill = colors)

chisq.test(counts)

sem_NA=na.omit(data_rna_LGG$age_at_index)

hist(sem_NA, breaks = seq(10, 90, 10), main = "Idades dos pacientes", xlab = "Faixa-etária", ylab = "Frequência")

summary(sem_NA)

shapiro.test(sem_NA)

merged_data <- data.frame(status_vital = (metadata_matriz_clean$vital_status), 
                          idade_diagnostico = (metadata_matriz_clean$age_at_index))

merged_data <- na.omit(merged_data)

index_not_reported <- which(apply(merged_data, 1, function(x) any(x == "Not Reported")))

merged_data_clean <- subset(merged_data, !row.names(merged_data) %in% index_not_reported)

boxplot(idade_diagnostico ~ status_vital, data = merged_data_clean, 
        main = "Idade de diagnóstico por status vital",
        xlab = "Status vital", ylab = "Idade de diagnóstico")

means <- tapply(merged_data_clean$idade_diagnostico, merged_data_clean$status_vital, mean)

points(means, col = "red", pch = 18, cex = 2, lwd = 2, add = TRUE)

dead_data <- subset(merged_data_clean, merged_data_clean$status_vital == "Dead") 

alive_data <- subset(merged_data_clean, merged_data_clean$status_vital == "Alive")

var_test <- var.test(alive_data$idade_diagnostico, dead_data$idade_diagnostico)  

resultado_teste <- t.test(alive_data$idade_diagnostico, dead_data$idade_diagnostico,var.equal = FALSE)

merged_data_3 <- data.frame(status_vital = metadata_matriz_clean$vital_status, 
                            idade_diagnostico = metadata_matriz_clean$age_at_index,
                            gender = metadata_matriz_clean$gender)

merged_data_3_clean <- na.omit(merged_data_3)
dim(merged_data_3_clean)

index_not_reported <- which(apply(merged_data_3_clean, 1, function(x) any(x == "Not Reported")))

merged_data_clean_3 <- subset(merged_data_3_clean, !row.names(merged_data) %in% index_not_reported)

dim(merged_data_clean_3)

tabela <- table(merged_data_clean_3$status_vital, merged_data_clean_3$idade_diagnostico, merged_data_clean_3$gender)

dead_data <- subset(merged_data_clean_3, merged_data_clean_3$status_vital == "Dead")

alive_data <- subset(merged_data_clean_3, merged_data_clean_3$status_vital == "Alive")

cores <- c("pink", "lightblue")

tabela_dead <- table(dead_data$gender, cut(dead_data$idade_diagnostico, breaks = seq(10, 90, 10)))

tabela_alive <- table(alive_data$gender, cut(alive_data$idade_diagnostico, breaks = seq(10, 90, 10)))

colnames(tabela_dead) <- paste("Dead", colnames(tabela_dead))

colnames(tabela_alive) <- paste("Alive", colnames(tabela_alive))

rownames(tabela_dead) <- c("Female", "Male")

rownames(tabela_alive) <- c("Female", "Male")

barplot(cbind(tabela_dead, tabela_alive), beside = TRUE, col = cores,
        main = "Comparação entre a frquencia  de pacientes mortos e vivos em função da faixa etária e tendo em atenção o genero dos mesmos",
        xlab = "Faixa etária", ylab = "Contagem",
        legend.text = c("Female", "Male"), args.legend = list(x = "topright"))


#Dados e sua filtragem

data_de <- data_rna_LGG[,!is.na(data_rna_LGG$paper_IDH.status)] 

ddsSE <- DESeqDataSet(data_de, design = ~ paper_IDH.status)

keep <- which(rowSums(counts(ddsSE) >= 10) >= 3) 

ddsSE_filtrado <- ddsSE[keep, ] 

ddsSE_norm <- DESeq(ddsSE_filtrado) 

resultsNames(ddsSE_norm)

res <- results(ddsSE_norm, name = "paper_IDH.status_WT_vs_Mutant") 

dea <- as.data.frame(res)

mcols(res, use.names = TRUE)

summary(res)

plotMA(res, main="DESeq2", ylim=c(-10,10))

hist(dea$pvalue, breaks=20,col = "grey", border = "white", xlab = "P-value",
     ylab = "Number of genes", main = "P-value value distribution")

genes_pvalue_fi <- rownames(res)[which(dea$pvalue < 0.05)]

pvalue_fi= sum(dea$pvalue < 0.05, na.rm=TRUE)

hist(dea$padj, breaks=20,col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

genes_padj_fi <- rownames(res)[which(dea$padj < 0.01)]

padj_fi= sum(dea$padj < 0.01, na.rm=TRUE)

#Análise da expressão diferencial:


topGene <- rownames(res)[which.min(res$padj)]

plotCounts(ddsSE_norm, gene = topGene, intgroup=c("paper_IDH.status"))

geneCounts <- plotCounts(ddsSE_norm, gene = topGene, intgroup = c("paper_IDH.status","vital_status"), 
                         returnData = TRUE)
ggplot(geneCounts, aes(x = paper_IDH.status, y = count, color = vital_status)) + 
  scale_y_log10() +  geom_beeswarm(cex = 3)

vsd <- varianceStabilizingTransformation(ddsSE_norm, blind = FALSE)

resOrdered <- res[order(res$padj),]

select <- rownames(head(resOrdered,20))

vsd.counts <- assay(vsd)[select,]

df <- as.data.frame(colData(ddsSE_norm)[,c("paper_IDH.status")])

anno <- as.data.frame(colData(vsd)[, c("paper_IDH.status", "vital_status")])

pheatmap(vsd.counts, show_colnames = F, annotation_col =anno , main="20 genes com maior diferença de expressão\n entre os mutantes e não mutantes")

##Enriquecimento 

get_entrez <- function(x){
  unlist(strsplit(x, split="[.]+"))[2]
}

ann <- select(org.Hs.eg.db,keys=sapply(rownames(res), get_entrez),columns=c("ENTREZID","SYMBOL","GENENAME"))
head(ann)

all_results.annotated <- cbind(res, ann)
head(all_results.annotated)

results.ord <- all_results.annotated[ order(-all_results.annotated[,"log2FoldChange"]), ]

ranks <- results.ord$log2FoldChange

names(ranks) <- results.ord$ENTREZID

pathways <- gmtPathways("C:/Users/rodri/OneDrive/Documentos/h.all.v7.4.entrez.gmt")

fgseaRes <- fgsea(pathways, ranks)

dim(fgseaRes)

head(fgseaRes[order(padj), ])

ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA")
