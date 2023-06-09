#Packages usados:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("genefilter")
BiocManager::install("fgsea")
BiocManager::install("limma")
BiocManager::install("factoextra")
install.packages("xfun")
install.packages("SummarizedExperiment") 
install.packages("factoextra")
install.packages("Rtsne")


library("DESeq2")
library("TCGAbiolinks")
library("Biobase")
library("DESeq2")
library("ggbeeswarm")
library("pheatmap")
library("org.Hs.eg.db")
library("fgsea")
library("ggplot2")
library("factoextra")
library("limma")
library("genefilter")
library("SummarizedExperiment")
library("factoextra")
library("Rtsne")

#ObtenC'C#o dos dados

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
str(data_rna_LGG,give.attr=FALSE)
summary(data_rna_LGG)
meta_LGG = colData(data_rna_LGG)
dim(meta_LGG)
colnames(meta_LGG)



#Metadados e sua anC!lise descritiva

cols_with_not_reported <- which(sapply(meta_LGG,function(x) sum(x == "not reported", na.rm = TRUE)) > 50)

cols_with_Not_Reported <- which(sapply(meta_LGG,function(x) sum(x == "Not Reported", na.rm = TRUE)) > 50)

cols_with_NA <- which(sapply(meta_LGG, function(x) sum(is.na(x))) > 50)

metadata_matriz_clean <- meta_LGG[, -c(cols_with_not_reported, cols_with_Not_Reported, cols_with_NA)]

dim(metadata_matriz_clean)

metadados_teste=as.data.frame(metadata_matriz_clean)

summary(metadados_teste)

freq_table <- table(metadata_matriz_clean$primary_diagnosis)

pie(freq_table,
    main = "DistribuiC'C#o de diagnC3sticos primC!rios",
    col = rainbow(length(freq_table)),
    border = "white",
    labels = paste(names(freq_table), ": ", freq_table, sep = ""))

counts <- table(metadata_matriz_clean$gender)
colors <- c("pink", "lightblue")
pie(counts, col = colors, main = "DistribuiC'C#o por gC)nero")
legend("topright", legend = names(counts), fill = colors)

chisq.test(counts)

sem_NA=na.omit(data_rna_LGG$age_at_index)

hist(sem_NA, breaks = seq(10, 90, 10), main = "Idades dos pacientes", xlab = "Faixa-etC!ria", ylab = "FrequC*ncia")

summary(sem_NA)

shapiro.test(sem_NA)

merged_data <- data.frame(status_vital = (metadata_matriz_clean$vital_status), 
                          idade_diagnostico = (metadata_matriz_clean$age_at_index))

merged_data <- na.omit(merged_data)

index_not_reported <- which(apply(merged_data, 1, function(x) any(x == "Not Reported")))

merged_data_clean <- subset(merged_data, !row.names(merged_data) %in% index_not_reported)

boxplot(idade_diagnostico ~ status_vital, data = merged_data_clean, 
        main = "Idade de diagnC3stico por status vital",
        xlab = "Status vital", ylab = "Idade de diagnC3stico")

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
        main = "ComparaC'C#o entre a frquencia  de pacientes mortos e vivos em funC'C#o da faixa etC!ria e tendo em atenC'C#o o genero dos mesmos",
        xlab = "Faixa etC!ria", ylab = "Contagem",
        legend.text = c("Female", "Male"), args.legend = list(x = "topright"))


#Dados e sua filtragem
sum(is.na(data_rna_LGG))
teste_h <- DESeqDataSetFromMatrix(countData = assay(data_rna_LGG), colData = colData(data_rna_LGG), design = ~ 1)
teste_h <- DESeq(teste_h)
res_t <- results(teste_h)
summary(res_t)
normalized_counts <- counts(teste_h)
mean_value <- mean(normalized_counts[1, ])
normali

data_de <- data_rna_LGG[,!is.na(data_rna_LGG$paper_IDH.status)] #filtragem das amostras

ddsSE <- DESeqDataSet(data_de, design = ~ paper_IDH.status)

keep <- which(rowSums(counts(ddsSE) >= 10) >= 3) #filtragem tendo em conta uma contagem minima

ddsSE_filtrado <- ddsSE[keep, ] 
dim(ddsSE_filtrado)

ddsSE_norm <- DESeq(ddsSE_filtrado)
resultsNames(ddsSE_norm)
res <- results(ddsSE_norm, name = "paper_IDH.status_WT_vs_Mutant") 
dea <- as.data.frame(res)
summary(dea)
desvio_padrao <- apply(dea, 1, sd)
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

dim(genes_padj_fi)
#AnC!lise da expressC#o diferencial:


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

pheatmap(vsd.counts, show_colnames = F, annotation_col =anno , main="20 genes com maior diferenC'a de expressC#o\n entre os mutantes e nC#o mutantes")

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

##pathways <- gmtPathways("C:/Users/rodri/OneDrive/Documentos/h.all.v7.4.entrez.gmt")
###pathways <- gmtPathways("C:/Users/Karyna/Desktop/Github/Extracao_de_Dados/Trabalho_R/h.all.v7.4.entrez.gmt")
pathways <- gmtPathways("C:/Users/guilh/OneDrive/Documentos/GitHub/Extracao_de_Dados/Trabalho_R/h.all.v7.4.entrez.gmt")

fgseaRes <- fgsea(pathways, ranks)

dim(fgseaRes)

head(fgseaRes[order(padj), ])

ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA")

#Clustering e outros métodos não supervisionados

##Redução de dimensionalidade

##PCA

###A reduçção de dimensionalidade pode ser uma estratégia para determiar o número de genes que representam uma certa parte da variação do modelo.
dim (ddsSE_norm)
data_rna_LGG_matrix <- as.matrix(assay(ddsSE_norm))

# Transpor a matriz
data_rna_LGG_transposed <- t(data_rna_LGG_matrix)
dim(data_rna_LGG_transposed)
# Realizar o PCA
pcares <- prcomp(data_rna_LGG_transposed, scale. = TRUE)

# Verificar os desvios padrão dos componentes principais
loadings=pcares$rotation

# Resumo estatístico dos componentes principais
summary(pcares)

# Gráfico da variância explicada pelo PCA
screeplot(pcares, type = "lines", main = "Variância-pcares")
plot(1:length(pcares$sdev), pcares$sdev^2, type = "line", xlab = "Componentes Principais", ylab = "Variância", main = "Variância-pcares")

# Adicionar legendas nos eixos
xlab("Componentes Principais")
ylab("Variância Explicada (%)")

min(which(summary(pcares)$importance[3,]>0.90))
var_total <- sum(pcares$sdev^2)
soma_variancias <- sum(head(pcares$sdev^2, 100))

cumulative_variance <- cumsum(variance)
 summary(pcares$rotation[,1:2])

###Foi realizada também uma análise de componentes principais (PCA) sobre estes dados de forma a visualizar os dados e efetuar uma redução de dimensionalidade. Os dados já se encontravam normalizados, e do PCA temos que a primeira dimensão agrega 20.7% da variabilidade da amostra, a segunda dimensão 6.9% e a terceira dimensão 3.6%, perfazendo um total de cerca de 31,2% de variabilidade cumulativa. Para perfazer mais de 90% da variabilidade total do dataset, seria necessário acumular 3 componentes.
plot(pcares$x, col = ddsSE_norm$paper_IDH.status, pch = 19)
legend("topright",legend=levels(ddsSE_norm$paper_IDH.status), col = 1:2, pch=19)
plot(pcares)
plot(pcares2)

library(scatterplot3d)
factor <- as.factor(ddsSE_norm$paper_IDH.status)

scatterplot3d(pcares$x[, 1], pcares$x[, 2], pcares$x[, 3], color = levels(factor),
              pch = 19, xlab = "PC1", ylab = "PC2", zlab = "PC3")
install.packages("factoextra")
library("factoextra")

fviz_famd_ind(pcares, geom = c("point"), col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              palette = "rainbow", addEllipses = FALSE, ellipse.type = "confidence",
              ggtheme = theme_minimal(), repel = TRUE, labels = F)  
#tsne


data_rna_LGG_matrix <- as.matrix(assay(ddsSE_norm))
data_rna_LGG_transposed <- t(data_rna_LGG_matrix)
Rtsne(data_rna_LGG_transposed)
data_rna_LGG_matrix_nd = data_rna_LGG_transposed[!duplicated(data_rna_LGG_transposed),]
dim(data_rna_LGG_matrix_nd)
res_tnse = Rtsne(data_rna_LGG_matrix_nd)
plot(res_tnse$Y, col = ddsSE_norm$paper_IDH.status, pch = 19)
install.packages("rgl")

#mds
data_rna_LGG_matrix <- as.matrix(assay(ddsSE_norm))

# Cálculo da matriz de distâncias euclidianas
dist_matrix <- dist(scale(data_rna_LGG_matrix))

# Aplicação do MDS
mds_result <- cmdscale(dist_matrix, k = 3, eig = TRUE)

# Cálculo da variância explicada por cada coordenada
var_pc1 <- sum(mds_result$eig[1]) / sum(abs(mds_result$eig)) * 100
var_pc2 <- sum(mds_result$eig[2]) / sum(abs(mds_result$eig)) * 100
var_pc3 <- sum(mds_result$eig[3]) / sum(abs(mds_result$eig)) * 100

# Legendas para os eixos
lab_pc1 <- paste("Coordinate 1 (var.", sprintf("%.3f", var_pc1), "%)")
lab_pc2 <- paste("Coordinate 2 (var.", sprintf("%.3f", var_pc2), "%)")

# Plotar o resultado do MDS
plot(mds_result$points[, c(1, 2)], xlab = lab_pc1, ylab = lab_pc2, main = "MDS Plot")

### Ao observar as linhagens representadas graficamente ao longo dos primeiro e segundo componentes, temos que ao colorir as linhagens pela sua qualidade de representação “cos2” que as linhagens mais próximas do 0 são aquelas cuja variação se encontra menos explicada pelos dois componenetes representados, enquanto que aquelas mais distantes ao longo do primeiro e segundo eixo são aquelas que se encontram melhor diferenciadas.

##Clustering Hierárquico
#por paciente##################################################
tt_mdr = rowttests(t(data_rna_LGG_matrix))
rank_de_mdr = order(tt_mdr$p.value)
genes_de_mdr = rank_de_mdr[1:30]
data_rna_LGG_rank = data_rna_LGG_transposed[genes_de_mdr,]

eucD = dist(data_rna_LGG_rank)
eucD
###complete
cl.hier <- hclust(eucD)
plot(cl.hier,xlab="", ylab="Distância", main="Dendograma da expressão dos 30 pacientes com menor p-value \nmétodo:complete, distância Euclidiana")

###single
cl.hier2 <- hclust(eucD, method="single")
plot(cl.hier2,xlab="", ylab="Distância", main="Dendograma da expressão dos 30 pacientes com menor p-value \nmétodo:single, distância Euclidiana")

###average
cl.hier3 <- hclust(eucD, method="average")
plot(cl.hier3,xlab="", ylab="Distância", main="Dendograma da expressão dos 30 pacientes com menor p-value \nmétodo:average, distância Euclidiana")

heatmap(data_rna_LGG_rank, labCol = F, main="Expressão dos 30 pacientes com menor p-value")

#por genes##################################################
tt_mdr_g = rowttests(data_rna_LGG_matrix)
rank_de_mdr_g = order(tt_mdr_g$p.value)
genes_de_mdr_g = rank_de_mdr_g[1:30]
data_rna_LGG_rank = data_rna_LGG_matrix[genes_de_mdr,]

eucD = dist(data_rna_LGG_rank)

###complete
cl.hier <- hclust(eucD)
plot(cl.hier,xlab="", ylab="Distância", main="Dendograma da expressão dos 30 genes com menor p-value \nmétodo:complete, distância Euclidiana")

###single
cl.hier2 <- hclust(eucD, method="single")
plot(cl.hier2,xlab="", ylab="Distância", main="Dendograma da expressão dos 30 genes com menor p-value \nmétodo:single, distância Euclidiana")

###average
cl.hier3 <- hclust(eucD, method="average")
plot(cl.hier3,xlab="", ylab="Distância", main="Dendograma da expressão dos 30 genes com menor p-value \nmétodo:average, distância Euclidiana")

heatmap(data_rna_LGG_rank, labCol = F, main="Expressão dos 30 genes com menor p-value")

### Escolher um metadado para fazer comparações?? (não tenho a certeza disto)

##k-means clustering

### optimal number of clusters 


fviz_nbclust(t(data_rna_LGG_matrix), kmeans, method = "silhouette")

##Para efetuar o clustering por k-means, de forma a efetuar uma classificação dos grupos observados no PCA, foi primeiro realizada uma siluette analysis sobre os dados logaritmizados com recurso a função fviz_nbclust, que nos indicou que a solução ótima residia em 2 clusters.

### Map of predicted clusters
install.packages("heatmaply")
library(heatmaply)

resKmeans <- kmeans(t(data_rna_LGG_matrix),centers=6)
resKmeans
centroides=resKmeans$cluster
table_result=table(centroides, ddsSE_filtrado$paper_IDH.status)
#plot(data_rna_LGG_transposed[,1],ddsSE_norm$paper_IDH.status, col = resKmeans$cluster, pch = 16)
library(ggplot2)
library(cluster)
#cotovelo
ofs <- c()

for (k in 2:10) {
  kmeans <- kmeans(t(data_rna_LGG_matrix), centers = k, nstart = 10)
  ofs <- c(ofs, kmeans$tot.withinss)
}

plot_data <- data.frame(num_clusters = 2:10, wss = ofs)

ggplot(plot_data, aes(x = num_clusters, y = wss)) +
  geom_line() +
  geom_point() +
  labs(x = "Num Clusters", y = "WSS") +
  theme_minimal()




# Adicionar legendas aos pontos
legend("topright", legend = 1:3, col = 1:3, pch = 8, cex = 1.2)

# Adicionar título e rótulos dos eixos
title("Gráfico de Clusters - K-means")
xlabel <- "Variável X"
ylabel <- "Variável Y"
xlab(xlabel)
ylab(ylabel)

#não consegui fazer o gráfico dos clusters k-mean porque está a faltar um metadado

#Análise supervisionada (Machine Learning)
#parece tudo bem- temos divisao data treino e teste
#classe 1 e 2 
#temos gráfico para mostrar a distribuição destas classe
set.seed(16718)
data_rna_LGG_data <- data.frame(data_rna_LGG_matrix)
ml_mutants <- as.data.frame(cbind(group = meta_LGG$paper_IDH.status, t(data_rna_LGG_data)))
ml_mutants_na <- na.omit(ml_mutants)
dim(ml_mutants_na)
ml_mutants_na$group = as.factor(ml_mutants_na$group)
ml_mutants_na$group
select_genee <- rownames(head(resOrdered,16000)) #top 16000 pelo p-adjust

# Crie um novo dataframe apenas com os genes selecionados
ml_mutants_selected <- ml_mutants_na[, c("group", select_genee)]

ml_mutants_selected$group
frequencia <- table(ml_mutants_na$group)
cores <- rainbow(length(frequencia))
pie(frequencia, col = cores)
# 1= MUTANT   2= WT
#Percentagem de exemplos corretamente classificados
pecc = function(obs,pred) sum(obs==pred)/length(obs)

# Raízquadradadamédiadoquadradodoserro
rmse = function(obs, pred) sqrt(mean((obs-pred)^2)) 

# Médiadosdesviosabsolutos
mad = function(obs, pred) mean(abs(obs-pred))

#Divisão em train e test em 70%, 30%
ind = sample(2, nrow(ml_mutants_selected), replace=TRUE, prob=c(0.7, 0.3)) 
trainData = ml_mutants_selected[ind==1,]
testData = ml_mutants_selected[ind==2,]
dim(trainData)
dim(testData)
table(trainData$group)
table(testData$group)
##############################################################
#Todavia isso nao garante a preservação das proporções das classe no treino e test
install.packages("rsample")
library(rsample)
library(caret)

set.seed(16718)
# Criar as dobras estratificadas usando o pacote rsample
folds <- rsample::vfold_cv(ml_mutants_selected, strata = "group", v = 10)

# Configurar o objeto trainControl com as dobras estratificadas
cv.control <- trainControl(method = "repeatedcv", number = 10, repeats = 5, index = folds$split)

# Executar o treinamento usando o trainControl modificado
set.seed(16718)
group_knn_cv <- train(group ~ ., data = trainData[,1:16001], method = "knn", tuneGrid = expand.grid(k = 1:10), trControl = cv.control)
best_k <- group_knn_cv$bestTune$k
best_k

pred_kn_cv <- predict(group_knn_cv, newdata = testData[,1:16001])
pred_kn_cv

#métricas
predictions <- pred_kn_cv


# Criar a matriz de confusão
confusion_matrix <- confusionMatrix(predictions, testData$group)

confusion_matrix
precision <- confusion_matrix$byClass["Pos Pred Value"]
recall <- confusion_matrix$byClass["Sensitivity"]
accuracy <- confusion_matrix$overall["Accuracy"]
f1_score <- confusion_matrix$byClass["F1"]
sensitivity <- confusion_matrix$byClass["Sensitivity"]

# Imprimir as métricas
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("Accuracy:", accuracy, "\n")
cat("F1-score:", f1_score, "\n")
cat("Sensitivity:", sensitivity, "\n")
######################################################
#naive bayes
group_nb_cv <- train(group ~ ., data = trainData[,1:16001], method = "nb", trControl = cv.control)
pred_nb_cv <- predict(group_nb_cv, newdata = testData[,1:16001])

#métricas
predictions <- pred_nb_cv


# Criar a matriz de confusão
confusion_matrix <- confusionMatrix(predictions, testData$group)

confusion_matrix
precision <- confusion_matrix$byClass["Pos Pred Value"]
recall <- confusion_matrix$byClass["Sensitivity"]
accuracy <- confusion_matrix$overall["Accuracy"]
f1_score <- confusion_matrix$byClass["F1"]
sensitivity <- confusion_matrix$byClass["Sensitivity"]

# Imprimir as métricas
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("Accuracy:", accuracy, "\n")
cat("F1-score:", f1_score, "\n")
cat("Sensitivity:", sensitivity, "\n")
#####################################################
# Decision Trees
set.seed(16718)

# Definir o grid de hiperparâmetros
grid <- expand.grid(cp = seq(0.01, 0.1, by = 0.01))

# Treinar o modelo usando cross-validation
group_tree_cv <- train(group ~ ., data = trainData[,1:16001], method = "rpart", tuneGrid = grid, trControl = cv.control)

# Obter o melhor hiperparâmetro
best_cp <- group_tree_cv$bestTune$cp

# Treinar o modelo final com o melhor hiperparâmetro
final_model <- rpart(group ~ ., data = trainData[,1:16001], method = "class", cp = best_cp)

# Fazer previsões em novos dados
pred_tree_cv <- predict(final_model, newdata = testData[,1:16001], type = "class")

# Métricas
confusion_matrix <- confusionMatrix(pred_tree_cv, testData$group)
precision <- confusion_matrix$byClass["Pos Pred Value"]
recall <- confusion_matrix$byClass["Sensitivity"]
accuracy <- confusion_matrix$overall["Accuracy"]
f1_score <- confusion_matrix$byClass["F1"]
sensitivity <- confusion_matrix$byClass["Sensitivity"]

# Imprimir as métricas
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("Accuracy:", accuracy, "\n")
cat("F1-score:", f1_score, "\n")
cat("Sensitivity:", sensitivity, "\n")

#k-Nearest Neighbors -Rodrigo

library(class)
#esta mal deveria ser em relação ao nosso dataset!
knn_pred = knn(trainData[,1:4], testData[,1:4], trainData$group) 
knn_pred

# 1= MUTANT   2= WT
t = table(knn_pred, testData$group)
t

pecc(knn_pred, testData$group) #0.7898089

#Naive Bayes 

library(e1071)

model = naiveBayes(group ~ ., trainData[,1:4]) #não sei para que serve o [,1:4], no entanto sem isso, não corre porque são muitos dados
nb_pred = predict(model, testData)
nb_pred

table(nb_pred, testData$group)

pecc(testData$group, nb_pred) #0.8535032

#Decision Trees

library(party)
ml_mutants_na_ctree <- ctree(group ~., data=trainData[,1:4])
print(ml_mutants_na_ctree)

plot(ml_mutants_na_ctree)

testPred <- predict(ml_mutants_na_ctree, testData)
testPred

table(testPred, testData$group)

pecc(testData$group, testPred) #0.8216561

#Regression Trees

library(rpart)
arvreg = rpart(group ~., data = trainData[,1:4])
plot(arvreg)
text(arvreg)
val_prev = predict(arvreg, testData)
val_prev

#não consigo fazer o que está em baixo porque a nossa variavel é qualitiativa e não quantitativa (eu acho)
rmse(val_prev, testData$group)
mad(val_prev, testData$group)

### Modelos funcionais lineares

#Regressão usando PLS

#não consigo correr isto porque estou a ter problemas com o package, mas penso que esteja bem
library(caret) 
pls.ml_mutants_na= plsda(trainData[,1:4], trainData[,5]) 
pred.pls.ml_mutants_na = predict(pls.ml_mutants_na, testData[1:4]) 
pecc(pred.pls.ml_mutants_na, testData$group)

table(pred.pls.ml_mutants_na, testData$group)

#Análise discriminante

library(MASS)
lda.model = lda(group ~., trainData[,1:4])
lda.model

test.lda = predict(lda.model, testData) 
test.lda$class

pecc(test.lda$class, testData$group) #0.8089172

#Regressão logística

x <- ml_mutants_na[sample(1:nrow(ml_mutants_na)),]
x$mutant <- x$group == "1" 
x$group <- NULL
model <- glm(mutant~ ., family = binomial(logit), data= x/10)
model

#Modelos não lineares

x$group

x$mutant

