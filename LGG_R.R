library(TCGAbiolinks) #interface com o tcga
query_LGG <- GDCquery(project = "TCGA-LGG", 
                      data.category = "Transcriptome Profiling", 
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "STAR - Counts") #uery de pesquisa para ir buscar o exato "pacote" de dados
GDCdownload(query_LGG)
data_rna_LGG <- GDCprepare(query_LGG, summarizedExperiment = TRUE) #guardamos os nossos dados e os metadados
data_rna_LGG <- GDCprepare(query_LGG, summarizedExperiment = TRUE, save = TRUE, save.filename = "TCGA_LGG.rda")
load("TCGA_LGG.rda") #para ser mais f?cil aceder aos dados
class(data_rna_LGG)
dim(data_rna_LGG) #ver o n? de linhas e colunas
names(data_rna_LGG) #cada linha ? um gene 
colnames(data_rna_LGG) #cada coluna ? uma amostra-paciente

library(Biobase)
meta_LGG = colData(data_rna_LGG) #buscar apenas os metadadaos
dim(meta_LGG) #ver as suas dimens?es onte as linhas s?o pacientes e as colunas

cols_with_not_reported <- which(sapply(meta_LGG,function(x) sum(x == "not reported", na.rm = TRUE)) > 50)
cols_with_Not_Reported <- which(sapply(meta_LGG,function(x) sum(x == "Not Reported", na.rm = TRUE)) > 50)
cols_with_NA <- which(sapply(meta_LGG, function(x) sum(is.na(x))) > 50)
metadata_matriz_clean <- meta_LGG[, -c(cols_with_not_reported, cols_with_Not_Reported, cols_with_NA)] #retirar as colunas que nao queremos
dim(metadata_matriz_clean) #vamos trabalahr a partir desta j? filtrada

freq_table <- table(metadata_matriz_clean$primary_diagnosis)
pie(freq_table,
    main = "Distribui??o de diagn?sticos prim?rios",
    col = rainbow(length(freq_table)),
    border = "white",
    labels = paste(names(freq_table), ": ", freq_table, sep = ""))


counts <- table(metadata_matriz_clean$gender)
colors <- c("pink", "lightblue")
pie(counts, col = colors, main = "Distribui??o por g?nero")
legend("topright", legend = names(counts), fill = colors) #gr?fico circulas para mostrar a distribui??o de genero

chisq.test(counts) # para ver estatisticamente se h? diferen?as entre male e female , sem bem que isso em termos bio nao diz muito nesta fase

sem_NA=na.omit(metadata_matriz_clean$age_at_index) #pr?-filtragem
hist(sem_NA, breaks = seq(10, 90, 10), main = "Idades dos pacientes", xlab = "Faixa et?ria", ylab = "Frequ?ncia") # facilitar a visualiza??o com um histograma
summary(sem_NA)
sd(sem_NA)
shapiro.test(sem_NA) #verificar a normalidade; pvalue <0.05 portanto nao tem distribui??o normal

merged_data <- data.frame(status_vital = (metadata_matriz_clean$vital_status), 
                          idade_diagnostico = (metadata_matriz_clean$age_at_index))
merged_data #jun??o das duas coluans de metadados que queremos analisar
merged_data <- na.omit(merged_data)
index_not_reported <- which(apply(merged_data, 1, function(x) any(x == "Not Reported")))
merged_data_clean <- subset(merged_data, !row.names(merged_data) %in% index_not_reported)
boxplot(idade_diagnostico ~ status_vital, data = merged_data_clean, 
        main = "Idade de diagn?stico por status vital",
        xlab = "Status vital", ylab = "Idade de diagn?stico")
means <- tapply(merged_data_clean$idade_diagnostico, merged_data_clean$status_vital, mean)
points(means, col = "red", pch = 18, cex = 2, lwd = 2, add = TRUE)

dead_data <- subset(merged_data_clean, merged_data_clean$status_vital == "Dead") #assocair aos mortos
alive_data <- subset(merged_data_clean, merged_data_clean$status_vital == "Alive") #associar apenas os vivo
var_test <- var.test(alive_data$idade_diagnostico, dead_data$idade_diagnostico) #ver se as variancia sao true or false ; s?o false mas muito perto de 0.05
resultado_teste <- t.test(alive_data$idade_diagnostico, dead_data$idade_diagnostico,var.equal = FALSE) #pvalue <0.05 que sugere que h? uma diferen?a real entre as idades de diagn?stico dos pacientes vivos e falecidos.


merged_data_3 <- data.frame(status_vital = metadata_matriz_clean$vital_status, 
                          idade_diagnostico = metadata_matriz_clean$age_at_index,
                          gender = metadata_matriz_clean$gender)
merged_data_3_clean <- na.omit(merged_data_3)
dim(merged_data_3_clean)
index_not_reported <- which(apply(merged_data_3_clean, 1, function(x) any(x == "Not Reported")))
merged_data_clean_3 <- subset(merged_data_3_clean, !row.names(merged_data) %in% index_not_reported)
dim(merged_data_clean_3)
tabela <- table(merged_data_clean_3$status_vital, merged_data_clean_3$idade_diagnostico, merged_data_clean_3$gender)

# Criar novo dataframe apenas com informa??es relevantes
dead_data <- subset(merged_data_clean_3, merged_data_clean_3$status_vital == "Dead")
alive_data <- subset(merged_data_clean_3, merged_data_clean_3$status_vital == "Alive")

# Configurar cores das barras
cores <- c("pink", "lightblue")

# Criar tabela de conting?ncia
tabela_dead <- table(dead_data$gender, cut(dead_data$idade_diagnostico, breaks = seq(10, 90, 10)))
tabela_alive <- table(alive_data$gender, cut(alive_data$idade_diagnostico, breaks = seq(10, 90, 10)))

# Definir nomes das colunas e linhas
colnames(tabela_dead) <- paste("Dead", colnames(tabela_dead))
colnames(tabela_alive) <- paste("Alive", colnames(tabela_alive))
rownames(tabela_dead) <- c("Female", "Male")
rownames(tabela_alive) <- c("Female", "Male")

# Criar gr?fico de barras
barplot(cbind(tabela_dead, tabela_alive), beside = TRUE, col = cores,
        main = "Compara??o entre pacientes mortos e vivos",
        xlab = "Sexo", ylab = "Contagem",
        legend.text = c("Female", "Male"), args.legend = list(x = "topright"))


# Parte Rodrigo (analise dados)

#criar a matriz:

# Extrai a matriz de expressão do objeto SummarizedExperiment
expr_matrix <- assay(data_rna_LGG)

# Converte a matriz em um data frame
expr_df <- as.data.frame(expr_matrix)

# Remove quaisquer linhas com valores NA ou NaN
expr_df <- expr_df[complete.cases(expr_df), ]
nrow(expr_df)
summary(expr_df)

# Converte o data frame em uma matriz
expr_matrix <- as.matrix(expr_df)
summary(expr_matrix)

# Filtra genes cujo rácio do máximo valor sobre o mínimo valor de expressão seja superior a 2
#maximos = apply(expr_matrix,1,max)
# = apply(expr_matrix,1,min)
#min_nonzero <- min(minimos[minimos != 0])
#minimos[minimos == 0] <- min_nonzero + 1e-10
#vl = maximos/min_nonzero >2
#ALLm2=expr_matrix[vl,]
#ALLm2


library(DESeq2)

data_de <- data_rna_lgg[,!is.na(data_rna_lgg$paper_IDH.status)]

ddsSE <- DESeqDataSet(data_de, design = ~ paper_IDH.status)

ddsSE


keep <- rowSums(counts(ddsSE)) >= 50


keep

ddsSE <- ddsSE[keep,]
ddsSE <- DESeq(ddsSE)

nrow(ddsSE)

-resultsNames(ddsSE)

res <- results(ddsSE, name = "paper_IDH.status_WT_vs_Mutant")
res
dea <- as.data.frame(res)
dea

summary(res)
plot(res)


head(assay(data_rna_LGG), )


seqinfo(rowRanges(data_rna_LGG))
















# t test para filtrar:

# Aplica o teste t a cada gene
t_test <- apply(expr_matrix, 1, function(x) t.test(x)$p.value)
t_test


# Define o limiar de significância
p_threshold <- 0.05

# Seleciona os genes com valores de p abaixo do limiar
sig_genes <- names(t_test[t_test > p_threshold])
sig_genes


# Filtra a matriz de dados com os genes selecionados
data_filtered <- expr_matrix[, sig_genes]


colnames(expr_matrix)
common_genes <- intersect(colnames(data_rna_LGG), sig_genes)

common_genes

