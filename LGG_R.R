library(TCGAbiolinks) #interface com o tcga
query_LGG <- GDCquery(project = "TCGA-LGG", 
                      data.category = "Transcriptome Profiling", 
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "STAR - Counts") #uery de pesquisa para ir buscar o exato "pacote" de dados
GDCdownload(query_LGG)
data_rna_LGG <- GDCprepare(query_LGG, summarizedExperiment = TRUE) #guardamos os nossos dados e os metadados
data_rna_LGG <- GDCprepare(query_LGG, summarizedExperiment = TRUE, save = TRUE, save.filename = "TCGA_LGG.rda")
load("TCGA_LGG.rda") #para ser mais fácil aceder aos dados
class(data_rna_LGG)
dim(data_rna_LGG) #ver o nº de linhas e colunas
names(data_rna_LGG) #cada linha é um gene 
colnames(data_rna_LGG) #cada coluna é uma amostra-paciente

library(Biobase)
meta_LGG = colData(data_rna_LGG) #buscar apenas os metadadaos
dim(meta_LGG) #ver as suas dimensões onte as linhas são pacientes e as colunas

cols_with_not_reported <- which(sapply(meta_LGG,function(x) sum(x == "not reported", na.rm = TRUE)) > 50)
cols_with_Not_Reported <- which(sapply(meta_LGG,function(x) sum(x == "Not Reported", na.rm = TRUE)) > 50)
cols_with_NA <- which(sapply(meta_LGG, function(x) sum(is.na(x))) > 50)
metadata_matriz_clean <- meta_LGG[, -c(cols_with_not_reported, cols_with_Not_Reported, cols_with_NA)] #retirar as colunas que nao queremos
dim(metadata_matriz_clean) #vamos trabalahr a partir desta já filtrada

freq_table <- table(metadata_matriz_clean$primary_diagnosis)
pie(freq_table,
    main = "Distribuição de diagnósticos primários",
    col = rainbow(length(freq_table)),
    border = "white",
    labels = paste(names(freq_table), ": ", freq_table, sep = ""))


counts <- table(metadata_matriz_clean$gender)
colors <- c("pink", "lightblue")
pie(counts, col = colors, main = "Distribuição por género")
legend("topright", legend = names(counts), fill = colors) #gráfico circulas para mostrar a distribuição de genero

chisq.test(counts) # para ver estatisticamente se há diferenças entre male e female , sem bem que isso em termos bio nao diz muito nesta fase

sem_NA=na.omit(metadata_matriz_clean$age_at_index) #pré-filtragem
hist(sem_NA, breaks = seq(10, 90, 10), main = "Idades dos pacientes", xlab = "Faixa etária", ylab = "Frequência") # facilitar a visualização com um histograma
summary(sem_NA)
sd(sem_NA)
shapiro.test(sem_NA) #verificar a normalidade; pvalue <0.05 portanto nao tem distribuição normal

merged_data <- data.frame(status_vital = (metadata_matriz_clean$vital_status), 
                          idade_diagnostico = (metadata_matriz_clean$age_at_index))
merged_data #junção das duas coluans de metadados que queremos analisar
merged_data <- na.omit(merged_data)
index_not_reported <- which(apply(merged_data, 1, function(x) any(x == "Not Reported")))
merged_data_clean <- subset(merged_data, !row.names(merged_data) %in% index_not_reported)
boxplot(idade_diagnostico ~ status_vital, data = merged_data_clean, 
        main = "Idade de diagnóstico por status vital",
        xlab = "Status vital", ylab = "Idade de diagnóstico")
means <- tapply(merged_data_clean$idade_diagnostico, merged_data_clean$status_vital, mean)
points(means, col = "red", pch = 18, cex = 2, lwd = 2, add = TRUE)

dead_data <- subset(merged_data_clean, merged_data_clean$status_vital == "Dead") #assocair aos mortos
alive_data <- subset(merged_data_clean, merged_data_clean$status_vital == "Alive") #associar apenas os vivo
var_test <- var.test(alive_data$idade_diagnostico, dead_data$idade_diagnostico) #ver se as variancia sao true or false ; sáo false mas muito perto de 0.05
resultado_teste <- t.test(alive_data$idade_diagnostico, dead_data$idade_diagnostico,var.equal = FALSE) #pvalue <0.05 que sugere que há uma diferença real entre as idades de diagnóstico dos pacientes vivos e falecidos.


merged_data_3 <- data.frame(status_vital = metadata_matriz_clean$vital_status, 
                          idade_diagnostico = metadata_matriz_clean$age_at_index,
                          gender = metadata_matriz_clean$gender)
merged_data_3_clean <- na.omit(merged_data_3)
dim(merged_data_3_clean)
index_not_reported <- which(apply(merged_data_3_clean, 1, function(x) any(x == "Not Reported")))
merged_data_clean_3 <- subset(merged_data_3_clean, !row.names(merged_data) %in% index_not_reported)
dim(merged_data_clean_3)
tabela <- table(merged_data_clean_3$status_vital, merged_data_clean_3$idade_diagnostico, merged_data_clean_3$gender)

# Criar novo dataframe apenas com informações relevantes
dead_data <- subset(merged_data_clean_3, merged_data_clean_3$status_vital == "Dead")
alive_data <- subset(merged_data_clean_3, merged_data_clean_3$status_vital == "Alive")

# Configurar cores das barras
cores <- c("pink", "lightblue")

# Criar tabela de contingência
tabela_dead <- table(dead_data$gender, cut(dead_data$idade_diagnostico, breaks = seq(10, 90, 10)))
tabela_alive <- table(alive_data$gender, cut(alive_data$idade_diagnostico, breaks = seq(10, 90, 10)))

# Definir nomes das colunas e linhas
colnames(tabela_dead) <- paste("Dead", colnames(tabela_dead))
colnames(tabela_alive) <- paste("Alive", colnames(tabela_alive))
rownames(tabela_dead) <- c("Female", "Male")
rownames(tabela_alive) <- c("Female", "Male")

# Criar gráfico de barras
barplot(cbind(tabela_dead, tabela_alive), beside = TRUE, col = cores,
        main = "Comparação entre pacientes mortos e vivos",
        xlab = "Sexo", ylab = "Contagem",
        legend.text = c("Female", "Male"), args.legend = list(x = "topright"))

library(DESeq2)
data_de <- data_rna_LGG[,!is.na(data_rna_LGG$paper_IDH.status)] #as amostras que possuem uma condição definida, excluindo aquelas com valor "NA" (que representa dados ausentes).
ddsSE <- DESeqDataSet(data_de, design = ~ paper_IDH.status) #estamos a comparar a expressão genética com os grupos existentes em design,neste caso mutante e wildtype
rowSums(counts(ddsSE))>=10 #A segunda linha cria um objeto DESeqDataSet, que contém as contagens de expressão gênica e o design experimental. A variável a ser testada é especificada. Neste caso, estamos comparando os grupos Mutant e Wildtype:
# Extrair as colunas onde a soma das linhas é <= 10
keep <- which(rowSums(counts(ddsSE) >= 10) >= 3) #filtragem aprofundada 
ddsSE_filtrado <- ddsSE[keep, ] #apenas colocar os genes que têm umas counts acima de 10  e tem de ter pelo menos 10,10,10 nao pode ter tudo zero e um 11

ddsSE_norm <- DESeq(ddsSE_filtrado) # Essa função realiza a normalização dos dados, estima os parâmetros do modelo e realiza o teste de hipótese para cada gene.
resultsNames(ddsSE_norm) # extrai o nome da coluna que contém as diferenças entre as condições Mutant e Wildtype:
res <- results(ddsSE_norm, name = "paper_IDH.status_WT_vs_Mutant") #vai de encontro ao grupo de metadados escolhido em cima com design
dea <- as.data.frame(res) # converte o objeto "res" em um data frame para facilitar a visualização dos resultados

# Seleciona as expressões dos genes diferencialmente expressos
de_genes <- rownames(res)[which(dea$padj < 0.05)] #probabilidade inferior a 5 % de encontrar um falso positivo
plotMA(res, main="DESeq2", ylim=c(-2,2))

