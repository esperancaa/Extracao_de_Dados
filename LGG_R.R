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
dim(sem_NA)
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

library(DESeq2)
data_de <- data_rna_LGG[,!is.na(data_rna_LGG$paper_IDH.status)] #as amostras que possuem uma condi??o definida, excluindo aquelas com valor "NA" (que representa dados ausentes).
ddsSE <- DESeqDataSet(data_de, design = ~ paper_IDH.status) #estamos a comparar a express?o gen?tica com os grupos existentes em design,neste caso mutante e wildtype
rowSums(counts(ddsSE))>=10 #A segunda linha cria um objeto DESeqDataSet, que cont?m as contagens de express?o g?nica e o design experimental. A vari?vel a ser testada ? especificada. Neste caso, estamos comparando os grupos Mutant e Wildtype:
# Extrair as colunas onde a soma das linhas ? <= 10
keep <- which(rowSums(counts(ddsSE) >= 10) >= 3) #filtragem aprofundada 
ddsSE_filtrado <- ddsSE[keep, ] #apenas colocar os genes que t?m umas counts acima de 10  e tem de ter pelo menos 10,10,10 nao pode ter tudo zero e um 11

ddsSE_norm <- DESeq(ddsSE_filtrado) # Essa fun??o realiza a normaliza??o dos dados, estima os par?metros do modelo e realiza o teste de hip?tese para cada gene.
resultsNames(ddsSE_norm) # extrai o nome da coluna que cont?m as diferen?as entre as condi??es Mutant e Wildtype:
res <- results(ddsSE_norm, name = "paper_IDH.status_WT_vs_Mutant") #vai de encontro ao grupo de metadados escolhido em cima com design
dea <- as.data.frame(res) # converte o objeto "res" em um data frame para facilitar a visualiza??o dos resultados

# Seleciona as express?es dos genes diferencialmente expressos
de_genes <- rownames(res)[which(dea$padj < 0.05)] #probabilidade inferior a 5 % de encontrar um falso positivo
plotMA(res, main="DESeq2", ylim=c(-2,2))


#PARTE NOVAAAAAAAAAA:
res
nrow(res)
mcols(res, use.names = TRUE)
#Adaptar este texto à nossa problemática:

#A primeira coluna, baseMean, é apenas a média dos valores de contagem normalizados, divididos pelos factores de tamanho, recolhidos em todas as amostras do conjunto de dados DESeq. 
#As restantes quatro colunas referem-se a um contraste específico, nomeadamente a comparação do nível trt sobre o nível untrt para a variável factor dex. A seguir descobriremos como obter 
#outros contrastes. A coluna log2FoldChange é a estimativa do tamanho do efeito. Ela diz-nos quanto a expressão do gene parece ter mudado devido ao tratamento com dexametasona em comparação 
#com amostras não tratadas. Este valor é reportado numa escala logarítmica à base 2: por exemplo, uma alteração de 1,5 da dobra log2 significa que a expressão do gene é aumentada por um 
#factor multiplicativo de 21.5≈2.82.Evidentemente, esta estimativa tem uma incerteza associada, que está disponível na coluna lfcSE, a estimativa de erro padrão para a estimativa da 
#mudança de log2 fold. Podemos também expressar a incerteza de uma estimativa de tamanho de efeito particular como resultado de um teste estatístico. O objectivo de um teste de 
#expressão diferencial é testar se os dados fornecem provas suficientes para concluir que este valor é realmente diferente de zero. DESeq2 executa para cada gene um teste de hipótese 
#ver se a evidência é suficiente para decidir contra a hipótese nula de que existe efeito zero do tratamento no gene e que a diferença observada entre tratamento e controlo foi meramente 
#causada pela variabilidade experimental (ou seja, o tipo de variabilidade que se pode esperar entre diferentes amostras no mesmo grupo de tratamento). Como é habitual nas estatísticas, 
#o resultado deste teste é reportado como um valor p, e é encontrado na coluna pvalue. Lembre-se que um valor de p indica a probabilidade de uma dobra mudar tão forte como a observada, 
#ou ainda mais forte, seria vista sob a situação descrita pela hipótese nula. Podemos também resumir os resultados com a seguinte linha de código, que reporta alguma informação adicional, 
#que será coberta em secções posteriores.

summary(res)

#Note-se que existem muitos genes com expressão diferencial devido ao tratamento com dexametasona ao nível de FDR de 10%. No entanto, há duas maneiras de ser mais rigoroso sobre que 
#conjunto de genes são considerados significativos: baixar o limiar da falsa taxa de descoberta (o limiar do padj na tabela de resultados) aumentar o limiar de mudança de log2 fold de 0 
#usando o argumento de resultados lfcThreshold. Se baixarmos o limiar da taxa de falsa descoberta, devemos também informar a função de resultados() sobre ela, para que a função possa 
#utilizar este limiar para a filtragem independente óptima que realiza:

hist(dea$pvalue, breaks=20,col = "grey", border = "white", xlab = "P-value",
     ylab = "Number of genes", main = "P-value value distribution") #hist para vizualizar o p.value

pvalue_fi= sum(dea$pvalue < 0.05, na.rm=TRUE)#filtragem por pvalue normal
pvalue_fi

#Existem 23323 genes com um valor de p inferior a 0,05 entre os 32648 genes para os quais o teste conseguiu reportar um valor de p:
  
#Agora, assumir por um momento que a hipótese nula é verdadeira para todos os genes. Então, pela definição do valor p, 
#esperamos que até 5% dos genes tenham um valor p inferior a 0,05. Isto equivale a 1632 genes. Se considerarmos apenas a lista de genes com um valor de p inferior a 0,05 como 
#diferentemente expresso, esta lista deverá, portanto, conter até 1632 / 23323 = 6% de falsos positivos.
#Assim, se considerarmos aceitável uma fracção de 1 percent de falsos positivos, podemos considerar todos os genes com um valor p ajustado inferior a 1% = 0,01 como significativos. 

hist(dea$padj, breaks=20,col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution") #hist para vizualizar o padj

padj_fi= sum(dea$padj < 0.01, na.rm=TRUE) #probabilidade inferior a 1 % de encontrar um falso positivo 
padj_fi

#ficamos então com 20126 genes



#analise de gráfico:


topGene <- rownames(res)[which.min(res$padj)] #gene com maior expressão
plotCounts(ddsSE_norm, gene = topGene, intgroup=c("paper_IDH.status")) #Gráfico que mostra as contagens do gene com maior expressão para a condição Mutante, e para a condição normal

install.packages("ggbeeswarm")
library("ggbeeswarm")

geneCounts <- plotCounts(ddsSE_norm, gene = topGene, intgroup = c("paper_IDH.status","vital_status"), # relacionar a mutação com o facto de estar vivo ou morto
                         returnData = TRUE)
ggplot(geneCounts, aes(x = paper_IDH.status, y = count, color = vital_status)) + 
  scale_y_log10() +  geom_beeswarm(cex = 3)  # criar o gráfico para vizualizar melhor. podemos explorar isto muito melhor

install.packages("apeglm")
library("apeglm")

# com todos os genes:
plotMA(ddsSE_norm, ylim = c(-10, 10)) #pontos a azul são diferencialmente expressos, a cinza não são. vai de encontro aos nossos testes


with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="green", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="green")
})

#apenas com o top 100
top_100_genes <- res[order(dea$padj)[1:100],]
plotMA(top_100_genes, main="DESeq2", ylim=c(-10,10))

with(top_100_genes[topGene, ], {
  points(baseMean, log2FoldChange, col="green", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="green")
})

#heatmap:
metadata_matriz_clean$paper_IDH.status

sem_NA2=na.omit(metadata_matriz_clean$paper_IDH.status)

sem_NA2

top_20_genes <- res[order(dea$padj)[1:20],]
nova_tabela <- subset(top_20_genes, select = c("padj", "pvalue"))
nova_tabela
nova_matriz <- data.matrix(nova_tabela)
nova_matriz

# Replace missing values with the minimum value in the matrix
min_val <- min(nova_matriz, na.rm = TRUE)
nova_matriz[is.na(nova_matriz)] <- min_val

# Create heatmap
heatmap(nova_matriz, Colv = NA, labCol = sem_NA2)


