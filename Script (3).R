library(bipartite)
library(iNEXT)
library(vegan)
library(ggplot2)
library(rethinking)
source("intramoddist.R")
library(lubridate)
library(tidyr)
library(reshape2)
library(ape)
library(gridExtra)
library(readxl)
library(V.PhyloMaker2)
library(picante)
library(geiger)
library(R.utils)
library(phytools)
library(colorspace)
library(ggridges)
library(dplyr)
library(forcats)
library(viridis)
library(RColorBrewer)
library(circlize)
library(ggbeeswarm)
library(forcats)
library(MASS)
library(emmeans)


DicnodesDic <- read.csv2("dic_MetA.csv", fileEncoding = "Windows-1252", stringsAsFactors = T)
Dic$Visitor <- gsub(pattern = "Hylephila phyleus1",
                                 replacement = "Hylephila phyleus", x = Dic$Visitor)
Dic$Visitor <- gsub(pattern = "Hylephila phyleus0",
                                 replacement = "Hylephila phyleus", x = Dic$Visitor)
Dic$Visitor <- gsub(pattern = "Chlosyne saundersii",
                    replacement = "Chlosyne lacinia saundersi", x = Dic$Visitor)
Dic$Visitor <- gsub(pattern = "Eueides dianasa",
                    replacement = "Eueides iasabella dianasa", x = Dic$Visitor)
Dic$Visitor <- gsub(pattern = "Lasaia sp.",
                    replacement = "Riodinidae sp.", x = Dic$Visitor)
Dic$Visitor <- gsub(pattern = "Urbamus proteus",
                    replacement = "Urbamus proteus proteus", x = Dic$Visitor)
#importando a tabela
Data.raw <- read.csv2("Mestrado Maringá (MetA).csv", fileEncoding = "Windows-1252", stringsAsFactors = T)
table(Data.raw$Plant_species)
summary(Data.raw)

#filtrar/remover 0
Data <- Data.raw [!Data.raw$Insect_species %in% c("0",""),]
Data$Insect_species <- droplevels(Data$Insect_species)
summary(Data)
table(Data$Insect_species)

# alterando nomes do dicionario

setdiff(Data$Insect_species, Dic$Morfo)
setdiff(Dic$Morfo, Data$Insect_species)

Data$Insect_species <- Dic$Visitor[match(Data$Insect_species, Dic$Morfo)]

table(Data$Insect_species)

Data$Month <- unlist(lapply(strsplit(as.character(Data$Date), split = " "),
                            FUN = function(x) x[2]))  #Aqui indica que é o segundo termo de Date e cria uma coluna.
Data$Month <- factor(Data$Month, levels = c("January", "February", "March", "April",
                                            "May", "June", "July", "August",
                                            "September", "October", "November", "December"))

## Atributos
Plants <- read.csv2("Plantas Mestrado Luis - Copia.csv", fileEncoding = "Windows-1252", stringsAsFactors = T, row.names = 1)
Mean_ab_comp <- read.csv2("Abertura e comprimento_novo.csv", fileEncoding = "Windows-1252", stringsAsFactors = T)
Mean_flowers <- aggregate(cbind(comprimento) ~ Espécie , data = Mean_ab_comp, FUN = mean)
match(rownames(Plants), Mean_flowers$Espécie)
Plants <- cbind(Plants, Mean_flowers[match(rownames(Plants), Mean_flowers$Espécie),2])
match(Mean_flowers$Espécie, rownames(Plants))
setdiff(rownames(Plants), Mean_flowers$Espécie)
setdiff(Mean_flowers$Espécie, rownames(Plants))
colnames(Plants)[2] <- "comprimento"

Plants <- Plants[rownames(Plants) %in% rownames(Web),]
Plants$Recurso_fornecido <- as.factor(Plants$Recurso_fornecido)

setdiff(rownames(Plants), rownames(Web))
setdiff(rownames(Web), rownames(Plants))

Data$Plant_species <- as.factor(Data$Plant_species)
Data$Insect_species <- as.factor(Data$Insect_species)

# Remova os níveis que não estão sendo usados
Data$Plant_species <- droplevels(Data$Plant_species)
Data$Insect_species <- droplevels(Data$Insect_species)
## Matrizes de interação

Web <- as.matrix(table(droplevels(Data$Plant_species), droplevels(Data$Insect_species)))
Web <- matrix(Web, ncol = ncol(Web), dimnames = dimnames(Web))
write.csv(Web, file = "Web.csv", row.names = TRUE)

#Soma as interações da planta em Web
#Dai adiciona em Plants uma coluna chamada Int e puxa Int

rede <- Web[,-1]
rede <- as.matrix(rede)
rede <- append(rede, list(rede))

Data.legit <- Data[Data$Interaction == 1, ]
Web.legit <- table(droplevels(Data.legit$Plant_species), droplevels(Data.legit$Insect_species))
Web.legit <- matrix(Web.legit, ncol = ncol(Web.legit), dimnames = dimnames(Web.legit))

soma_interacoes.legit <- rowSums(Web.legit)
Plants$Int.legit <- soma_interacoes.legit[match(rownames(Plants), rownames(Web.legit))]

Data.pill <- Data[Data$Interaction == 2, ]
Web.pill <- table(droplevels(Data.pill$Plant_species), droplevels(Data.pill$Insect_species))
Web.pill <- matrix(Web.pill, ncol = ncol(Web.pill), dimnames = dimnames(Web.pill))

write.csv(Web.legit, file = "Web_legit.csv", row.names = TRUE)

# Salvar Web.pill como CSV
write.csv(Web.pill, file = "Web_pill.csv", row.names = TRUE)

plotweb(Web)
plotweb(Web.legit)
plotweb(Web.pill)
colnames(Web)


########################### Number of interactions ##################

barplot(sort(apply(Web, 2, sum), decreasing = T))
sort(apply(Web, 2, sum), decreasing = T)
sort(apply(Web.legit, 2, sum), decreasing = T)
sort(apply(Web.pill, 2, sum), decreasing = T)

barplot(sort(apply(Web, 1, sum), decreasing = T))
sort(apply(Web, 1, sum), decreasing = T)
sort(apply(Web.legit, 1, sum), decreasing = T)
sort(apply(Web.pill, 1, sum), decreasing = T)

##################################### Modularity ###################

Modulos <- metaComputeModules(Web, N=10)
Modulos_legit <- metaComputeModules(Web.legit, N=10)
Modulos_pill <- metaComputeModules(Web.pill, N=10)


listModuleInformation(Modulos)
listModuleInformation(Modulos_legit)
listModuleInformation(Modulos_pill)

modulos.p <- vector("list", length = 10)
modulos.a <- vector("list", length = 10)

modulos.p.legit<-vector("list", length = 8)
modulos.a.legit <- vector("list", length = 8)

modulos.p.pill<-vector("list", length = 7)
modulos.a.pill <- vector("list", length = 7)

for (i in 1:10){
  modulos.p[[i]] <- listModuleInformation(Modulos)[[2]][[i]][[1]]
}
for (i in 1:10){
  modulos.a[[i]] <- listModuleInformation(Modulos)[[2]][[i]][[2]]
}

for (i in 1:8){
  modulos.p.legit[[i]] <- listModuleInformation(Modulos_legit)[[2]][[i]][[1]]
}
for (i in 1:8){
  modulos.a.legit[[i]] <- listModuleInformation(Modulos_legit)[[2]][[i]][[2]]
}

for (i in 1:7){
  modulos.p.pill[[i]] <- listModuleInformation(Modulos_pill)[[2]][[i]][[1]]
}
for (i in 1:7){
  modulos.a.pill[[i]] <- listModuleInformation(Modulos_pill)[[2]][[i]][[2]]
}

Plantas.Mod <- data.frame(sp.planta = unlist(modulos.p),
                      modulo = rep(1:10, unlist(lapply(modulos.p, length))))

Insetos.Mod <- data.frame(sp.inseto = unlist(modulos.a),
                          modulo = rep(1:10,
                                       unlist(lapply(modulos.a, length))))

Plantas.Mod.legit <- data.frame(sp.planta = unlist(modulos.p.legit),
                          modulo = rep(1:8, unlist(lapply(modulos.p.legit, length))))

Insetos.Mod.legit <- data.frame(sp.inseto = unlist(modulos.a.legit),
                          modulo = rep(1:8, unlist(lapply(modulos.a.legit, length))))

Plantas.Mod.pill <- data.frame(sp.planta = unlist(modulos.p.pill),
                                modulo = rep(1:7, unlist(lapply(modulos.p.pill, length))))

Insetos.Mod.pill <- data.frame(sp.inseto = unlist(modulos.a.pill),
                                modulo = rep(1:7, unlist(lapply(modulos.a.pill, length))))


###################### Overall network #################################
plotModuleWeb(Modulos)
printoutModuleInformation(Modulos)

CZ_Polin <- czvalues(Modulos, level="higher")
CZ_Plants <- czvalues(Modulos, level="lower")
df_Web <- as.data.frame(CZ_Plants)

nulls<- vaznull(web = Web, N = 500)
module.nulls <- sapply(nulls, metaComputeModules, N = 10)

like.nulls<- sapply(module.nulls, function(x) x@likelihood)
quantile(unlist(like.nulls), 0.95)
rethinking::dens(like.nulls, xlim = c(0, .5))
abline(v = Modulos@likelihood, col = 2)

Reps.over<-500
P.value <- min(c(sum(c(like.nulls,Modulos@likelihood) >= Modulos@likelihood)/(Resps.over+1),
                 sum(c(like.nulls,Modulos@likelihood) <= Modulos@likelihood)/(Resps.over+1)))*2

czvalues(Modulos, level="higher")
czvalues(Modulos, level="lower")
Valores_Web<- czvalues(Modulos, level="lower")
df_Web <- as.data.frame(Valores_Web)
print(df_Web)


null.cz.h <- lapply(module.nulls, czvalues, level="higher")
null.cz.l <- lapply(module.nulls, czvalues, level="lower")


null.cs.h <- sapply(null.cz.h, function(x) x$c)
null.cs.l <- sapply(null.cz.l, function(x) x$c)

null.zs.h <- sapply(null.cz.h, function(x) x$z)
null.zs.l <- sapply(null.cz.l, function(x) x$z)

Z_critical_l <- quantile(null.zs.l, 0.95, na.rm=T)

C_critical_l <- quantile(null.cs.l, 0.95, na.rm=T)

df_Web$role <- NA
df_Web$role[df_Web$c >= C_critical_l & df_Web$z >= Z_critical_l] <- "Hub"
df_Web$role[df_Web$c >= C_critical_l & df_Web$z <= Z_critical_l] <- "Connector"
df_Web$role[df_Web$c <= C_critical_l & df_Web$z >= Z_critical_l] <- "Module Hub"
df_Web$role[df_Web$c <= C_critical_l & df_Web$z <= Z_critical_l] <- "Peripheral"
df_Web$role <- as.factor(df_Web$role)
df_web_clean <- df_Web[complete.cases(df_Web), ]
par(mar = c(4,4,1,1))

svg("Speciesrole.web.1.svg")
plot(df_web_clean$c, df_web_clean$z, col = df_web_clean$role, pch = 16, xlab = "C",
     ylab = "Z")
abline(h = Z_critical_l, lty = 3)
abline(v = C_critical_l, lty = 3)
dev.off()
Hubs_Web <- rownames(df_web_clean)[df_web_clean$role == "Hub"]
Connectors_Web <- rownames(df_web_clean)[df_web_clean$role == "Connector"]
ModHubs_Web <- rownames(df_web_clean)[df_web_clean$role == "Module Hub"]
Periph_Web <- rownames(df_web_clean)[df_web_clean$role == "Peripheral"]

Plant.dist[rownames(Plant.dist) %in% Connectors_Web, colnames(Plant.dist) %in% Connectors_Web]
Plant.dist[rownames(Plant.dist) %in% ModHubs_Web, colnames(Plant.dist) %in% ModHubs_Web]

Modulos_SPs <- lapply(listModuleInformation(Modulos)[[2]], FUN = function(x)x[[1]])


df_Web$module <- NA
for(i in 1:length(Modulos_SPs)){
  df_Web$module[rownames(df_Web) %in% Modulos_SPs[[i]]] <- i
}
svg("Speciesrole.svg")
par(mfrow = c(3, 1))
par(mar = c(4,4,1,1))
svg("Speciesrole.web.2 modl.svg")
plot(df_Web$c, df_Web$z, col = df_Web$module, pch = 16, xlab = "C",
     ylab = "Z", cex = 1.3)
abline(h = Z_critical_l, lty = 3)
abline(v = C_critical_l, lty = 3)
dev.off()

###################################### Legit network ############################

plotModuleWeb(Modulos_legit)
printoutModuleInformation(Modulos_legit)

CZ_Plants_legit <- czvalues(Modulos_legit, level="lower")
df_Web_legit <- as.data.frame(CZ_Plants_legit)

nulls_legit<- vaznull(web = Web.legit, N = 500)
module.nulls_legit <- sapply(nulls_legit, metaComputeModules, N =10)

like.nulls_legit<- sapply(module.nulls_legit, function(x) x@likelihood)
quantile(unlist(like.nulls_legit), 0.95)
rethinking::dens(like.nulls_legit, xlim = c(0, .5))
abline(v = Modulos_legit@likelihood, col = 2)

Reps.legit <- 500
P.value.legit <- min(c(sum(c(like.nulls_legit,Modulos_legit@likelihood) >= Modulos_legit@likelihood)/(Reps.legit+1),
                      sum(c(like.nulls_legit,Modulos_legit@likelihood) <= Modulos_legit@likelihood)/(Reps.legit+1)))*2

czvalues(Modulos_legit, level="higher")
czvalues(Modulos_legit, level="lower")

Valores_Web_legit<- czvalues(Modulos_legit, level="lower")
df_Web_legit <- as.data.frame(Valores_Web_legit)
print(df_Web_legit)

null.cz.l_legit <- lapply(module.nulls_legit, czvalues, level="lower")

# compute 95% CI for c and z:

# c-values across all species in nulls
null_web.cs.l_legit <- sapply(null.cz.l_legit, function(x) x$c)

#

null_web.zs.l_legit <- sapply(null.cz.l_legit, function(x) x$z)

C_critical_l_legit <- quantile(null_web.cs.l_legit, 0.95, na.rm=T)
Z_critical_l_legit <- quantile(null_web.zs.l_legit, 0.95, na.rm=T)

df_Web_legit$role <- NA
df_Web_legit$role[df_Web_legit$c > C_critical_l_legit & df_Web_legit$z > Z_critical_l_legit] <- "Hub"
df_Web_legit$role[df_Web_legit$c > C_critical_l_legit & df_Web_legit$z < Z_critical_l_legit] <- "Connector"
df_Web_legit$role[df_Web_legit$c < C_critical_l_legit & df_Web_legit$z > Z_critical_l_legit] <- "Module Hub"
df_Web_legit$role[df_Web_legit$c < C_critical_l_legit & df_Web_legit$z < Z_critical_l_legit] <- "Peripheral"
df_Web_legit$role <- as.factor(df_Web_legit$role)
df_web_clean_legit <- df_Web_legit[complete.cases(df_Web_legit), ]

par(mar = c(4,4,1,1))
svg("Species_role.web_legit.2 modl.svg")
plot(df_web_clean_legit$c, df_web_clean_legit$z, col = df_web_clean_legit$role, pch = 16, xlab = "C",
     ylab = "Z")
abline(h = Z_critical_l_legit, lty = 3)
abline(v = C_critical_l_legit, lty = 3)
dev.off()

Connectors_Web_legit <- rownames(df_web_clean_legit)[df_web_clean_legit$role == "Connector"]
ModHubs_Web_legit <- rownames(df_web_clean_legit)[df_web_clean_legit$role == "Module Hub"]
Periph_Web_legit <- rownames(df_web_clean_legit)[df_web_clean_legit$role == "Peripheral"]

Plant.dist[rownames(Plant.dist) %in% Connectors_Web_legit, colnames(Plant.dist) %in% Connectors_Web_legit]

Plant.dist[rownames(Plant.dist) %in% ModHubs_Web_legit, colnames(Plant.dist) %in% ModHubs_Web_legit]

###################################### Plot c and z legit network###############################
Modulos_SPs_legit <- lapply(listModuleInformation(Modulos_legit)[[2]], FUN = function(x)x[[1]])
df_web_clean_legit$module <- NA
for(i in 1:length(Modulos_SPs_legit)){
  df_web_clean_legit$module[rownames(df_web_clean_legit) %in% Modulos_SPs_legit[[i]]] <- i
}

par(mar = c(4,4,1,1))
svg("Species_role.web_legit.2 modl.svg")
plot(df_web_clean_legit$c, df_web_clean_legit$z, col = df_web_clean_legit$module, pch = 16, xlab = "C",
     ylab = "Z", cex = 1.3)
abline(h = Z_critical_l_legit, lty = 3)
abline(v = C_critical_l_legit, lty = 3)
dev.off()

##################################### Illegit network ############################################

plotModuleWeb(Modulos_pill)
printoutModuleInformation(Modulos_pill)

CZ_Polin_pill <- czvalues(Modulos_pill, level="higher")
CZ_Plants_pill <- czvalues(Modulos_pill, level="lower")
df_Web_pill <- as.data.frame(CZ_Plants_pill)


nulls_pill<- vaznull(web = Web.pill, N = 500)
module.nulls_pill <- sapply(nulls_pill, metaComputeModules, N=10)

like.nulls_pill<- sapply(module.nulls_pill, function(x) x@likelihood)
quantile(unlist(like.nulls_pill), 0.95)
rethinking::dens(like.nulls_pill, xlim = c(0, .5))
abline(v = Modulos_pill@likelihood, col = 2)

Reps.pill <- 500
P.value.pill <- min(c(sum(c(like.nulls_pill,Modulos_pill@likelihood) >= Modulos_pill@likelihood)/(Reps.pill+1),
                 sum(c(like.nulls_pill,Modulos_pill@likelihood) <= Modulos_pill@likelihood)/(Reps.pill+1)))*2

P.value.pill*500

czvalues(Modulos_pill, level="higher")
czvalues(Modulos_pill, level="lower")
Valores_Web_pill<- czvalues(Modulos_pill, level="lower")
df_Web_pill <- as.data.frame(Valores_Web_pill)
print(df_Web_pill)


null.cz.h_pill <- lapply(module.nulls_pill, czvalues, level="higher")
null.cz.l_pill <- lapply(module.nulls_pill, czvalues, level="lower")

# compute 95% CI for c and z:


null.cs.h_pill <- sapply(null.cz.h_pill, function(x) x$c) # c-values across all species in nulls
null.cs.l_pill <- sapply(null.cz.l_pill, function(x) x$c)

#
null.zs.h_pill <- sapply(null.cz.h_pill, function(x) x$z)
null.zs.l_pill <- sapply(null.cz.l_pill, function(x) x$z)

Z_critical_l_pill <- quantile(null.zs.l_pill, 0.95, na.rm=T)

C_critical_l_pill <- quantile(null.cs.l_pill, 0.95, na.rm=T)

df_Web_pill$role <- NA
df_Web_pill$role[df_Web_pill$c >= C_critical_l_pill & df_Web_pill$z >= Z_critical_l_pill] <- "Hub"
df_Web_pill$role[df_Web_pill$c >= C_critical_l_pill & df_Web_pill$z <= Z_critical_l_pill] <- "Connector"
df_Web_pill$role[df_Web_pill$c <= C_critical_l_pill & df_Web_pill$z >= Z_critical_l_pill] <- "Module Hub"
df_Web_pill$role[df_Web_pill$c <= C_critical_l_pill & df_Web_pill$z <= Z_critical_l_pill] <- "Peripheral"
df_Web_pill$role <- as.factor(df_Web_pill$role)
df_Web_pill$role <- as.factor(df_Web_pill$role)
df_web_clean_pill <- df_Web_pill[complete.cases(df_Web_pill[, c("c", "z", "role")]), ]
par(mar = c(4,4,1,1))
svg("Species_role.web_pill.2 modl.svg")
plot(df_web_clean_pill$c, df_web_clean_pill$z, col = df_web_clean_pill$role, pch = 16, xlab = "C",
     ylab = "Z")
abline(h = Z_critical_l_pill, lty = 3)
abline(v = C_critical_l_pill, lty = 3)
dev.off()

Hubs_Web_pill <- rownames(df_web_clean_pill)[df_web_clean_pill$role == "Hub"]
Connectors_Web_pill <- rownames(df_web_clean_pill)[df_web_clean_pill$role == "Connector"]
ModHubs_Web_pill <- rownames(df_web_clean_pill)[df_web_clean_pill$role == "Module Hub"]
Periph_Web_pill <- rownames(df_web_clean_pill)[df_web_clean_pill$role == "Peripheral"]
plotModuleWeb(Modulos_pill)

Plant.dist[rownames(Plant.dist) %in% ModHubs_Web_pill, colnames(Plant.dist) %in% Periph_Web_pill]

Modulos_SPs_pill <- lapply(listModuleInformation(Modulos_pill)[[2]], FUN = function(x)x[[1]])
df_web_clean_pill$module <- NA
for(i in 1:length(Modulos_SPs_pill)){
  df_web_clean_pill$module[rownames(df_web_clean_pill) %in% Modulos_SPs_pill[[i]]] <- i
}

par(mar = c(4,4,1,1))
svg("Species_role.web_pill.2 modl.svg")
plot(df_web_clean_pill$c, df_web_clean_pill$z, col = df_web_clean_pill$module, pch = 16, xlab = "C",
     ylab = "Z", cex = 1.3)
abline(h = Z_critical_l_pill, lty = 3)
abline(v = C_critical_l_pill, lty = 3)
dev.off()

#################### Phylogeny #########################

Regional <- read_xlsx("Regional.xlsx")
Filo.Regional <- phylo.maker(sp.list=Regional)
Filo.Regional <- Filo.Regional$scenario.3
View(Filo.Regional)
plot.phylo(Filo.Regional, use.edge.length = TRUE, cex = 0.8)

tree_lido <- read.nexus("arvore_Phylo_Regional.nex")
file.rename("arvore_Phylo_Regional.nex", "Perugini_SX_phylogeny.nex")
plot.phylo(tree_lido)

Filo.dist <- cophenetic(Filo.Regional)
rownames(Filo.dist) <- gsub(pattern = "_", replacement = " ", x = rownames(Filo.dist))
colnames(Filo.dist) <- gsub(pattern = "_", replacement = " ", x = colnames(Filo.dist))


Filo.dist.legit <- Filo.dist[rownames(Filo.dist) %in% rownames(Web.legit),colnames(Filo.dist) %in% rownames(Web.legit)]

Filo.dist.pill <- Filo.dist[rownames(Filo.dist) %in% rownames(Web.pill),colnames(Filo.dist) %in% rownames(Web.pill)]

################ Plant distances ####################

Plant.dist <- as.matrix(cluster::daisy(Plants, metric = "gower"))
Plant.dist_legit <- as.matrix(cluster::daisy(Plants, metric = "gower"))
Plant.dist_pill <- as.matrix(cluster::daisy(Plants, metric = "gower"))

################################# Species roles
####### Morpho similarity ##################
#ModHub_WEB

Dist_ModHubs_Web <- mean(Plant.dist[ModHubs_Web, ModHubs_Web])
null_ModHubs_Web <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_ModHubs_Web)){
  Aleatorio <- sample(1:nrow(Plant.dist), size = length(ModHubs_Web))
  null_ModHubs_Web[i] <- mean(Plant.dist[Aleatorio, Aleatorio])
}
quantile(null_ModHubs_Web, probs = c(0.025, 0.975))
rethinking::dens(null_ModHubs_Web, xlim = c(0, 1))
abline(v = Dist_ModHubs_Web, col = 2, lwd = 2)

Reps.ModHubs_Web<-10000
P.value.ModHubs_Web <- min(c(sum(c(null_ModHubs_Web,Dist_ModHubs_Web) >= Dist_ModHubs_Web)/(Reps.ModHubs_Web+1),
                 sum(c(null_ModHubs_Web,Dist_ModHubs_Web) <= Dist_ModHubs_Web)/(Reps.ModHubs_Web+1)))*2

#Connectors_Web
Dist_Connectors_Web <- mean(Plant.dist[Connectors_Web, Connectors_Web])
null_Connectors_Web <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Connectors_Web)){
  Aleatorio <- sample(1:nrow(Plant.dist), size = length(Connectors_Web))
  null_Connectors_Web[i] <- mean(Plant.dist[Aleatorio, Aleatorio])
}
quantile(null_Connectors_Web, probs = c(0.025, 0.975))
rethinking::dens(null_Connectors_Web, xlim = c(0, 1))
abline(v = Dist_Connectors_Web, col = 2, lwd = 2)

Reps.Connectors_Web<-10000
P.value.Connectors_Web <- min(c(sum(c(null_Connectors_Web,Dist_Connectors_Web) >= Dist_Connectors_Web)/(Reps.Connectors_Web+1),
                             sum(c(null_Connectors_Web,Dist_Connectors_Web) <= Dist_Connectors_Web)/(Reps.Connectors_Web+1)))*2

#Periferic_Web
Dist_Peri_Web <- mean(Plant.dist[Periph_Web, Periph_Web])
null_Periph_Web <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Periph_Web)){
  Aleatorio <- sample(1:nrow(Plant.dist), size = length(Periph_Web))
  null_Periph_Web[i] <- mean(Plant.dist[Aleatorio, Aleatorio])
}
quantile(null_Periph_Web, probs = c(0.025, 0.975))
rethinking::dens(null_Periph_Web, xlim = c(0, 1))
abline(v = Dist_Peri_Web, col = 2, lwd = 2)

Reps.Peri_Web<-10000
P.value.Peri_Web <- min(c(sum(c(null_Periph_Web,Dist_Peri_Web) >= Dist_Peri_Web)/(Reps.Peri_Web+1),
                                sum(c(null_Periph_Web,Dist_Peri_Web) <= Dist_Peri_Web)/(Reps.Peri_Web+1)))*2

#ModHub de legit


Dist_ModHubs_Web_legit <- mean(Plant.dist_legit[ModHubs_Web_legit, ModHubs_Web_legit])
null_ModHubs_Web_legit <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_ModHubs_Web_legit)){
  Aleatorio <- sample(1:nrow(Plant.dist_legit), size = length(ModHubs_Web_legit))
  null_ModHubs_Web_legit[i] <- mean(Plant.dist_legit[Aleatorio, Aleatorio])
}
quantile(null_ModHubs_Web_legit, probs = c(0.025, 0.975))
rethinking::dens(null_ModHubs_Web_legit, xlim = c(0, 1))
abline(v = Dist_ModHubs_Web_legit, col = 2, lwd = 2)

Reps.ModHubs_Web_legit<-10000
P.value.ModHubs_Web_legit <- min(c(sum(c(null_ModHubs_Web_legit,Dist_ModHubs_Web_legit) >= Dist_ModHubs_Web_legit)/(Reps.ModHubs_Web_legit+1),
                          sum(c(null_ModHubs_Web_legit,Dist_ModHubs_Web_legit) <= Dist_ModHubs_Web_legit)/(Reps.ModHubs_Web_legit+1)))*2

#Connectors_Web_legit

Dist_Connectors_Web_legit <- mean(Plant.dist_legit[Connectors_Web_legit, Connectors_Web_legit])

null_Connectors_Web_legit <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Connectors_Web_legit)){
  Aleatorio <- sample(1:nrow(Plant.dist_legit), size = length(Connectors_Web_legit))
  null_Connectors_Web_legit[i] <- mean(Plant.dist_legit[Aleatorio, Aleatorio])
}
quantile(null_Connectors_Web_legit, probs = c(0.025, 0.975))
rethinking::dens(null_Connectors_Web_legit, xlim = c(0, 1))
abline(v = Dist_Connectors_Web_legit, col = 2, lwd = 2)

Reps.Connectors_Web_legit<-10000
P.value.Connectors_Web_legit <- min(c(sum(c(null_Connectors_Web_legit,Dist_Connectors_Web_legit) >= Dist_Connectors_Web_legit)/(Reps.Connectors_Web_legit+1),
                                   sum(c(null_Connectors_Web_legit,Dist_Connectors_Web_legit) <= Dist_Connectors_Web_legit)/(Reps.Connectors_Web_legit+1)))*2

#Periferic_Web_legit

Dist_Peri_Web_legit <- mean(Plant.dist_legit[Periph_Web_legit, Periph_Web_legit])
null_Periph_Web_legit <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Periph_Web_legit)){
  Aleatorio <- sample(1:nrow(Plant.dist_legit), size = length(Periph_Web_legit))
  null_Periph_Web_legit[i] <- mean(Plant.dist_legit[Aleatorio, Aleatorio])
}
quantile(null_Periph_Web_legit, probs = c(0.025, 0.975))
rethinking::dens(null_Periph_Web_legit, xlim = c(0, 1))
abline(v = Dist_Peri_Web_legit, col = 2, lwd = 2)

Reps.Peri_Web_legit<-10000
P.value.Peri_Web_legit <- min(c(sum(c(null_Periph_Web_legit,Dist_Peri_Web_legit) >= Dist_Peri_Web_legit)/(Reps.Peri_Web_legit+1),
                                      sum(c(null_Periph_Web_legit,Dist_Peri_Web_legit) <= Dist_Peri_Web_legit)/(Reps.Peri_Web_legit+1)))*2

#ModHubs_Web_pill

Plant.dist_pill <- as.matrix(cluster::daisy(Plants, metric = "gower"))

Dist_ModHubs_Web_pill <- mean(Plant.dist_pill[ModHubs_Web_pill, ModHubs_Web_pill])
null_ModHubs_Web_pill <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_ModHubs_Web_pill)){
  Aleatorio <- sample(1:nrow(Plant.dist), size = length(ModHubs_Web_pill))
  null_ModHubs_Web_pill[i] <- mean(Plant.dist_pill[Aleatorio, Aleatorio])
}
quantile(null_ModHubs_Web_pill, probs = c(0.025, 0.975))
rethinking::dens(null_ModHubs_Web_pill, xlim = c(0, 1))
abline(v = Dist_ModHubs_Web_pill, col = 2, lwd = 2)

Reps.ModHubs_Web_pill<-10000
P.value.ModHubs_Web_pill <- min(c(sum(c(null_ModHubs_Web_pill,Dist_ModHubs_Web_pill) >= Dist_ModHubs_Web_pill)/(Reps.ModHubs_Web_pill+1),
                                sum(c(null_ModHubs_Web_pill,Dist_ModHubs_Web_pill) <= Dist_ModHubs_Web_pill)/(Reps.ModHubs_Web_pill+1)))*2

#Connec_Web_pill
Dist_Connec_Web_pill <- mean(Plant.dist_pill[Connectors_Web_pill, Connectors_Web_pill])
null_Connec_Web_pill <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Connec_Web_pill)){
  Aleatorio <- sample(1:nrow(Plant.dist), size = length(Connectors_Web_pill))
  null_Connec_Web_pill[i] <- mean(Plant.dist_pill[Aleatorio, Aleatorio])
}
quantile(null_Connec_Web_pill, probs = c(0.025, 0.975))
rethinking::dens(null_Connec_Web_pill, xlim = c(0, 1))
abline(v = Dist_Connec_Web_pill, col = 2, lwd = 2)

Reps.Connec_Web_pill<-10000
P.value.Connec_Web_pill <- min(c(sum(c(null_Connec_Web_pill,Dist_Connec_Web_pill) >= Dist_Connec_Web_pill)/(Reps.Connec_Web_pill+1),
                                  sum(c(null_Connec_Web_pill,Dist_Connec_Web_pill) <= Dist_Connec_Web_pill)/(Reps.Connec_Web_pill+1)))*2

#Periph_Web_pill

Dist_Peri_Web_pill <- mean(Plant.dist_pill[Periph_Web_pill, Periph_Web_pill])
null_Periph_Web_pill <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Periph_Web_pill)){
  Aleatorio <- sample(1:nrow(Plant.dist_pill), size = length(Periph_Web_pill))
  null_Periph_Web_pill[i] <- mean(Plant.dist_pill[Aleatorio, Aleatorio])
}
quantile(null_Periph_Web_pill, probs = c(0.025, 0.975))
rethinking::dens(null_Periph_Web_pill, xlim = c(0, 1))
abline(v = Dist_Peri_Web_pill, col = 2, lwd = 2)

Reps.Peri_Web_pill<-10000
P.value.Peri_Web_pill <- min(c(sum(c(null_Periph_Web_pill,Dist_Peri_Web_pill) >= Dist_Peri_Web_pill)/(Reps.Peri_Web_pill+1),
                                 sum(c(null_Periph_Web_pill,Dist_Peri_Web_pill) <= Dist_Peri_Web_pill)/(Reps.Peri_Web_pill+1)))*2

################################# Phylogeny and species role #################################

################################# Overall ###########################################
#Connector_Web_filo_global
Filo_Connec_Web_global <- mean(Filo.dist[Connectors_Web, Connectors_Web])
null_Filo_connec_Web_global <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Filo_connec_Web_global)){
  Aleatorio <- sample(1:nrow(Filo.dist), size = length(Connectors_Web))
  null_Filo_connec_Web_global[i] <- mean(Filo.dist[Aleatorio, Aleatorio])
}
quantile(null_Filo_connec_Web_global, probs = c(0.025, 0.975))
rethinking::dens(null_Filo_connec_Web_global, xlim = c(0, 300))
abline(v = Filo_Connec_Web_global, col = 2, lwd = 2)

Reps._Filo.Connec_Web_global<-10000
P.value.Filo_Connec_Web_global <- min(c(sum(c(null_Filo_connec_Web_global,Filo_Connec_Web_global) >= Filo_Connec_Web_global)/(Reps.Connec_Web_global+1),
                               sum(c(null_Filo_connec_Web_global,Filo_Connec_Web_global) <= Filo_Connec_Web_global)/(Reps.Connec_Web_global+1)))*2

#ModHub_Web_filo_global
Filo_ModHub_Web_global <- mean(Filo.dist[ModHubs_Web, ModHubs_Web])
null_Filo_ModHub_Web_global <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Filo_ModHub_Web_global)){
  Aleatorio <- sample(1:nrow(Filo.dist), size = length(ModHubs_Web))
  null_Filo_ModHub_Web_global[i] <- mean(Filo.dist[Aleatorio, Aleatorio])
}
quantile(null_Filo_ModHub_Web_global, probs = c(0.025, 0.975))
rethinking::dens(null_Filo_ModHub_Web_global, xlim = c(0, 300))
abline(v = Filo_ModHub_Web_global, col = 2, lwd = 2)

Reps._Filo.ModHub_Web_global<-10000
P.value.Filo_Connec_Web_global <- min(c(sum(c(null_Filo_ModHub_Web_global,Filo_ModHub_Web_global) >= Filo_ModHub_Web_global)/(Reps._Filo.ModHub_Web_global+1),
                                   sum(c(null_Filo_ModHub_Web_global,Filo_ModHub_Web_global) <= Filo_ModHub_Web_global)/(Reps._Filo.ModHub_Web_global+1)))*2

#Periph_filo_global
Filo_Periph_Web_global <- mean(Filo.dist[Periph_Web, Periph_Web])
null_Filo_Periph_Web_global <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Filo_Periph_Web_global)){
  Aleatorio <- sample(1:nrow(Filo.dist), size = length(Periph_Web))
  null_Filo_Periph_Web_global[i] <- mean(Filo.dist[Aleatorio, Aleatorio])
}
quantile(null_Filo_Periph_Web_global, probs = c(0.025, 0.975))
rethinking::dens(null_Filo_Periph_Web_global, xlim = c(0, 300))
abline(v = Filo_Periph_Web_global, col = 2, lwd = 2)

Reps._Filo.Periph_Web_global<-10000
P.value.Filo_Periph_Web_global <- min(c(sum(c(null_Filo_Periph_Web_global,Filo_Periph_Web_global) >= Filo_Periph_Web_global)/(Reps._Filo.Periph_Web_global+1),
                                   sum(c(null_Filo_Periph_Web_global,Filo_Periph_Web_global) <= Filo_Periph_Web_global)/(Reps._Filo.Periph_Web_global+1)))*2

############################## Legitimate #########################################
#Filo_Connec_Web_legit
Filo_Connec_Web_legit <- mean(Filo.dist.legit[Connectors_Web_legit, Connectors_Web_legit])
null_Filo_connec_Web_legit <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Filo_connec_Web_legit)){
  Aleatorio <- sample(1:nrow(Filo.dist.legit), size = length(Connectors_Web_legit))
  null_Filo_connec_Web_legit[i] <- mean(Filo.dist.legit[Aleatorio, Aleatorio])
}
quantile(null_Filo_connec_Web_legit, probs = c(0.025, 0.975))
rethinking::dens(null_Filo_connec_Web_legit, xlim = c(0, 300))
abline(v = Filo_Connec_Web_legit, col = 2, lwd = 2)

Reps._Filo.Connec_Web_legit<-10000
P.value.Filo_Connec_Web_legit <- min(c(sum(c(null_Filo_connec_Web_legit,Filo_Connec_Web_legit) >= Filo_Connec_Web_legit)/(Reps._Filo.Connec_Web_legit+1),
                                   sum(c(null_Filo_connec_Web_legit,Filo_Connec_Web_legit) <= Filo_Connec_Web_legit)/(Reps._Filo.Connec_Web_legit+1)))*2

#Filo_ModHub_Web_legit
Filo_ModHub_Web_legit <- mean(Filo.dist.legit[ModHubs_Web_legit, ModHubs_Web_legit])
null_Filo_ModHub_Web_legit <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Filo_ModHub_Web_legit)){
  Aleatorio <- sample(1:nrow(Filo.dist.legit), size = length(ModHubs_Web_legit))
  null_Filo_ModHub_Web_legit[i] <- mean(Filo.dist.legit[Aleatorio, Aleatorio])
}
quantile(null_Filo_ModHub_Web_legit, probs = c(0.025, 0.975))
rethinking::dens(null_Filo_ModHub_Web_legit, xlim = c(0, 300))
abline(v = Filo_ModHub_Web_legit, col = 2, lwd = 2)

Reps._Filo.ModHub_Web_legit<-10000
P.value.Filo_ModHub_Web_legit <- min(c(sum(c(null_Filo_ModHub_Web_legit,Filo_ModHub_Web_legit) >= Filo_ModHub_Web_legit)/(Reps._Filo.ModHub_Web_legit+1),
                                  sum(c(null_Filo_ModHub_Web_legit,Filo_ModHub_Web_legit) <= Filo_ModHub_Web_legit)/(Reps._Filo.ModHub_Web_legit+1)))*2

#Filo_Periph_Web_legit
Filo_Periph_Web_legit <- mean(Filo.dist.legit[Periph_Web_legit, Periph_Web_legit])
null_Filo_Periph_Web_legit <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Filo_Periph_Web_legit)){
  Aleatorio <- sample(1:nrow(Filo.dist.legit), size = length(Periph_Web_legit))
  null_Filo_Periph_Web_legit[i] <- mean(Filo.dist.legit[Aleatorio, Aleatorio])
}
quantile(null_Filo_Periph_Web_legit, probs = c(0.025, 0.975))
rethinking::dens(null_Filo_Periph_Web_legit, xlim = c(0, 300))
abline(v = Filo_Periph_Web_legit, col = 2, lwd = 2)

Reps._Filo_Periph_Web_legit<-10000
P.value.Filo_Periph_Web_legit <- min(c(sum(c(null_Filo_Periph_Web_legit,Filo_Periph_Web_legit) >= Filo_Periph_Web_legit)/(Reps._Filo_Periph_Web_legit+1),
                                       sum(c(null_Filo_Periph_Web_legit,Filo_Periph_Web_legit) <= Filo_Periph_Web_legit)/(Reps._Filo_Periph_Web_legit+1)))*2


#################################### Illegitimate ########################################
#Connector_Web_Pill_Filo
Filo_Connec_Web_pill <- mean(Filo.dist.pill[Connectors_Web_pill, Connectors_Web_pill])
null_Filo_connec_Web_pill <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Filo_connec_Web_pill)){
  Aleatorio <- sample(1:nrow(Filo.dist.pill), size = length(Connectors_Web_pill))
  null_Filo_connec_Web_pill[i] <- mean(Filo.dist.pill[Aleatorio, Aleatorio])
}
quantile(null_Filo_connec_Web_pill, probs = c(0.025, 0.975))
rethinking::dens(null_Filo_connec_Web_pill, xlim = c(0, 300))
abline(v = Filo_Connec_Web_pill, col = 2, lwd = 2)

Reps._Filo_Connec_Web_pill<-10000
P.value.Filo_Connec_Web_pill <- min(c(sum(c(null_Filo_connec_Web_pill,Filo_Connec_Web_pill) >= Filo_Connec_Web_pill)/(Reps._Filo_Connec_Web_pill+1),
                                       sum(c(null_Filo_connec_Web_pill,Filo_Connec_Web_pill) <= Filo_Connec_Web_pill)/(Reps._Filo_Connec_Web_pill+1)))*2

#Modhub_Web_Pill_Filo
Filo_Modhub_Web_pill <- mean(Filo.dist.pill[ModHubs_Web_pill, ModHubs_Web_pill])
null_Filo_modhub_Web_pill <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Filo_modhub_Web_pill)){
  Aleatorio <- sample(1:nrow(Filo.dist.pill), size = length(ModHubs_Web_pill))
  null_Filo_modhub_Web_pill[i] <- mean(Filo.dist.pill[Aleatorio, Aleatorio])
}
quantile(null_Filo_modhub_Web_pill, probs = c(0.025, 0.975))
rethinking::dens(null_Filo_modhub_Web_pill, xlim = c(0, 300))
abline(v = Filo_Modhub_Web_pill, col = 2, lwd = 2)

Reps._Filo_Modhub_Web_pill<-10000
P.value.Filo_Modhub_Web_pill <- min(c(sum(c(null_Filo_modhub_Web_pill,Filo_Modhub_Web_pill) >= Filo_Modhub_Web_pill)/(Reps._Filo_Modhub_Web_pill+1),
                                      sum(c(null_Filo_modhub_Web_pill,Filo_Modhub_Web_pill) <= Filo_Modhub_Web_pill)/(Reps._Filo_Modhub_Web_pill+1)))*2

#Periph_Web_Pill_Filo
Filo_Peri_Web_pill <- mean(Filo.dist.pill[Periph_Web_pill, Periph_Web_pill])
null_Filo_Periph_Web_pill <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Filo_Periph_Web_pill)){
  Aleatorio <- sample(1:nrow(Filo.dist.pill), size = length(Periph_Web_pill))
  null_Filo_Periph_Web_pill[i] <- mean(Filo.dist.pill[Aleatorio, Aleatorio])
}
quantile(null_Filo_Periph_Web_pill, probs = c(0.025, 0.975))
rethinking::dens(null_Filo_Periph_Web_pill, xlim = c(0, 300))
abline(v = Filo_Peri_Web_pill, col = 2, lwd = 2)

Reps._Filo_Peri_Web_pill<-10000
P.value.Filo_Peri_Web_pill <- min(c(sum(c(null_Filo_Periph_Web_pill,Filo_Peri_Web_pill) >= Filo_Peri_Web_pill)/(Reps._Filo_Peri_Web_pill+1),
                                      sum(c(null_Filo_Periph_Web_pill,Filo_Peri_Web_pill) <= Filo_Peri_Web_pill)/(Reps._Filo_Peri_Web_pill+1)))*2

####################################### Phenology e sp. role #####
###Global ###

#Fenol_connec_Web
Fenol_connec_Web <- mean(Fenol.dist_web[Connectors_Web, Connectors_Web])
null_Fenol_connec_Web <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Fenol_connec_Web)){
  Aleatorio <- sample(1:nrow(Fenol.dist_web), size = length(Connectors_Web))
  null_Fenol_connec_Web[i] <- mean(Fenol.dist_web[Aleatorio, Aleatorio])
}
quantile(null_Fenol_connec_Web, probs = c(0.025, 0.975))
rethinking::dens(null_Fenol_connec_Web, xlim = c(0, 1))
abline(v = Fenol_connec_Web, col = 2, lwd = 2)

Reps._Fenol_connec_Web<-10000
P.value.Fenol_connec_Web <- min(c(sum(c(null_Fenol_connec_Web,Fenol_connec_Web) >= Fenol_connec_Web)/(Reps._Fenol_connec_Web+1),
                                    sum(c(null_Fenol_connec_Web,Fenol_connec_Web) <= Fenol_connec_Web)/(Reps._Fenol_connec_Web+1)))*2
Fenol_connec_Web.obs<-(Fenol_connec_Web - mean(null_Fenol_connec_Web))/sd(null_Fenol_connec_Web)

#Fenol_modhub_Web
Fenol_modhub_Web <- mean(Fenol.dist_web[ModHubs_Web, ModHubs_Web])
null_Fenol_modhub_Web <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Fenol_modhub_Web)){
  Aleatorio <- sample(1:nrow(Fenol.dist_web), size = length(ModHubs_Web))
  null_Fenol_modhub_Web[i] <- mean(Fenol.dist_web[Aleatorio, Aleatorio])
}
quantile(null_Fenol_modhub_Web, probs = c(0.025, 0.975))
rethinking::dens(null_Fenol_modhub_Web, xlim = c(0, 1))
abline(v = Fenol_modhub_Web, col = 2, lwd = 2)

Reps._Fenol_modhub_Web<-10000
P.value.Fenol_modhub_Web <- min(c(sum(c(null_Fenol_modhub_Web,Fenol_modhub_Web) >= Fenol_modhub_Web)/(Reps._Fenol_connec_Web+1),
                                  sum(c(null_Fenol_modhub_Web,Fenol_modhub_Web) <= Fenol_modhub_Web)/(Reps._Fenol_connec_Web+1)))*2

Fenol_modhub_Web.obs<-(Fenol_modhub_Web - mean(null_Fenol_modhub_Web))/sd(null_Fenol_modhub_Web)

#Fenol_periph_web
Fenol_periph_Web <- mean(Fenol.dist_web[Periph_Web, Periph_Web])
null_Fenol_periph_Web <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Fenol_periph_Web)){
  Aleatorio <- sample(1:nrow(Fenol.dist_web), size = length(Periph_Web))
  null_Fenol_periph_Web[i] <- mean(Fenol.dist_web[Aleatorio, Aleatorio])
}
quantile(null_Fenol_periph_Web, probs = c(0.025, 0.975))
rethinking::dens(null_Fenol_periph_Web, xlim = c(0, 1))
abline(v = Fenol_periph_Web, col = 2, lwd = 2)

Reps._Fenol_periph_Web<-10000
P.value.Fenol_periph_Web <- min(c(sum(c(null_Fenol_periph_Web,Fenol_periph_Web) >= Fenol_periph_Web)/(Reps._Fenol_periph_Web+1),
                                  sum(c(null_Fenol_periph_Web,Fenol_periph_Web) <= Fenol_periph_Web)/(Reps._Fenol_periph_Web+1)))*2

####################################Legit##########################################

#Fenol_modhub_Web_legit
Fenol_modhub_Web_legit <- mean(Fenol.dist_legit[ModHubs_Web_legit, ModHubs_Web_legit])
null_Fenol_modhub_Web_legit <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Fenol_modhub_Web_legit)){
  Aleatorio <- sample(1:nrow(Fenol.dist_legit), size = length(ModHubs_Web_legit))
  null_Fenol_modhub_Web_legit[i] <- mean(Fenol.dist_legit[Aleatorio, Aleatorio])
}
quantile(null_Fenol_modhub_Web_legit, probs = c(0.025, 0.975))
rethinking::dens(null_Fenol_modhub_Web_legit, xlim = c(0, 1))
abline(v = Fenol_modhub_Web_legit, col = 2, lwd = 2)

Reps._Fenol_modhub_Web_legit<-10000
P.value.Fenol_modhub_Web_legit <- min(c(sum(c(null_Fenol_modhub_Web_legit,Fenol_modhub_Web_legit) >= Fenol_modhub_Web_legit)/(Reps._Fenol_modhub_Web_legit+1),
                                  sum(c(null_Fenol_modhub_Web_legit,Fenol_modhub_Web_legit) <= Fenol_modhub_Web_legit)/(Reps._Fenol_modhub_Web_legit+1)))*2

#Fenol_connec_Web_legit
Fenol_connec_Web_legit <- mean(Fenol.dist_legit[Connectors_Web_legit, Connectors_Web_legit])
null_Fenol_connec_Web_legit <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Fenol_connec_Web_legit)){
  Aleatorio <- sample(1:nrow(Fenol.dist_legit), size = length(Connectors_Web_legit))
  null_Fenol_connec_Web_legit[i] <- mean(Fenol.dist_legit[Aleatorio, Aleatorio])
}
quantile(null_Fenol_connec_Web_legit, probs = c(0.025, 0.975))
rethinking::dens(null_Fenol_connec_Web_legit, xlim = c(0, 1))
abline(v = Fenol_connec_Web_legit, col = 2, lwd = 2)

Reps._Fenol_connec_Web_legit<-10000
P.value.Fenol_connec_Web_legit <- min(c(sum(c(null_Fenol_connec_Web_legit,Fenol_connec_Web_legit) >= Fenol_connec_Web_legit)/(Reps._Fenol_connec_Web_legit+1),
                                        sum(c(null_Fenol_connec_Web_legit,Fenol_connec_Web_legit) <= Fenol_connec_Web_legit)/(Reps._Fenol_connec_Web_legit+1)))*2

#Fenol_Peri_Web_legit
Fenol_Peri_Web_legit <- mean(Fenol.dist_legit[Periph_Web_legit, Periph_Web_legit])
null_Fenol_Periph_Web_legit <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Fenol_Periph_Web_legit)){
  Aleatorio <- sample(1:nrow(Fenol.dist_legit), size = length(Periph_Web_legit))
  null_Fenol_Periph_Web_legit[i] <- mean(Fenol.dist_legit[Aleatorio, Aleatorio])
}
quantile(null_Fenol_Periph_Web_legit, probs = c(0.025, 0.975))
rethinking::dens(null_Fenol_Periph_Web_legit, xlim = c(0, 1))
abline(v = Fenol_Peri_Web_legit, col = 2, lwd = 2)

Reps._Fenol_Peri_Web_legit<-10000
P.value.Fenol_Peri_Web_legit <- min(c(sum(c(null_Fenol_Periph_Web_legit,Fenol_Peri_Web_legit) >= Fenol_Peri_Web_legit)/(Reps._Fenol_Peri_Web_legit+1),
                                        sum(c(null_Fenol_Periph_Web_legit,Fenol_Peri_Web_legit) <= Fenol_Peri_Web_legit)/(Reps._Fenol_Peri_Web_legit+1)))*2

Fenol_modhub_Web.obs<-(Fenol_Peri_Web_legit - mean(null_Fenol_Periph_Web_legit))/sd(null_Fenol_Periph_Web_legit)


####################################Illegit#######################################
#Connec_Periph_Web_Pill_Fenol
Fenol_Connec_Web_pill <- mean(Fenol.dist_web_pill[Connectors_Web_pill, Connectors_Web_pill])
null_Fenol_Connec_Web_pill <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Fenol_Connec_Web_pill)){
  Aleatorio <- sample(1:nrow(Fenol.dist_web_pill), size = length(Connectors_Web_pill))
  null_Fenol_Connec_Web_pill[i] <- mean(Fenol.dist_web_pill[Aleatorio, Aleatorio])
}
quantile(null_Fenol_Connec_Web_pill, probs = c(0.025, 0.975))
rethinking::dens(null_Fenol_Connec_Web_pill, xlim = c(0, 1))
abline(v = Fenol_Connec_Web_pill, col = 2, lwd = 2)

Reps._Fenol_Connec_Web_pill<-10000
P.value.Fenol_Connec_Web_pill <- min(c(sum(c(null_Fenol_Connec_Web_pill,Fenol_Connec_Web_pill) >= Fenol_Connec_Web_pill)/(Reps._Fenol_Connec_Web_pill+1),
                                      sum(c(null_Fenol_Connec_Web_pill,Fenol_Connec_Web_pill) <= Fenol_Connec_Web_pill)/(Reps._Fenol_Connec_Web_pill+1)))*2

#ModHub_Web_Pill_Fenol
Fenol_Modhub_Web_pill <- mean(Fenol.dist_web_pill[ModHubs_Web_pill, ModHubs_Web_pill])
null_Fenol_ModHubs_Web_pill <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Fenol_ModHubs_Web_pill)){
  Aleatorio <- sample(1:nrow(Fenol.dist_web_pill), size = length(ModHubs_Web_pill))
  null_Fenol_ModHubs_Web_pill[i] <- mean(Fenol.dist_web_pill[Aleatorio, Aleatorio])
}
quantile(null_Fenol_ModHubs_Web_pill, probs = c(0.025, 0.975))
rethinking::dens(null_Fenol_ModHubs_Web_pill, xlim = c(0, 1))
abline(v = Fenol_Modhub_Web_pill, col = 2, lwd = 2)

Reps._Fenol_Modhub_Web_pill<-10000
P.value.Fenol_Modhub_Web_pill <- min(c(sum(c(null_Fenol_ModHubs_Web_pill,Fenol_Modhub_Web_pill) >= Fenol_Modhub_Web_pill)/(Reps._Fenol_Modhub_Web_pill+1),
                                       sum(c(null_Fenol_ModHubs_Web_pill,Fenol_Modhub_Web_pill) <= Fenol_Modhub_Web_pill)/(Reps._Fenol_Modhub_Web_pill+1)))*2

#Periph_Web_Pill_Fenol
Fenol_Peri_Web_pill <- mean(Fenol.dist_web_pill[Periph_Web_pill, Periph_Web_pill])
null_Fenol_Periph_Web_pill <- vector(mode = "numeric", length = 10000)
for(i in 1:length(null_Fenol_Periph_Web_pill)){
  Aleatorio <- sample(1:nrow(Fenol.dist_web_pill), size = length(Periph_Web_pill))
  null_Fenol_Periph_Web_pill[i] <- mean(Fenol.dist_web_pill[Aleatorio, Aleatorio])
}
quantile(null_Fenol_Periph_Web_pill, probs = c(0.025, 0.975))
rethinking::dens(null_Fenol_Periph_Web_pill, xlim = c(0, 1))
abline(v = Fenol_Peri_Web_pill, col = 2, lwd = 2)

Reps._Fenol_Peri_Web_pill<-10000
P.value.Fenol_Peri_Web_pill <- min(c(sum(c(null_Fenol_Periph_Web_pill,Fenol_Peri_Web_pill) >= Fenol_Peri_Web_pill)/(Reps._Fenol_Peri_Web_pill+1),
                                       sum(c(null_Fenol_Periph_Web_pill,Fenol_Peri_Web_pill) <= Fenol_Peri_Web_pill)/(Reps._Fenol_Peri_Web_pill+1)))*2

sum(c(null_Fenol_Periph_Web_pill,Fenol_Peri_Web_pill) >= Fenol_Peri_Web_pill)/(Reps._Fenol_Peri_Web_pill+1)

sum(c(null_Fenol_Periph_Web_pill,Fenol_Peri_Web_pill) <= Fenol_Peri_Web_pill)/(Reps._Fenol_Peri_Web_pill+1)
#Como salvar objetos do R
save(Web, Plants, file = "Objetos.RData")

######################################## Similarity module test Morpho and Phylo #######################################

Mod.Web <- IntraMod.dist(Modulos, dist = Plant.dist)
Mod.Web.Legit<- IntraMod.dist(Modulos_legit, dist = Plant.dist_legit)
Mod_Web_pill<- IntraMod.dist(Modulos_pill, dist = Plant.dist_pill)


Mod.Filo.Web <- IntraMod.dist(Modulos, dist = Filo.dist)
Mod.Filo.Web.Legit<- IntraMod.dist(Modulos_legit, dist = Filo.dist.legit)
Mod.Filo.Web.Pill<- IntraMod.dist(Modulos_pill, dist = Filo.dist.pill)


##################################### Phenology ####################
library(spaa)

################################Overall Phenology#################
Data.raw_unique <- Data.raw[!duplicated(Data.raw[, c("Plant_species", "Date")]),]
Data.raw_unique$Month <- unlist(lapply(strsplit(as.character(Data.raw_unique$Date), split = " "),
                                       FUN = function(x) x[2]))
occurrence_table <- table(Data.raw_unique$Month, Data.raw_unique$Plant_species)

Fenol.dist_web <- as.matrix(niche.overlap(occurrence_table, method = c("czech")))
Mod.Fenol.Web <- IntraMod.dist(Modulos, dist = Fenol.dist_web)
dens(Mod.Fenol.Web$Null.mean, xlim = c(0, 1))
abline(v = Mod.Fenol.Web$OBS.mean, col = 2, lwd = 2)
Data.raw$Date

############################################# Legit Pheno. ####################
Fenol.dist_legit <- Fenol.dist_web[rownames(Fenol.dist_web) %in% rownames(Web.legit),
                                   colnames(Fenol.dist_web) %in% rownames(Web.legit)]
Mod.Fenol.Web_legit <- IntraMod.dist(Modulos_legit, dist = Fenol.dist_legit)

dens(Mod.Fenol.Web_legit$Null, xlim = c(0, 1))
abline(v = Mod.Fenol.Web_legit$OBS, col = 2, lwd = 2)
Data.legit_unique$Date

###############################################Illegit Pheno. ######################
Fenol.dist_web_pill <- Fenol.dist_web[rownames(Fenol.dist_web) %in% rownames(Web.pill),
                                      colnames(Fenol.dist_web) %in% rownames(Web.pill)]
Mod.Fenol.Web_pill <- IntraMod.dist(Modulos_pill, dist = Fenol.dist_web_pill)
dens(Mod.Fenol.Web_pill$Null.mean, xlim = c(0, 1))
abline(v = Mod.Fenol.Web_pill$OBS.mean, col = 2, lwd = 2)
Data.pill_unique$Date

################## Trait mean in each module ##############################
## Overall ##
soma_interacoes <- rowSums(Web)
Plants$Int <- soma_interacoes[match(rownames(Plants), rownames(Web))]


CWMs <- data.frame(N = unlist(lapply(Modulos_SPs, length)),
                   row.names = 1:length(Modulos_SPs))
for(i in 1:length(Modulos_SPs)){
  CWMs$comprimento[i] <- mean(Plants$comprimento[rownames(Plants) %in% Modulos_SPs[[i]]])
  CWMs$recurso[i] <- mean((as.numeric(Plants$Recurso_fornecido)-1)[rownames(Plants) %in% Modulos_SPs[[i]]])
}

for(i in 1:length(Modulos_SPs)){
  Plants.mod <- Plants[rownames(Plants) %in% Modulos_SPs[[i]],]
  CWMs$comprimento.pond[i] <- weighted.mean(Plants.mod$comprimento, w = Plants.mod$Int)
  CWMs$recurso.pond[i] <- weighted.mean(as.numeric(Plants.mod$Recurso_fornecido) - 1, w = Plants.mod$Int)
}

write.xlsx(CWMs, file = "Mean_morpho_dif.xlsx")
########################### Modules dif. model #####################################

Plants$Modulo <- NA

for (i in seq_along(Modulos_SPs)) {

  especies_modulo <- Modulos_SPs[[i]]

  Plants$Modulo[rownames(Plants) %in% especies_modulo] <- i
}

head(Plants)
sum(Plants$Int)
Plants$Modulo <- as.factor(Plants$Modulo)

####### Dist. Gama ###########

Plants$comprimento[Plants$comprimento == 0] <- 0.0001

model_gama <- glm(comprimento ~ Modulo, data = Plants, weights = Int, family = Gamma)

summary(model_gama)
summary(aov(model_gama))
emmeans_gamma<-emmeans(model_gama, ~ Modulo, type = "response")
summary(emmeans_gamma)
summary(model_gama)$coefficients[, "Std. Error"]

######## species role phenology model #####################

model_nb_role <- glm.nb(flowering_duration ~ role, data = duration_table_filtered)
summary(model_nb_role)
emmeans_nb<-emmeans(model_nb_role, pairwise ~ role, type = "response")
summary(emmeans_nb)

####################### Figures ####################
############################# Network Plots #############################

png("rede_modulos.png", width = 2400, height = 1800, res = 300)

par(cex = 0.7)
plotModuleWeb(Modulos, displayAlabels = T, displayBlabels = T, weighted = T, labsize = 0.5)
dev.off()

svg("rede_modulos_legit.svg")
plotModuleWeb(Modulos_legit, displayAlabels = F, displayBlabels = F, weighted = F)
dev.off()

svg("rede_modulos_pill.svg")
plotModuleWeb(Modulos_pill, displayAlabels = F, displayBlabels = F, weighted = F)
dev.off()

# Criar paletas de cores com trans  parência ajustada
color_palette <- colorspace::qualitative_hcl(10, palette = "Set2")
color_palette_alpha <- adjustcolor(color_palette, alpha.f = 1)

color_palette.legit <- colorspace::qualitative_hcl(8, palette = "Set2")
color_palette_alpha.legit <- adjustcolor(color_palette.legit, alpha.f = 1)

color_palette.pill <- colorspace::qualitative_hcl(9, "Set2")
color_palette_alpha.pill <- adjustcolor(color_palette.pill, alpha.f = 1)

# Ajustar cores de interação para módulos
Col.Int <- t(outer(Plantas.Mod$modulo, Insetos.Mod$modulo,
                   FUN = function(X, Y) {
                     ifelse(X == Y,
                            color_palette_alpha[Plantas.Mod$modulo],
                            "#99999996")
                   }))

Col.Int.legit <- t(outer(Plantas.Mod.legit$modulo, Insetos.Mod.legit$modulo,
                         FUN = function(X, Y) {
                           ifelse(X == Y,
                                  color_palette_alpha.legit[Plantas.Mod.legit$modulo],
                                  "#99999996")
                         }))

Col.Int.pill <- t(outer(Plantas.Mod.pill$modulo, Insetos.Mod.pill$modulo,
                        FUN = function(X, Y) {
                          ifelse(X == Y,
                                 color_palette_alpha.pill[Plantas.Mod.pill$modulo],
                                 "#99999996")
                        }))

# Gerar plots SVG e PNG
svg("redes_plotadas.svg", width = 24, height = 36)
png("redes_plotadas.png", width = 2000, height = 3000, res = 300)
par(mfrow = c(3, 1), mar = c(0, 4, 0, 2))

# Plot 1: Rede geral
plotweb(Modulos@moduleWeb, method = "normal",
        col.high = color_palette[Insetos.Mod$modulo],
        col.low = color_palette[Plantas.Mod$modulo],
        col.interaction = c(Col.Int),
        high.lablength = 0,
        low.lablength = 0,
        bor.col.interaction = rgb(0,0,0,0),
        bor.col.high = rgb(0,0,0,0),
        bor.col.low = rgb(0,0,0,0),
        low.spacing = 0.001, high.spacing = 0.001)
mtext("A", side = 3, line = 0.5, adj = 0, cex = 1, font = 1.5)
mtext("Animal species", side = 3, line = -1.5, adj = 0, cex = 0.9)
mtext("Plant species", side = 1, line = -1.5, adj = 0, cex = 0.9)

# Plot 2: Rede legítima
plotweb(Modulos_legit@moduleWeb, method = "normal",
        col.high = color_palette.legit[Insetos.Mod.legit$modulo],
        col.low = color_palette.legit[Plantas.Mod.legit$modulo],
        col.interaction = c(Col.Int.legit),
        high.lablength = 0,
        low.lablength = 0,
        bor.col.interaction = rgb(0,0,0,0),
        bor.col.high = rgb(0,0,0,0),
        bor.col.low = rgb(0,0,0,0),
        low.spacing = 0.001, high.spacing = 0.001)
mtext("B", side = 3, line = 0.5, adj = 0, cex = 1, font = 1.5)
mtext("Animal species", side = 3, line = -1.5, adj = 0, cex = 0.9)
mtext("Plant species", side = 1, line = -1.5, adj = 0, cex = 0.9)

# Plot 3: Rede pill
plotweb(Modulos_pill@moduleWeb, method = "normal",
        col.high = color_palette.pill[Insetos.Mod.pill$modulo],
        col.low = color_palette.pill[Plantas.Mod.pill$modulo],
        col.interaction = c(Col.Int.pill),
        high.lablength = 0,
        low.lablength = 0,
        bor.col.interaction = rgb(0,0,0,0),
        bor.col.high = rgb(0,0,0,0),
        bor.col.low = rgb(0,0,0,0),
        low.spacing = 0.001, high.spacing = 0.001)
mtext("C", side = 3, line = 0.5, adj = 0, cex = 1, font = 1.5)
mtext("Animal species", side = 3, line = -1.5, adj = 0, cex = 0.9)
mtext("Plant species", side = 1, line = -1.5, adj = 0, cex = 0.9)

dev.off()

################################### Joy plot ##############################

occurrence_table.graphic <- as.data.frame(occurrence_table)
colnames(occurrence_table.graphic) <- c("Mes", "Especie", "Frequency")

combined_data <- merge(occurrence_table.graphic, mod_df, by.x = "Especie", by.y = "Especie")

combined_data$Mes <- factor(combined_data$Mes, levels = month.name)

combined_data <- combined_data %>%
  mutate(Especie = factor(Especie, levels = unique(Especie[order(Modulo)])))
ordered_species <- unique(combined_data$Especie)

combined_data <- combined_data %>%
  complete(Mes = factor(month.name, levels = month.name), Especie, fill = list(Frequency = 0))

combined_data <- combined_data %>%
  complete(Mes = factor(month.name, levels = month.name), Especie, fill = list(Frequency = 0))
combined_data %>%
  group_by(Mes) %>%
  summarize(Soma_Frequency = sum(Frequency))

ggplot(combined_data, aes(x = Mes, y = Especie, height = Frequency, fill = factor(Modulo), group = Especie)) +
  geom_density_ridges(
    stat = "identity",
    alpha = 0.5) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Month", y = "Plant Species", fill = "Module") +
  theme_ridges(font_size = 20, grid = F) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 20),
    axis.text.y = element_blank()
  )

# Salvar o gráfico
ggsave("joyplot_fixed.svg", width = 12, height = 12)
ggsave("joyplot_fixed.png", width = 12, height = 12)
#### o problem está na hora de ver as frequncias, não está pegando correto do combined_data, exlcuindo algumas olhe para frequncia da ultima que é althernanthera tenela, mes dezembro e julho está vazio



especies <- c()
modulos <- c()

for (i in seq_along(Modulos_SPs)) {
  especies <- c(especies, Modulos_SPs[[i]])
  modulos <- c(modulos, rep(i, length(Modulos_SPs[[i]])))
}

mod_df <- data.frame(Especie = especies, Modulo = modulos)

colnames(occurrence_table)
occurrence_table.graphic <- as.data.frame(occurrence_table)
colnames(occurrence_table.graphic) <- c("Mes", "Especie", "Frequency")

combined_data <- merge(occurrence_table.graphic, mod_df, by.x = "Especie", by.y = "Especie")

combined_data$Mes <- factor(combined_data$Mes, levels = month.name)

combined_data <- combined_data %>%
  group_by(Especie) %>%
  mutate(Total_Frequency = sum(Frequency)) %>%
  ungroup() %>%
  mutate(Especie = factor(Especie, levels = names(sort(tapply(Total_Frequency, Especie, sum), decreasing = TRUE))))

combined_data$Modulo <- factor(combined_data$Modulo)


n_modulos <- length(unique(combined_data$Modulo))
paleta_modulos <- rainbow(n_modulos)


combined_data <- combined_data %>%
  mutate(ModuloCor = paleta_modulos[Modulo])


generate_gradient <- function(color) {
  scales::col_numeric(palette = c("white", color), domain = c(0, max(combined_data$Frequency)))
}


n_modulos <- length(unique(combined_data$Modulo))
paleta_modulos <- rainbow(n_modulos)

combined_data <- combined_data %>%
  mutate(ModuloCor = paleta_modulos[Modulo])

generate_gradient <- function(color) {
  scale_fill_gradient(low = "white", high = color)
}

plots <- list()

for(modulo in unique(combined_data$Modulo)) {
  data_modulo <- combined_data %>% filter(Modulo == modulo)

  p <- ggplot(data_modulo, aes(x = Mes, y = Especie, fill = Frequency)) +
    geom_tile(color = "white", size = 0.1) +
    generate_gradient(unique(data_modulo$ModuloCor)) +
    labs(title = paste("Frequency Distribution of Species by Month - Module", modulo),
         x = "Month",
         y = "Species",
         fill = "Frequency") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 10, face = "bold"))

  plots[[modulo]] <- p  # Armazena o gráfico na lista
}

plots_set1 <- plots[1:4]
plots_set2 <- plots[5:8]
plots_set3 <- plots[9:10]

do.call(grid.arrange, c(plots_set1, ncol = 2))

do.call(grid.arrange, c(plots_set2, ncol = 2))

do.call(grid.arrange, c(plots_set3, ncol = 2))

png("plots_set1.png", width = 1500, height = 1000)
grid.arrange(grobs = plots_set1, ncol = 2)
dev.off()

png("plots_set2.png", width = 1500, height = 1000)
grid.arrange(grobs = plots_set2, ncol = 2)
dev.off()

png("plots_set3.png", width = 1500, height = 1000)
grid.arrange(grobs = plots_set3, ncol = 2)
dev.off()

######################################### Box plot #########################
library(dplyr)
library(ggbeeswarm)
library(forcats)

emmeans_df <- as.data.frame(emmeans_gamma)
emmeans_df <- emmeans_df %>%
  mutate(
    has_inf = is.infinite(upper.CL),
    upper.CL = ifelse(has_inf, NA, upper.CL))


Plants.mod.graphic <- Plants %>%
  mutate(comprimento_vis = ifelse(abs(comprimento - 0.0001) < 1e-6, 0.1, comprimento))

ggplot() +
  geom_beeswarm(data = Plants.mod.graphic, aes(x = as.factor(Modulo), y = comprimento_vis, color = as.factor(Modulo)),
                alpha = 0.7, size = 2) +

  geom_point(data = emmeans_df, aes(x = as.factor(Modulo), y = response),
             color = "black", size = 2, shape = 17) +

  geom_errorbar(data = emmeans_df %>% filter(!has_inf),
                aes(x = as.factor(Modulo),
                    ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.2, color = "black") +

  geom_text(data = emmeans_df %>% filter(has_inf),
            aes(x = as.factor(Modulo), y = response + 5, label = "*"),
            color = "black", size = 6, vjust = -0.5) +

  scale_y_log10(
    breaks = c(0.1, 1, 10, 100),
    labels = c("<0.1", "1", "10", "100")
  ) +

  xlab("Modules") +
  ylab("Corolla Depth (mm)") +

  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15)
  )
ggsave("geom_point_with_emmeans_and_CI.png", width = 6, height = 4, dpi = 300)
############################## Fig. species role #####################

png("multiplot.png", width = 5, height = 14, units = "in", res = 600)

par(mfrow = c(3, 1))

for (i in 1:3) {
  plot(df_web_clean$c, df_web_clean$z, col = df_web_clean$role, pch = 16,
       xlab = "Among-module connectivity, C", ylab = "Within-module degree, Z",
       cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, xlim = c(0, 1.1))

  # Linhas críticas em todos os gráficos
  abline(h = Z_critical_l, lty = 3)
  abline(v = C_critical_l, lty = 3)

  if (i == 1) {
      # Texto azul e verde, já com as alturas ajustadas
      text(0.02, Z_critical_l - 0.2, "Peripheral", col = 4, font = 2, adj = 0, cex = 1.6)
      text(0.02, max(df_web_clean$z) + 0.05, "Module Hub", col = "green3", font = 2, adj = 0, cex = 1.6)
      text(C_critical_l + 0.04, Z_critical_l - 0.2, "Connector", col = "black", font = 2, adj = 0, cex = 1.6) #+0.04 (direita e esquerd) - 0.2 cima e baixo
      text(C_critical_l + 0.03, max(df_web_clean$z) + 0.02, "Network Hub", col = 2, font = 2, adj = 0, cex = 1.6) #cex tamanho coo sempre

  }
}

dev.off()


########################### Matriz ##############################
png("rede_modulos.png", width = 2000, height = 2000, res = 300)
plotModuleWeb(Modulos, displayAlabels = TRUE, displayBlabels = TRUE, weighted = TRUE, labsize = 0.5)
dev.off()

# Gráfico 2: Modulos_legit
png("rede_modulos_legit.png", width = 2000, height = 2000, res = 300)
plotModuleWeb(Modulos_legit, displayAlabels = TRUE, displayBlabels = TRUE, weighted = TRUE, labsize = 0.5)
dev.off()

# Gráfico 3: Modulos_pill
png("rede_modulos_pill.png", width = 2000, height = 2000, res = 300)
plotModuleWeb(Modulos_pill, displayAlabels = TRUE, displayBlabels = TRUE, weighted = TRUE, labsize = 0.5)
dev.off()

################### Supplementary Material #####################

#Polliantors################
library(tibble)
library(dplyr)
library(purrr)

modulos_pol_WEB <- map2_dfr(
  listModuleInformation(Modulos)[[2]],
  seq_along(listModuleInformation(Modulos)[[2]]),
  ~ tibble(polinizador = .x[[2]], modulo_main = .y)
)

modulos_pol_WEB_legit <- map2_dfr(
  listModuleInformation(Modulos_legit)[[2]],
  seq_along(listModuleInformation(Modulos_legit)[[2]]),
  ~ tibble(polinizador = .x[[2]], modulo_legit = .y)
)

modulos_pol_WEB_pill <- map2_dfr(
  listModuleInformation(Modulos_pill)[[2]],
  seq_along(listModuleInformation(Modulos_pill)[[2]]),
  ~ tibble(polinizador = .x[[2]], modulo_pill = .y)
)
padronizar_polinizador <- function(df) {
  df %>%
    mutate(polinizador = trimws(tolower(polinizador)))
}

library(dplyr)
Supplementary_Material_7 <- read.csv2("Supplementary Material 7.csv")
Supplementary_Material_7 <- Supplementary_Material_7 %>%
  mutate(pollinator_species = trimws(tolower(pollinator_species))) %>%
  left_join(modulos_pol_WEB      %>% mutate(polinizador = trimws(tolower(polinizador))),
            by = c("pollinator_species" = "polinizador")) %>%
  left_join(modulos_pol_WEB_legit %>% mutate(polinizador = trimws(tolower(polinizador))),
            by = c("pollinator_species" = "polinizador")) %>%
  left_join(modulos_pol_WEB_pill  %>% mutate(polinizador = trimws(tolower(polinizador))),
            by = c("pollinator_species" = "polinizador"))


write.csv2(as.matrix(Supplementary_Material_7),
           file = "Supplementary_Material_7.csv",
           row.names = TRUE)

dados_web <- list(
  "web" = Web,
  "web.legit" = Web.legit,
  "web.pill" = Web.pill
)
write.xlsx(dados_web, file = "dados_web.xlsx")
dados_CZ <- list(
  "df_web_clean" = df_web_clean,
  "df_web_clean_legit" = df_web_clean_legit,
  "df_web_clean_pill" = df_web_clean_pill
)

# Salvar em um arquivo Excel
library(openxlsx)
write.xlsx(dados_CZ, file = "CZ_Values.xlsx")

fenologia_dados <- list(
  "Fenol_Dist_Legit" = Fenol.dist_legit,
  "Fenol_Dist_Web" = Fenol.dist_web,
  "Fenol_Dist_Web_Pill" = Fenol.dist_web_pill
)

# Salvar como Excel
library(openxlsx)
write.xlsx(fenologia_dados, file = "Fenol_Dist.xlsx")


filogenia_dados <- list(
  "Filo_dist_web" = Filo.dist,
  "Filo.dist.legit" = Filo.dist.legit,
  "Filo.dist.pill" = Filo.dist.pill
)

# Salvar como Excel
library(openxlsx)
write.xlsx(filogenia_dados, file = "phylo_Dist.xlsx")


morpho_dados <- list(
  "Morpho.dist" = Plant.dist,
  "Morpho.dist_legit" = Plant.dist_legit,
  "Morpho.dist_pill" = Plant.dist_pill
)

# Salvar como Excel
library(openxlsx)
write.xlsx(morpho_dados, file = "morpho_dist.xlsx")


legitimate_interactions <- subset(Data, Interaction == "1")
illegitimate_interactions <- subset(Data, Interaction == "2")
legitimate_matrix <- table(legitimate_interactions$Plant_species, legitimate_interactions$Insect_species)
illegitimate_matrix <- table(illegitimate_interactions$Plant_species, illegitimate_interactions$Insect_species)
write.csv2(as.matrix(legitimate_matrix),
           file = "Supplementary Material 8.csv",
           row.names = TRUE)

# Salvando a matriz ilegítima
write.csv2(as.matrix(illegitimate_matrix),
           file = "Supplementary Material 9.csv",
           row.names = TRUE)



library(dplyr)
Data_filtered <- Data %>%
  select(Plant_species, Month)

Data_with_family <- Data_filtered %>%
  left_join(Regional, by = c("Plant_species" = "Especie"))

Plants_with_species <- Plants
Plants_with_species$Especie <- rownames(Plants)
final_data <- Data_with_family %>%
  left_join(Plants_with_species, by = c("Plant_species" = "Especie"))
final_data <- final_data %>% select(-Genero, -Int)
final_data_unique <- final_data %>% distinct()
final_data_combined <- final_data %>%
  group_by(Plant_species, Familia, Recurso_fornecido, comprimento, Modulo) %>%
  summarise(Months = paste(unique(Month), collapse = ", ")) %>%
  ungroup()
colnames(final_data_combined)[colnames(final_data_combined) == "comprimento (mm, x̄)"] <- "Depth (mm, x̄)"
colnames(final_data_combined)[colnames(final_data_combined) == "Familia"] <- "Family"
colnames(final_data_combined)[colnames(final_data_combined) == "Recurso_fornecido"] <- "Floral resource"
colnames(final_data_combined)[colnames(final_data_combined) == "Modulo"] <- "Module"

write.csv2(final_data_combined, "Supplementary Material 6.csv", row.names = FALSE)

# Criando cópias dos dataframes iniciais para preservar os originais
df_web_clean_copy <- df_web_clean
df_web_clean_legit_copy <- df_web_clean_legit
df_web_clean_pill_copy <- df_web_clean_pill

# Adicionando a coluna 'Species' a partir dos rownames nas cópias
df_web_clean_copy$Species <- rownames(df_web_clean)
df_web_clean_legit_copy$Species <- rownames(df_web_clean_legit)
df_web_clean_pill_copy$Species <- rownames(df_web_clean_pill)

# Realizando os merges com base na coluna 'Species'
merged_df <- merge(df_web_clean_copy, df_web_clean_legit_copy, by = "Species", all = TRUE)
merged_df <- merge(merged_df, df_web_clean_pill_copy, by = "Species", all = TRUE)

# Renomeando colunas para indicar a origem dos valores
colnames(merged_df) <- c("Species",
                         "c_web", "z_web", "role_web", "module_web",
                         "c_legit", "z_legit", "role_legit",
                         "c_pill", "z_pill", "role_pill", "module_pill")

# Removendo as colunas 'module_web' e 'module_pill' do dataframe final
merged_df$module_web <- NULL
merged_df$module_pill <- NULL

# Salvando o dataframe resultante como CSV
write.csv2(merged_df, file = "Supplementary Material 5.csv", row.names = FALSE)

######## Phylogenetic tree plot #######
png("arvore_alta.png", width = 2000, height = 5000, res = 300)
plot.phylo(Filo.Regional, use.edge.length = TRUE, cex = 0.8)
dev.off()

######################### phenology ####################

df_web_clean_dispersion <- tibble::rownames_to_column(df_web_clean, "Espécie")
colnames(df_web_clean_dispersion)[colnames(df_web_clean_dispersion) == "Espécie"] <- "Especie"
df_web_clean_dispersion <- df_web_clean_dispersion[, c("Especie", "role")]
df_web_clean_dispersion <- merge(occurrence_table.graphic, df_roles, by = "Especie")
df_web_clean_dispersion$Mes <- factor(df_web_clean_dispersion$Mes,
                                      levels = c("January", "February", "March", "April", "May", "June",
                                                 "July", "August", "September", "October", "November", "December"),
                                      ordered = TRUE)
df_web_clean_dispersion$Mes_num <- as.numeric(df_web_clean_dispersion$Mes)
library(dplyr)
dispersion_table <- df_web_clean_dispersion %>%
  filter(Frequency > 0) %>%
  group_by(Especie, role) %>%
  summarise(feno_dispersion = var(Mes), .groups = "drop")

duration_table <- df_web_clean_dispersion %>%
  filter(Frequency > 0) %>%
  group_by(Especie, role) %>%
  summarise(flowering_duration = n_distinct(Mes_num), .groups = "drop")

duration_table_filtered <- duration_table %>%
  group_by(role) %>%
  filter(n() > 1)

duration_table_filtered <- duration_table_filtered %>%
  left_join(select(Plants_with_species, Especie, Int), by = "Especie")

