#############################################################                                                  #
#  Statistiques sur R                                       #        
#############################################################


data<- read.csv("data2/data_temperature.csv", row.names=1)
fix(data)



temp<-data$temperature
mean(temp)
hist(temp)

############### test de conformité
# Test de normalité
shapiro.test(temp) 
#l'hypothèse 0 : les données sont issues d'une population 
#normalement distribuée

t.test(temp, mu = 24.3) 
#test de student, comparaison de moyenne

# on peut modifier IC
t.test(temp, mu = 24.3,conf.level=0.9)

# test unilatéral
t.test(temp, mu = 24.3, alternative="less")
t.test(temp, mu = 24.3, alternative="greater")



############### égalité des moyennes non appariées
daur<- read.csv("data2/daurades.csv", row.names=1)
fix(daur)

# Test de normalité
shapiro.test(daur[,1])
shapiro.test(daur[,2])
# Conclusion test de Shapiro: p-value = 0.08785 
#et 0.53.66>0.05 ==> normalité des données

# Test d'homoscédasticité: test F de Fisher Snedecor 
#(Homogénéité des variances lorsque la distribution est normale)
var.test(daur[,1], daur[,2])

# Conclusion sur l'homoscédasticité: pvalue>0.05 ==> garder H0 : homoscédasticité des données
#(H1 étant : true ratio of variances is not equal to 1 i.e. pas d'homoscédasticité)
  
# Test de Student sur échantillons non appariés (vrai test t)
t.test(daur[,1],daur[,2], var.equal = TRUE, paired = FALSE)
#H0 du t.test les deux moyennes sont égales
  


############### égalité des moyennes appariés
  
  cd<- read.csv("data2/cadmium.csv", row.names=1) #concentration de Cm dans le muscle et les gonades par individu
  fix(cd)
  
  # Test de normalité sur les différences
  #Pour le test de student apparié, il faut vérifier 
  #la différence de (= x-y) suit une loi normale
  
  shapiro.test(cd[,1]-cd[,2])
  
  # Conclusion test de Shapiro: p-value >0.05 ==> normalité des différences
  
  # Test de Student sur échantillons appariés 
  t.test(cd[,1], cd[,2], paired = TRUE)
  
  #Conclusion: pvalue>0.05 ==> garder H0, H0 étant l'égalité des moyennes ==> pas de différence entre les
  #moyennes ==> pas de différence significative entre les concentrations en Cd dans les 2 organes
  #(H1 étant true difference in means is not equal to 0 ==> pas d'égalité des moyennes)
    
    # En cas de non homoscédasticité : test de Wilcoxon:
    #wilcox.test(cd[,1],cd[,2], paired = TRUE)
    
    
    ############### corrélation
    
    crev<- read.csv("crevette.csv", row.names=1)
    
    # Test de normalité
    shapiro.test(crev$Taille)
    shapiro.test(crev$Antenne)
    
    taille<-crev$Taille
    antenne<-crev$Antenne
    
    cor.test(antenne, taille, method = "pearson")
    
    model<-lm(taille~antenne)
    model
    
    plot(taille~antenne)
    
    abline(lm(taille ~ antenne),col="red")
    abline(model, col="red")
    summary(model)
    
  boxplot(crev)
    
    
  
  
  ############### Analyses exploratoires ###############
    spe <- read.csv("data2/DoubsSpe.csv", row.names=1)
    env <- read.csv("data2/DoubsEnv.csv", row.names=1)
    spa <- read.csv("data2/DoubsSpa.csv", row.names=1)
    
    
    # Fonctions de base
    ----------------------
    # Dimension de la matrice
    dim(spe)
    # Label des colonnes
    #colnames(spe)
    #names(spe)
    # Label des lignes
    #rownames(spe)
    # Ouvrir le fichier
    fix(spe)
    
    
    # Statistiques descriptives des colonnes
    summary(env)
    
    
    # Exploration des données
    
    # Plot des stations d'échantillonnage
    plot(spa, xlab="Longitude (km)", ylab="Latitude (km)", 
         main="Localisation géographique des sites")
    
    # Ajouter une ligne bleue qui relie les stations
    lines(spa, col="light blue")
    
    # Ajouter le né des stations
    text(spa, row.names(spa), cex=0.7, col="steelblue", pos=1)
    
    # Ajouter "Amont" "Aval"
    text(50, 10, "Amont", cex=1, col="black")
    text(30, 120, "Aval", cex=1, col="black")
    
    # Profil des variables environnementales
    fix(env)
    boxplot(env$alt, main="Altitude")
    boxplot(env$oxy, main="Oxygène")
    boxplot(env$nit, main="Nitrate")
    
    boxplot(env$nit, env$oxy, env$alt, main="Paramètres abiotiques")
    
    abio_data<-data.frame(env$alt,env$oxy,env$nit)
    boxplot(scale(abio_data))
    
    
    # Plots de cartes
   
      # Plot des abondances par espèce
      
      par(mfrow=c(2,2))
    
    plot(spa, cex=spe$TRU, col="brown", main="Truite", 
         xlab="Longitude (km)", ylab="Latitude (km)")
    lines(spa, col="light blue")
    
    plot(spa, cex=spe$Gou, col="brown", main="Goujon", 
         xlab="Longitude (km)", ylab="Latitude (km)")
    lines(spa, col="light blue")
    
    plot(spa, cex=spe$BAR, col="brown", main="Barbeau", 
         xlab="Longitude  km)", ylab="Latitude (km)")
    lines(spa, col="light blue")
    
    plot(spa, cex=spe$ROT, col="brown", main="Rotengle", 
         xlab="Longitude (km)", ylab="Latitude (km)")
    lines(spa, col="light blue")
    
    
    # Plot des variables environnementales
    plot(spa, main="Altitude", pch=21, col="white", 
         bg="red", cex=5*env$alt/max(env$alt), xlab="Longitude (km)",
         ylab="Latitude (km)")
    lines(spa, col="light blue")
    
    
    # Calcul d'indices écologiques
  
      library(vegan)
    
    # RS
    RS <- apply(spe > 0, 1, sum) #spe > 0 pour compter les occurences
    
    # Indice d'entropie de Shannon
    H <- diversity(spe,index="shannon") 
    H
    
    # Indice de diversité de Simpson
    N<-diversity(spe,index="simpson")
    N
    
    # Indice d'équitabilité de Pielou
    J<-H/log(RS)
    J
    
    # Plot des indices 
    plot(spa, main="carte de RS", pch=21, col="white", bg="brown", 
         cex=5*RS/max(RS), xlab="Longitude", ylab="Latitude")
    lines(spa, col="light blue")
    
    plot(spa, main="indice de Shannon", pch=21, col="white", 
         bg="brown", cex=H, xlab="Longitude", ylab="Latitude")
    lines(spa, col="light blue")
    
    
    # ACP
    ---------------------
      # avec FactoMine R
      
    library (FactoMineR)
    acp.spe <- PCA(spe, scale=TRUE)
    acp.spe
    acp.spe$eig
    #Hierarchical Clustering on Principle Components (HCPC)
    HCPC(spe, nb.clust=-1)
    
 ###########################end######################""   
    
    # Classification hiérarchique
    library(ade4)

    spe.norm1<- decostand(spe, "standardize") #standardiser les données
    spe.norm2<- vegdist(spe.norm1, "euc") #calculer les distances
    
    dendro<-hclust(spe.norm2, method="average")
    plot(dendro, xlab="stations", ylab="similarité", main="dendrogramme")
    
    rect.hclust(dendro, k=5)
    
    
    
    