###MULLER's paper
##February-March 2022

##Hypothesis
#Phylogeny will be the most important factor in explaining variability on the mineral nutrient concentration of plant organs.
#Although non-phylogenetic relationships among species, such as affinity for gypsum soils, will also play a relevant role in certain nutrient concentration, but not soil nutrient concentration.
#In particular, we expect wide gypsophiles to be significant related to S, Ca, and Mg concentrations in all organs, especially in leaves, while gypsovags and narrow gypsophiles will not have a significant relationship with any mineral nutrient concentration of plant organs. 

##How to test this hypothesis?
#using MCMCglmm
#First....analyzing the effect of phylogeny, non-phylogentic relationships and intraspecific variation in nutrient concentration
#Second...including affinity to gypsum soils as fixed effect, as well as, with soil concentration

setwd("D:/Clare Muller/article")
library(readxl)
library(ggplot2)
library(Rmisc)
library(scales)
library(car)
library(ggrepel)
library(devtools)
library(grid)
library(reshape2)
library(tidyr)
library(lmtest)
library(Hmisc)
library(ggrepel)
library(ape)
library(MCMCglmm)
data<- read_excel("data_all_subset_11022022.xlsx", sheet="data_1st")


###LEAVES

#tree of target species
tree <- read.tree("Clare_tree36_tree.tre", keep.multi = FALSE)
plot(tree)


##MCMCglmmm models
#to prepare phylogeny as random factor
phylo.inv<-inverseA(tree,nodes="TIPS",scale=TRUE)
prior2<-list(G=list(G1=list(V=1,nu=0.02),G2=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))


##effect of phylogeny and other relationships
data2<-data[which(data$fractfix=="leaves"), ]
data2$species<-as.factor(data2$species)

#select nutrient
test.data<-data.frame(variable=data2$Ca, phylo=data2$phylo, species=data2$species)

#model without fixed effects
model1<-MCMCglmm(variable~1, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model1)

lambda <- model1$VCV[,'phylo']/
  (model1$VCV[,'phylo']+model1$VCV[,'species']+
     model1$VCV[,'units'])

mean(lambda)
HPDinterval(lambda)

gamma<-model1$VCV[,'species']/
  (model1$VCV[,'phylo']+model1$VCV[,'species']+
     model1$VCV[,'units'])

mean(gamma)
HPDinterval(gamma)

rho<-model1$VCV[,'units']/
  (model1$VCV[,'phylo']+model1$VCV[,'species']+
     model1$VCV[,'units'])

mean(rho)
HPDinterval(rho)

#model per each gypsum affinity + soil
test.data<-data.frame(variable=data2$N, country=data2$country,factor=data2$statfix, factor2=data2$statfix2, soil=data2$soilN,phylo=data2$phylo, species=data2$species)
model2<-MCMCglmm(variable~-1+factor+soil, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model2)

#model between gypsum affinity + soil
model3<-MCMCglmm(variable~factor+soil, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model3)

model4<-MCMCglmm(variable~factor2+soil, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model4)

#model per each region
#wide
data3<-test.data[which(test.data$factor=="wide"), ]
model5<-MCMCglmm(variable~-1+country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=data3,nitt=150000,burnin=1000,thin=50) 
summary(model5)

model6<-MCMCglmm(variable~country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model6)

#narrow
data3<-test.data[which(test.data$factor=="narrow"), ]
model5<-MCMCglmm(variable~-1+country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=data3,nitt=150000,burnin=1000,thin=50) 
summary(model5)

model6<-MCMCglmm(variable~country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model6)

#gypsovags
data3<-test.data[which(test.data$factor=="gypsovag"), ]
model5<-MCMCglmm(variable~-1+country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=data3,nitt=150000,burnin=1000,thin=50) 
summary(model5)

model6<-MCMCglmm(variable~country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model6)


###STEMS

#tree of target species
tree <- read.tree("Clare_tree36_tree.tre", keep.multi = FALSE)
plot(tree)

phylo.inv<-inverseA(tree,nodes="TIPS",scale=TRUE)
prior2<-list(G=list(G1=list(V=1,nu=0.02),G2=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))

data2<-data[which(data$fractfix=="stems"), ]
data2$species<-as.factor(data2$species)

test.data<-data.frame(variable=data2$N, phylo=data2$phylo, species=data2$species)

model1<-MCMCglmm(variable~1, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model1)

lambda <- model1$VCV[,'phylo']/
  (model1$VCV[,'phylo']+model1$VCV[,'species']+
     model1$VCV[,'units'])

mean(lambda)
HPDinterval(lambda)

gamma<-model1$VCV[,'species']/
  (model1$VCV[,'phylo']+model1$VCV[,'species']+
     model1$VCV[,'units'])

mean(gamma)
HPDinterval(gamma)

rho<-model1$VCV[,'units']/
  (model1$VCV[,'phylo']+model1$VCV[,'species']+
     model1$VCV[,'units'])

mean(rho)
HPDinterval(rho)

#model per each gypsum affinity + soil
test.data<-data.frame(variable=data2$N, country=data2$country, factor=data2$statfix, factor2=data2$statfix2,soil=data2$soilN, phylo=data2$phylo, species=data2$species)
model2<-MCMCglmm(variable~-1+factor+soil, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model2)

#model between gypsum affinity + soil
model3<-MCMCglmm(variable~factor+soil, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model3)

model4<-MCMCglmm(variable~factor2+soil, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model4)

#model per each region
#wide
data3<-test.data[which(test.data$factor=="wide"), ]
model5<-MCMCglmm(variable~-1+country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=data3,nitt=150000,burnin=1000,thin=50) 
summary(model5)

model6<-MCMCglmm(variable~country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model6)

#narrow
data3<-test.data[which(test.data$factor=="narrow"), ]
model5<-MCMCglmm(variable~-1+country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=data3,nitt=150000,burnin=1000,thin=50) 
summary(model5)

model6<-MCMCglmm(variable~country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model6)

#gypsovags
data3<-test.data[which(test.data$factor=="gypsovag"), ]
model5<-MCMCglmm(variable~-1+country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=data3,nitt=150000,burnin=1000,thin=50) 
summary(model5)

model6<-MCMCglmm(variable~country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model6)


###COARSE ROOTS
tip_coarse<-c("Bouteloua_breviseta","Bouteloua_curtipendula","Gaillardia_pulchella","Sporobolus_nealleyi","Sporobolus_sp")
tree_coarse<-drop.tip(tree,tip_coarse)
plot(tree_coarse)

phylo.inv<-inverseA(tree_coarse,nodes="TIPS",scale=TRUE)
prior2<-list(G=list(G1=list(V=1,nu=0.02),G2=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))

data2<-data[which(data$fractfix=="coarse"), ]
data2$species<-as.factor(data2$species)

test.data<-data.frame(variable=data2$N, phylo=data2$phylo, species=data2$species)

model1<-MCMCglmm(variable~1, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model1)

lambda <- model1$VCV[,'phylo']/
  (model1$VCV[,'phylo']+model1$VCV[,'species']+
     model1$VCV[,'units'])

mean(lambda)
HPDinterval(lambda)

gamma<-model1$VCV[,'species']/
  (model1$VCV[,'phylo']+model1$VCV[,'species']+
     model1$VCV[,'units'])

mean(gamma)
HPDinterval(gamma)

rho<-model1$VCV[,'units']/
  (model1$VCV[,'phylo']+model1$VCV[,'species']+
     model1$VCV[,'units'])

mean(rho)
HPDinterval(rho)

#model per each gypsum affinity + soil
test.data<-data.frame(variable=data2$N,country=data2$country, factor=data2$statfix, factor2=data2$statfix2,soil=data2$soilN, phylo=data2$phylo, species=data2$species)
model2<-MCMCglmm(variable~-1+factor+soil, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model2)

#model between gypsum affinity + soil
model3<-MCMCglmm(variable~factor+soil, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model3)

model4<-MCMCglmm(variable~factor2+soil, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model4)

#model per each region
#wide
data3<-test.data[which(test.data$factor=="wide"), ]
model5<-MCMCglmm(variable~-1+country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=data3,nitt=150000,burnin=1000,thin=50) 
summary(model5)

model6<-MCMCglmm(variable~country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model6)

#narrow
data3<-test.data[which(test.data$factor=="narrow"), ]
model5<-MCMCglmm(variable~-1+country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=data3,nitt=150000,burnin=1000,thin=50) 
summary(model5)

model6<-MCMCglmm(variable~country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model6)

#gypsovags
data3<-test.data[which(test.data$factor=="gypsovag"), ]
model5<-MCMCglmm(variable~-1+country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=data3,nitt=150000,burnin=1000,thin=50) 
summary(model5)

model6<-MCMCglmm(variable~country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model6)

###FINE ROOTS
tip_fine<-c("Abronia_angustifolia","Abronia_fragrans","Acleisanthes_longiflora","Anulocaulis_eriosolenus","Anulocaulis_leiosolenus_gyp","Anulocaulis_leiosolenus_las","Nerisyrenia_linearifolia","Nyctaginia_capitata","Ononis_tridentata")
tree_fine<-drop.tip(tree,tip_fine)
plot(tree_fine)

phylo.inv<-inverseA(tree_fine,nodes="TIPS",scale=TRUE)
prior2<-list(G=list(G1=list(V=1,nu=0.02),G2=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))

data2<-data[which(data$fractfix=="fine"), ]
data2$species<-as.factor(data2$species)

test.data<-data.frame(variable=data2$N, phylo=data2$phylo, species=data2$species)

model1<-MCMCglmm(variable~1, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model1)

lambda <- model1$VCV[,'phylo']/
  (model1$VCV[,'phylo']+model1$VCV[,'species']+
     model1$VCV[,'units'])

mean(lambda)
HPDinterval(lambda)

gamma<-model1$VCV[,'species']/
  (model1$VCV[,'phylo']+model1$VCV[,'species']+
     model1$VCV[,'units'])

mean(gamma)
HPDinterval(gamma)

rho<-model1$VCV[,'units']/
  (model1$VCV[,'phylo']+model1$VCV[,'species']+
     model1$VCV[,'units'])

mean(rho)
HPDinterval(rho)

#model per each gypsum affinity + soil
test.data<-data.frame(variable=data2$N, country=data2$country, factor=data2$statfix, factor2=data2$statfix2,soil=data2$soilP, phylo=data2$phylo, species=data2$species)
model2<-MCMCglmm(variable~-1+factor+soil, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model2)

#model between gypsum affinity + soil
model3<-MCMCglmm(variable~factor+soil, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model3)

model4<-MCMCglmm(variable~factor2+soil, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model4)

#model per each region
#wide
data3<-test.data[which(test.data$factor=="wide"), ]
model5<-MCMCglmm(variable~-1+country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=data3,nitt=150000,burnin=1000,thin=50) 
summary(model5)

model6<-MCMCglmm(variable~country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model6)

#narrow
data3<-test.data[which(test.data$factor=="narrow"), ]
model5<-MCMCglmm(variable~-1+country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=data3,nitt=150000,burnin=1000,thin=50) 
summary(model5)

model6<-MCMCglmm(variable~country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model6)

#gypsovags
data3<-test.data[which(test.data$factor=="gypsovag"), ]
model5<-MCMCglmm(variable~-1+country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=data3,nitt=150000,burnin=1000,thin=50) 
summary(model5)

model6<-MCMCglmm(variable~country, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model6)


###ORGANS


#model with organs
test.data<-data.frame(variable=data$N, factor=data$fractfix,factor2=data$fractfix2, factor3=data$fractfix3,phylo=data$phylo, species=data$species)
model5<-MCMCglmm(variable~-1+factor, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model5)

model6<-MCMCglmm(variable~factor, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model6)

model7<-MCMCglmm(variable~factor2, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model7)

model8<-MCMCglmm(variable~factor3, random=~phylo+species,family="gaussian", ginverse=list(phylo=phylo.inv$Ainv),prior=prior2,data=test.data,nitt=150000,burnin=1000,thin=50) 
summary(model8)

##PLOTS

##posterior mean
##graph Figure 1
data <- read_excel("results.xlsx", sheet="results_regionsWV")
data<-data.frame(data)
data$organ <- factor(data$organ, levels = c("leaf", "stem", "coarse", "fine"))


est <- function(data) {
  ggplot(data=data, mapping=aes(x=type2, y=mean, ymin=upper, ymax=lower))+
    geom_pointrange(color="black", fill="black", shape=19)+
    geom_hline(yintercept = 0)+
    scale_x_discrete(limits=c("WM","WN","VM","VN"), labels=c("Med","N.America","Med","N.America"))+
    theme(
      legend.position="none",
      strip.text.x = element_text(size=15),
      strip.background = element_rect(colour="black", fill="white"),
      axis.text=element_text(size=10, colour="black"),
      axis.title.x = element_blank(),
      axis.title.y=element_blank(),
      legend.title = element_blank(), 
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.background = element_rect(fill = "white", colour = "grey50"))+
    facet_wrap(~organ, nrow=4, strip.position="right")+
    coord_flip()
}

data2<-data[which(data$variable=="S"), ]
sulphur<-est(data2)

data2<-data[which(data$variable=="Ca"), ]
calcium<-est(data2)

data2<-data[which(data$variable=="Mg"), ]
magnesium<-est(data2)

data2<-data[which(data$variable=="K"), ]
potassium<-est(data2)

data2<-data[which(data$variable=="P"), ]
phosphorus<-est(data2)

data2<-data[which(data$variable=="N"), ]
nitrogen<-est(data2)

require(gridExtra)
grid.arrange(calcium, magnesium, sulphur,potassium, phosphorus, nitrogen, nrow=1,ncol=6)



###phylogenetic signal per elements

library(picante)
library(tidyverse)

data1<-melt(data, measure.vars=c("Ca","N","S","Mg","P","K"))
data2<-data.frame(data1$species,data1$type, data1$fractfix,data1$variable,data1$value)
colnames(data2)<-c("species","type","fractfix","variable","value")
tgc1 <- summarySE(na.omit(data2), measurevar=c("value"), groupvars=c("variable","species","fractfix"))
tgc1<-data.frame(tgc1)
tgc2<-dcast(tgc1, species ~ variable+fractfix)

rownames(tgc2)<-tgc2$species

#tree of target species
tree <- read.tree("Clare_tree36_tree.tre", keep.multi = FALSE)
plot(tree)



#dendogram
library(phytools)
value<-setNames(tgc2$S_coarse,rownames(tgc2))
fitEB<-anc.ML(tree,value)
fitEB
a<-contMap(tree,value, method="user",anc.states=fitEB$ace,plot=FALSE)
a<-setMap(a, colors=c("gray90","darkred"))
plot(a)

#drop tip for coarse roots and fine roots
library(ape)
tip_coarse<-c("Bouteloua_breviseta","Bouteloua_curtipendula","Gaillardia_pulchella","Sporobolus_nealleyi","Sporobolus_sp")
tree2<-drop.tip(tree,tip_coarse)
plot(tree2)
value3<-na.omit(value2<-data.frame(tgc2$S_coarse,tgc2$species))
colnames(value3)<-c("coarse","species")
rownames(value3)<-value3$species
value<-setNames(value3$coarse,rownames(value3))

fitEB<-anc.ML(tree2,value)
fitEB
a<-contMap(tree2,value, method="user",anc.states=fitEB$ace,plot=FALSE)
a<-setMap(a, colors=c("gray90","darkred"))
plot(a)

tip_fine<-c("Abronia_angustifolia","Abronia_fragrans","Acleisanthes_longiflora","Anulocaulis_eriosolenus","Anulocaulis_leiosolenus_gyp","Anulocaulis_leiosolenus_las","Nerisyrenia_linearifolia","Nyctaginia_capitata","Ononis_tridentata")
tree3<-drop.tip(tree,tip_fine)
plot(tree3)
value3<-na.omit(value2<-data.frame(tgc2$N_fine,tgc2$species))
colnames(value3)<-c("fine","species")
rownames(value3)<-value3$species
value<-setNames(value3$fine,rownames(value3))

fitEB<-anc.ML(tree3,value)
fitEB
a<-contMap(tree3,value, method="user",anc.states=fitEB$ace,plot=FALSE)
a<-setMap(a, colors=c("gray90","darkred"))
plot(a)
