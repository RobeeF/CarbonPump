#####Pour le rapport :
#Dataframe
Results <- data.frame(groupe = c(rep("Apports",2), rep("Production",3),rep("Respiration",3),"Export"),
                      Name = c("POC","DOC","Procaryotes libres","Procaryotes fixés","Zooplancton","Procaryotes libres","Procaryotes fixés","Zooplancton","POC"),
                      Valeurs = c(POC,DOC,ProdBf,ProdBatt,ProdZoo,RespBf,RespBatt,respZoo,bathy))
head(Results)

Results$Name <- factor(Results$Name, levels = c("POC","DOC","Procaryotes libres","Procaryotes fixés","Zooplancton"))
#cette etape sert juste e les placer dans l'ordre  qu'on veut pour GGPLOT qui va les placer par ordfe alphabetique si cest juste des chaine de charact

#Choisir tranche profondeur:
tranche = "28-500m"


library(ggplot2)

###Histo empile

ggplot(data=Results, aes(x=groupe, y=Valeurs, fill = Name)) +
  geom_bar(stat="identity",color="black")+
  ggtitle("Carbon budget for the mesopelagic zone") +
  ylab(expression(paste(" Flux en mg C ", m^-2, d^-1)))+
  xlab(tranche)+
  theme_bw()+
  ylim(0,1050)+
  scale_x_discrete(limits=c("Apports","Production","Respiration"))+#met dans l'ordre les barres
  scale_fill_manual(values=c("green3","chartreuse4","palevioletred1","deeppink4","lemonchiffon4"))
theme(plot.title = element_text(size = 10, face = "plain"),
      axis.text.x= element_text(size=10),
      axis.text.y= element_text(size=10),
      axis.title.x = element_blank(),
      legend.title=element_blank())