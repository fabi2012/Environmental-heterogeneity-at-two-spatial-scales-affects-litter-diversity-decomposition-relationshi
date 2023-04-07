###Main analysis of  Environmental heterogeneity at two spatial scales affects litter diversity-decomposition relationships

library(vegan)
library(tidyverse)
library(car)
library(lmerTest)


####read the data
setwd("C:/Users/sriva/Dropbox/grad students and PDFs/Fabiola/Decomposition")
macroscale <- read.csv(file.choose(),head=T, sep=",") #spatialscales

macroscale <-  mutate(macroscale, decomposition=(weightinitial - weightend)*100/weightinitial,
                      tasa= -(log10(weightend/weightinitial))/(120/365),
                      alnus=(Alnusinitial - Alnusend)/Alnusinitial, 
                      piper= (Piperinitial - Piperend)/Piperinitial,
                      croton=(Crotoninitial - Crotonend)/Crotoninitial,
                      block = sample+12*(elevation == "Middle")+24*(elevation=="High"),
                      block = as.factor(block),
                      litter_rich = (Alnusinitial>0)+ (Piperinitial>0)+(Crotoninitial>0))
str(macroscale)
#######Decomposition-------------

#litter decomposition 

decomposition.glmm<- lmer(decomposition~ elevation*treatment*habitat+(1|block), data=macroscale,REML=TRUE)
Anova(decomposition.glmm ,3, test="Chi")

decomposition.glmm2<-lmer(decomposition~litter_rich*elevation*habitat+(1|block)+(1|treatment), 
                          data=macroscale,REML=TRUE)
Anova(decomposition.glmm2 ,3, test="Chi")


##Alnus
treatment.alnus <- filter(macroscale, Alnusinitial>0)

alnus.glmm<- lmer(decomposition~ elevation*treatment*habitat+(1|block), 
                  data=treatment.alnus,REML=TRUE)
Anova(alnus.glmm ,3, test="Chi")#effects of litter composition , elevation, habitat

alnus.glmm2<-lmer(decomposition~litter_rich*elevation*habitat+(1|block)+(1|treatment), data=treatment.alnus,REML=TRUE)
Anova(alnus.glmm2 ,3, test="Chi") #no effect of litter richness on Alnus decomposition

#Piper
treatment.piper <- filter(macroscale, Piperinitial>0)

piper.glmm<- lmer(decomposition~ elevation*treatment*habitat+(1|block), 
                  data=treatment.piper ,REML=TRUE)
Anova(piper.glmm ,3, test="Chi")#effects of litter composition , elevation, habitat
#higher decomp with croton than alnus, lower all three litter together
#but depends on elevation, habitat

piper.glmm2<-lmer(decomposition~litter_rich*elevation*habitat+(1|block)+(1|treatment), 
                  data=treatment.piper ,REML=TRUE)
Anova(piper.glmm2 ,3, test="Chi") #litter richness x habitat affect Piper decomp

##croton
treatment.croton <- filter(macroscale, Crotoninitial>0)

croton.glmm<- lmer(decomposition~ elevation*treatment*habitat+(1|block), 
                   data=treatment.croton ,REML=TRUE)
Anova(croton.glmm ,3, test="Chi")#effects of litter composition, elevation

croton.glmm2<-lmer(decomposition~litter_rich*elevation*habitat+(1|block)+(1|treatment), 
                   data=treatment.croton ,REML=TRUE)
Anova(croton.glmm2 ,3, test="Chi") #no effects of litter rich on Croton decomp


### environmental variables --------
#temperature

tem.glmm<- lmer(temperature~ elevation*treatment*habitat+(1|block), 
                data=macroscale ,REML=TRUE)
Anova(tem.glmm ,3, test="Chi")#habitat only

tem.glmm2<-lmer(temperature~litter_rich*elevation*habitat+(1|block)+(1|treatment), 
                data=macroscale,REML=TRUE)
Anova(tem.glmm2 ,3, test="Chi") #habitat only

#pH
ph.glmm<- lmer(pH~ elevation*treatment*habitat+(1|block), 
               data=macroscale ,REML=TRUE)
Anova(ph.glmm ,3, test="Chi") #habitat

ph.glmm2<-lmer(pH~litter_rich*elevation*habitat+(1|block)+(1|treatment),
               data=macroscale ,REML=TRUE)
Anova(ph.glmm2 ,3, test="Chi") #habitat, habitat*richness


#total solid dissolved

conductivity.glmm<- lmer(conductivity~ elevation*treatment*habitat+(1|block),
                         data=macroscale ,REML=TRUE)
Anova(conductivity.glmm ,3, test="Chi") #habitat

conductivity.glmm2<- lmer(conductivity~ litter_rich*elevation*habitat+(1|block)+(1|treatment), 
                          data=macroscale ,REML=TRUE)
Anova(conductivity.glmm2 ,3, test="Chi") #habitat only


########Diversity effects on decomposition (selection and complementary effects)----------
#arrange the data for analysis

data.diversity<- macroscale[,c(1:4,18:21)] %>% 
  gather(species,decompositionvalue,alnus:croton, na.rm=F)

Y.data <- filter(data.diversity,treatment=="mixA"|treatment=="mixB"|treatment=="mixC"|treatment=="mixD")
names(Y.data)[7]<-"Y"

M.data <- filter(data.diversity,treatment=="nomixA"|treatment=="nomixB"|treatment=="nomixC")
names(M.data)[7]<-"M"

df <- (rbind(cbind(filter(Y.data,treatment=="mixA"),filter(M.data,treatment=="nomixA")),
             cbind(filter(Y.data,treatment=="mixA"),filter(M.data,treatment=="nomixB")),
             cbind(filter(Y.data,treatment=="mixB"),filter(M.data,treatment=="nomixA")),
             cbind(filter(Y.data,treatment=="mixB"),filter(M.data,treatment=="nomixC")),
             cbind(filter(Y.data,treatment=="mixC"),filter(M.data,treatment=="nomixB")),
             cbind(filter(Y.data,treatment=="mixC"),filter(M.data,treatment=="nomixC")),
             cbind(filter(Y.data,treatment=="mixD"),filter(M.data,treatment=="nomixA")),
             cbind(filter(Y.data,treatment=="mixD"),filter(M.data,treatment=="nomixB")),
             cbind(filter(Y.data,treatment=="mixD"),filter(M.data,treatment=="nomixC")))[,c(1:7,14)])
df <- na.omit(df)
df<-df %>% 
  mutate(sprich = 2+(treatment=="mixD"),
         percent.M = M, 
         percent.Y = Y,
         M = M*3, #must convert to absolute g litter lost, not relative %
         Y = Y*(3/sprich)) %>% #must convert to absolute g litter lost, not relative %
  unite("hstb", c(elevation,habitat, treatment,block), remove = F) %>% 
  mutate(hstb = as.factor(hstb))

#first a function that calculates SE and CE in treehole subsets individually

funCESE <- function(df,n,tr,sh) {
  test <- df %>% filter(block == n, treatment==tr, habitat == sh)
  print(test)
  test$RYe <- 1/nrow(test) 
  test$RYo <- test$Y/test$M
  test$chRY <- test$RYo-test$RYe
  CE = mean(test$M)*mean(test$chRY)*nrow(test)
  SE = cov(test$M,test$chRY)*(nrow(test)-1)
  outdf<-as.data.frame(cbind(CE,SE))
  return(outdf)} 

funCESE(df=df,"1", "mixD", "edge")
funCESE(df=df,"2", "mixA", "edge")
funCESE(df=df,"25", "mixA", "edge")

#now this function for entire dataset
funCESE2 <- function(test) {
  test$RYe <- 1/nrow(test) 
  test$RYo <- test$Y/test$M
  test$chRY <- test$RYo-test$RYe
  CE = mean(test$M)*mean(test$chRY)*nrow(test)
  SE = cov(test$M,test$chRY)*(nrow(test)-1)
  NBE = sum(test$M*test$chRY)
  outdf<- cbind(CE,SE,NBE)
  return(outdf)} 


#calculate the diversity effects ( selections  and complementarity)
####script of Isbell F, Cowles J, Dee LE, et al. Quantifying effects of biodiversity on ecosystem functioning across times and places. Ecol Lett. 2018;21(6):763–778. doi:10.1111/ele.12928######

funcCE <- function(df) {return(data.frame(CE = mean(df$M)*mean(df$chRY)*nrow(df)))} 
funcSE <- function(df) {return(data.frame(SE = cov(df$M,df$chRY)*(nrow(df)-1)))}

fune <- function(df) {
  df$RYe <- 1/df$sprich
  df$RYo <- df$Y/df$M
  df$chRY <- df$RYo-df$RYe
  return(data.frame(
    NBE = sum(df$M*df$chRY), #net biodiversity effect
    TC = mean(df$M)*mean(df$chRY)*nrow(df), #total complementarity effect
    total.se = cov(df$M,df$chRY)*(nrow(df)-1), #total selection effect
  ))} 


#result.diversity <- ddply(df,.(habitat,subhabitat, treatment,block), fune)### lo logreeee.jijiji
#the above works, but requires plyr which conflicts with dplyr
#write_csv(result.diversity, "result_diversity.csv")

#the following avoid plyr and gives identical output:

result.diversity2 <- df %>% 
  group_by(hstb) %>% 
  do(data.frame(funCESE2(.))) %>% 
  ungroup() %>% 
  separate(col = hstb, into = c("elevation", "habitat", "treatment", "block"), sep = "_")

meanx<-function(x) {mean(x, na.rm=TRUE)}


result.diversity2 <- result.diversity2%>% mutate(litter_rich=
                                                   case_when(treatment =="mixA" ~ 2, 
                                                             treatment =="mixB" ~ 2,
                                                             treatment =="mixC" ~ 2,
                                                             treatment =="mixD" ~ 3 ))


tapply(result.diversity2$CE, result.diversity2$treatment, meanx)
#large positive mixA, negative mixC

tapply(result.diversity2$SE, result.diversity2$treatment, meanx)
#everything negative esp mixA, mixD

tapply(result.diversity2$NBE, result.diversity2$treatment, meanx)
#mixA very positive, mixC very negative

#selection effect

Selection2<-lmer(SE~elevation*treatment*habitat+(1|block),data=result.diversity2)
Anova(Selection2,2,test="Chi") #treatment, habitat, treat x habitat

tapply(result.diversity2$SE, result.diversity2$habitat, meanx)
tapply(result.diversity2$SE, result.diversity2$treatment, meanx)

#complementary effect

Complementary2 <-lmer(CE~elevation*treatment*habitat +(1|block),
                      data=result.diversity2)
Anova(Complementary2,2,test="Chi") #treatment, elevation x treat, helevation x habitat...now treatment x subhabitat, habitat

par(mfrow=c(2,2))
plot(Complementary2)

tapply(result.diversity2$CE, result.diversity2$habitat, meanx)
tapply(result.diversity2$CE, result.diversity2$treatment, meanx)

#net effect

Net2 <-lmer(NBE~elevation* treatment*habitat +(1|block),data=result.diversity2)
Anova(Net2,2,test="Chi")#treatment,habitat,elevation x treatment, habitat x elevation...now treatment x subhabitat

tapply(result.diversity2$NBE, result.diversity2$elevation, meanx)
#     High         Low      Middle 
# 0.19270175 -0.30140244  0.05947464 

tapply(result.diversity2$NBE, result.diversity2$treatment, meanx)
#mixA        mixB        mixC        mixD 
#0.31273438  0.07154930 -0.34428571 -0.03018229 

##### partitioning biodiversity effects to look at spatial insurance

#first lets look at three-species treatment

dfD<-df %>% 
  filter(treatment=="mixD") %>% 
  unite("sub_species", c(habitat,elevation, species), remove = FALSE) %>% 
  mutate(sub_species = as.factor(sub_species))

df2D<-dfD %>% 
  group_by(sub_species) %>% 
  summarise(sub_species = first(sub_species),
            habitat = first(habitat),
            elevation = first(elevation),
            species = first(species),
            Y = mean(Y, na.rm=TRUE),
            M = mean(M, na.rm=TRUE),
            sprich = mean(sprich, na.rm=TRUE))


fune3sp <- function(df) {
  #add columns to df that are needed for quantifying biodiversity effects df$RYo <- df$Y/df$M
  df$RYe <- 1/df$sprich 
  df$RYo <- df$Y/df$M
  df$chRY <- df$RYo-df$RYe
  jdf <- aggregate.data.frame(df$Y,by=list(df$elevation,df$habitat),sum)
  names(jdf) <- c("elevation","habitat","Yo")
  df <- left_join(df,jdf) 
  df$p <- df$Y/df$Yo
  df$chRYo <- df$RYo-df$p
  df$dp <- df$p-df$RYe
  jdf <- aggregate.data.frame(df$M,by=list(df$species,df$elevation),mean) 
  names(jdf) <- c("species","elevation","Mst") #here elevation is taking the place of "time" in original Isbell script
  df <- left_join(df,jdf)
  jdf <- aggregate.data.frame(df$M,by=list(df$species,df$habitat),mean) 
  names(jdf) <- c("species","habitat","Msp")#here habitat stands in for place
  df <- left_join(df,jdf)
  jdf <- aggregate.data.frame(df$dp,by=list(df$species,df$elevation),mean)
  names(jdf) <- c("species","elevation","pst")
  df <- left_join(df,jdf)
  jdf <- aggregate.data.frame(df$dp,by=list(df$species,df$habitat),mean) 
  names(jdf) <- c("species","habitat","psp")
  df <- left_join(df,jdf)
  df <- df[order(df$elevation,df$habitat,df$species),]
  return(data.frame(
    NBE = sum(df$M*df$chRY), #net biodiversity effect
    TC = mean(df$M)*mean(df$chRY)*nrow(df), #total complementarity effect
    total.se = cov(df$M,df$chRY)*(nrow(df)-1), #total selection effect
    NO = cov(df$M,df$chRYo)*(nrow(df)-1), #nonrandom overyielding effect
    IE = cov(df$M,df$dp)*(nrow(df)-1), #total insurance effect: IE = AS+TI+SI+ST
    AS = sum((tapply(df$M,df$species,mean)-mean(df$M))*(tapply(df$dp,df$species,mean)-
                                                          mean(df$dp)))*2*3, #average selection effect (note 2 = P = habitat, 3 = T = elevation)
    TI = sum((as.vector(tapply(df$M,list(df$species,df$elevation),mean))-
                rep(tapply(df$M,df$species,mean),3))*(as.vector(tapply(df$dp,list(df$species,df$elevation),
                                                                       mean))-rep(tapply(df$dp,df$species,mean),3)))*2, #elevation insurance effect (note 2 = P)
    SI = sum((as.vector(tapply(df$M,list(df$species,df$habitat),mean))-
                rep(tapply(df$M,df$species,mean),2))*(as.vector(tapply(df$dp,list(df$species,df$habitat),
                                                                       mean))-rep(tapply(df$dp,df$species,mean),2)))*3, #habitat insurance effect..first and second "2" is the number of habitats (P), and the last "3" = T (three elevations)
    ST = sum((df$M-df$Mst-df$Msp+rep(tapply(df$M,df$species,mean),6)+rep(mean(df$M),18))*(df$dp-df$pst-df$psp+rep(tapply(df$dp,df$species,mean),6)+rep(mean(df$dp),18))) ))} #spatiotemporal insurance effect


fune3sp(df2D)

#  NBE       TC   total.se         NO         IE         AS
# 0.2391385 1.149621 -0.9104824 -0.6929399 -0.2175425 -0.3434348

# TI           SI         ST
# 0.1229727 -0.002439746 0.00535935

dfA<-df %>% 
  filter(treatment=="mixA") %>% 
  unite("sub_species", c(habitat,elevation, species), remove = FALSE) %>% 
  mutate(sub_species = as.factor(sub_species)) %>% 
  group_by(sub_species) %>% 
  summarise(sub_species = first(sub_species),
            habitat = first(habitat),
            elevation = first(elevation),
            species = first(species),
            Y = mean(Y, na.rm=TRUE),
            M = mean(M, na.rm=TRUE),
            sprich = mean(sprich, na.rm=TRUE))

dfB<-df %>% 
  filter(treatment=="mixB") %>% 
  unite("sub_species", c(habitat,elevation, species), remove = FALSE) %>% 
  mutate(sub_species = as.factor(sub_species)) %>% 
  group_by(sub_species) %>% 
  summarise(sub_species = first(sub_species),
            habitat = first(habitat),
            elevation = first(elevation),
            species = first(species),
            Y = mean(Y, na.rm=TRUE),
            M = mean(M, na.rm=TRUE),
            sprich = mean(sprich, na.rm=TRUE))

dfC<-df %>% 
  filter(treatment=="mixC") %>% 
  unite("sub_species", c(habitat,elevation, species), remove = FALSE) %>% 
  mutate(sub_species = as.factor(sub_species)) %>% 
  group_by(sub_species) %>% 
  summarise(sub_species = first(sub_species),
            habitat = first(habitat),
            elevation = first(elevation),
            species = first(species),
            Y = mean(Y, na.rm=TRUE),
            M = mean(M, na.rm=TRUE),
            sprich = mean(sprich, na.rm=TRUE))

fune2sp <- function(df) {
  #add columns to df that are needed for quantifying biodiversity effects df$RYo <- df$Y/df$M
  df$RYe <- 1/df$sprich 
  df$RYo <- df$Y/df$M
  df$chRY <- df$RYo-df$RYe
  jdf <- aggregate.data.frame(df$Y,by=list(df$elevation,df$habitat),sum)
  names(jdf) <- c("elevation","habitat","Yo")
  df <- left_join(df,jdf) 
  df$p <- df$Y/df$Yo
  df$chRYo <- df$RYo-df$p
  df$dp <- df$p-df$RYe
  jdf <- aggregate.data.frame(df$M,by=list(df$species,df$elevation),mean) 
  names(jdf) <- c("species","elevation","Mst") #here elevationis taking the place of "time" in original Isbell script
  df <- left_join(df,jdf)
  jdf <- aggregate.data.frame(df$M,by=list(df$species,df$habitat),mean) 
  names(jdf) <- c("species","habitat","Msp")#here habitat stands in for place
  df <- left_join(df,jdf)
  jdf <- aggregate.data.frame(df$dp,by=list(df$species,df$elevation),mean)
  names(jdf) <- c("species","elevation","pst")
  df <- left_join(df,jdf)
  jdf <- aggregate.data.frame(df$dp,by=list(df$species,df$habitat),mean) 
  names(jdf) <- c("species","habitat","psp")
  df <- left_join(df,jdf)
  df <- df[order(df$elevation,df$habitat,df$species),]
  return(data.frame(
    NBE = sum(df$M*df$chRY), #net biodiversity effect
    TC = mean(df$M)*mean(df$chRY)*nrow(df), #total complementarity effect
    total.se = cov(df$M,df$chRY)*(nrow(df)-1), #total selection effect
    NO = cov(df$M,df$chRYo)*(nrow(df)-1), #nonrandom overyielding effect
    IE = cov(df$M,df$dp)*(nrow(df)-1), #total insurance effect: IE = AS+TI+SI+ST
    AS = sum((tapply(df$M,df$species,mean)-mean(df$M))*(tapply(df$dp,df$species,mean)-
                                                          mean(df$dp)))*2*3, #average selection effect (note 2 = P = habitat, 3 = T = elevation)
    TI = sum((as.vector(tapply(df$M,list(df$species,df$elevation),mean))-
                rep(tapply(df$M,df$species,mean),3))*(as.vector(tapply(df$dp,list(df$species,df$elevation),
                                                                       mean))-rep(tapply(df$dp,df$species,mean),3)))*2, #elevation insurance effect (note 2 = P)
    SI = sum((as.vector(tapply(df$M,list(df$species,df$habitat),mean))-
                rep(tapply(df$M,df$species,mean),2))*(as.vector(tapply(df$dp,list(df$species,df$habitat),
                                                                       mean))-rep(tapply(df$dp,df$species,mean),2)))*3, #habitat insurance effect..first and second "2" is the number of habitats (P), and the last "3" = T (three elevations)
    ST = sum((df$M-df$Mst-df$Msp+rep(tapply(df$M,df$species,mean),6)+rep(mean(df$M),12))*(df$dp-df$pst-df$psp+rep(tapply(df$dp,df$species,mean),6)+rep(mean(df$dp),12))) ))} #spatiotemporal insurance effect


fune2sp(dfA) #mixA = AP

# NBE       TC      total.se         NO         IE         AS
# 1.962828 3.042309 -1.079481 -0.9303962 -0.1490845 -0.1408265
# TI           SI           ST
# 0.002960246 -0.009650038 -0.001568209

fune2sp(dfB) #mixB = AC

# NBE       TC    total.se         NO        IE        AS
# 0.812697 1.469749 -0.6570516 -0.8603312 0.2032796 0.2103557

#TI           SI           ST
# -0.001842886 -0.001068182 -0.004165009

fune2sp(dfC) #mixC = CP

# NBE        TC     total.se         NO        IE        AS
# -1.979687 -1.788436 -0.191251 -0.3614089 0.1701579 0.2495491

#TI          SI           ST
# -0.03527782 -0.04202639 -0.002086953

###Figure litter richness------------

library(ggplot2)

my_sum<- macroscale%>%
  group_by(litter_rich,elevation,habitat) %>%
  summarise( 
    n=n(),
    mean=mean(decomposition,na.rm=T),
    sd=sd(decomposition,na.rm = T),
    mean_alnus=mean(alnus*100,na.rm=T),
    sd_alnus=sd(alnus*100,na.rm=T),
    mean_piper=mean(piper*100,na.rm=T),
    sd_piper=sd(piper*100,na.rm=T),
    mean_croton=mean(croton*100,na.rm=T),
    sd_croton=sd(croton*100,na.rm=T),
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))%>%
  mutate( se_alnus=sd_alnus/sqrt(n))  %>%
  mutate( ic_alnus=se_alnus * qt((1-0.05)/2 + .5, n-1))%>%
  mutate( se_croton=sd_croton/sqrt(n))  %>%
  mutate( ic_croton=se_croton * qt((1-0.05)/2 + .5, n-1))%>%
  mutate( se_piper=sd_piper/sqrt(n))  %>%
  mutate( ic_piper=se_piper * qt((1-0.05)/2 + .5, n-1))%>%
  mutate(elevation.2=
           case_when(elevation =="High" ~ 3,
                     elevation =="Low" ~ 1,
                     elevation =="Middle" ~ 2))

figure1A <- ggplot(my_sum,
                   aes(x=as.factor(litter_rich),
                       y=mean,
                       color = as.factor(elevation.2),
                       shape=habitat)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  xlab("Litter diversity") +
  ylab("Litter mass loss (%)") +
  scale_y_continuous(limit=c(0,100),breaks = seq(0, 100, by = 10)) +
  scale_shape_manual(name = "Habitat  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Elevation",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))

figure1B <- ggplot(my_sum,
                   aes(x=as.factor(litter_rich),
                       y=mean_alnus,
                       color = as.factor(elevation.2),
                       shape=habitat)) +
  geom_errorbar(aes(ymin=mean_alnus-se_alnus, ymax=mean_alnus+se_alnus),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  xlab("Litter diversity") +
  ylab(expression(paste("Litter mass loss of"," ",italic ("A. acuminata")," ","(%)")))+
  scale_y_continuous(limit=c(0,100),breaks = seq(0, 100, by = 10)) +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))

figure1C <- ggplot(my_sum,
                   aes(x=as.factor(litter_rich),
                       y=mean_piper,
                       color = as.factor(elevation.2),
                       shape=habitat)) +
  geom_errorbar(aes(ymin=mean_piper-se_piper, ymax=mean_piper+se_piper),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  xlab("Litter diversity") +
  ylab(expression(paste("Litter mass loss of"," ",italic ("P. imperiale")," ","(%)")))+
  scale_y_continuous(limit=c(0,100),breaks = seq(0, 100, by = 10)) +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))

figure1D <- ggplot(my_sum,
                   aes(x=as.factor(litter_rich),
                       y=mean_croton,
                       color = as.factor(elevation.2),
                       shape=habitat)) +
  geom_errorbar(aes(ymin=mean_croton-se_croton, ymax=mean_croton+se_croton),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  xlab("Litter diversity") +
  ylab(expression(paste("Litter mass loss of"," ",italic ("C. magdalenensis")," ","(%)")))+
  scale_y_continuous(limit=c(0,100),breaks = seq(0, 100, by = 10)) +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))

tiff("Figure 1corregida.tiff",
     height = 9,
     width = 7,
     unit="in",
     res = 600)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend<-g_legend(figure1A+theme(legend.position="bottom"))

library(gridExtra)
library(gridGraphics)

grid.arrange(arrangeGrob(figure1A+
                           ggtitle("a") +
                           theme(legend.position = "none",
                                 legend.title=element_text(size=rel(1)),
                                 legend.text=element_text(size=rel(1)),
                                 axis.title.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.text.y=element_text(size=rel(1)),
                                 plot.title = element_text(size=rel(1),face = "bold")),
                         figure1B+
                           ggtitle("b")+
                           theme(legend.position = "none",
                                 legend.title=element_text(size=rel(1)),
                                 legend.text=element_text(size=rel(1)),
                                 axis.title.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.text.y=element_text(size=rel(1)),
                                 plot.title = element_text(size=rel(1),face = "bold")),
                         figure1C+
                           ggtitle("c")+
                           theme(legend.position = "none",
                                 axis.title.x=element_text(size=rel(1),family="Arial"),
                                 axis.text.x=element_text(size=rel(1),family="Arial"),
                                 axis.text.y=element_text(size=rel(1)),
                                 plot.title = element_text(size=rel(1),face = "bold")),
                         figure1D+
                           ggtitle("d")+
                           theme(legend.position = "none",
                                 axis.title.x=element_text(size=rel(1)),
                                 axis.text.x=element_text(size=rel(1)),
                                 axis.text.y=element_text(size=rel(1)),
                                 plot.title = element_text(size=rel(1),face = "bold")),ncol = 2),
             mylegend,heights=c(10, 1))

dev.off()

   
##### Figure environmental conditions

my_sum_env<- macroscale%>%
  group_by(litter_rich,elevation,habitat) %>%
  summarise( 
    n=n(),
    mean_T =mean(temperature  ,na.rm=T),
    sd_T =sd(temperature  ,na.rm = T),
    mean_PH=mean(pH ,na.rm=T),
    sd_PH=sd(pH ,na.rm=T),
    mean_C=mean(conductivity,na.rm=T),
    sd_C=sd(conductivity,na.rm=T)) %>%
  mutate( se_T=sd_T/sqrt(n))  %>%
  mutate( se_PH=sd_PH/sqrt(n))  %>%
  mutate( se_C=sd_C/sqrt(n))  %>%
  mutate(elevation.2=
           case_when(elevation =="High" ~ 3,
                     elevation =="Low" ~ 1,
                     elevation =="Middle" ~ 2))

figureSA <- ggplot(my_sum_env,
                   aes(x=as.factor(litter_rich),
                       y=mean_T,
                       color = as.factor(elevation.2),
                       shape=habitat)) +
  geom_errorbar(aes(ymin=mean_T-se_T, ymax=mean_T+se_T),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  ylab(expression(paste("Water temperature"," ","("^o,"C",")")))+
  scale_y_continuous(limits=c(10,17),breaks = seq(0, 17, by = 1)) +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))

figureSB <- ggplot(my_sum_env,
                   aes(x=as.factor(litter_rich),
                       y=mean_PH,
                       color = as.factor(elevation.2),
                       shape=habitat)) +
  geom_errorbar(aes(ymin=mean_PH-se_PH, ymax=mean_PH+se_PH),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  ylab(expression(paste("Water pH")))+
  scale_y_continuous(limits=c(5,7),breaks = seq(5, 7, by = 0.5)) +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))

figureSC <- ggplot(my_sum_env,
                   aes(x=as.factor(litter_rich),
                       y=mean_C,
                       color = as.factor(elevation.2),
                       shape=habitat)) +
  geom_errorbar(aes(ymin=mean_C-se_C, ymax=mean_C+se_C),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  xlab("Litter diversity") +
  ylab(expression(paste("Water total dissolved solid "," ", "(","mgl"^"-1",")")))+
  scale_y_continuous(limits=c(0,600),breaks = seq(0, 600, by = 100)) +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))



tiff("FigureS1.tiff",
     height = 9,
     width = 7,
     unit="in",
     res = 600)

grid.arrange(arrangeGrob(figureSA+
                           ggtitle("A") +
                           theme(legend.position = "none",
                                 legend.title=element_text(size=rel(1)),
                                 legend.text=element_text(size=rel(1)),
                                 axis.title.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.text.y=element_text(size=rel(1)),
                                 plot.title = element_text(size=rel(1),face = "bold")),
                         figureSB+
                           ggtitle("B")+
                           theme(legend.position = "none",
                                 legend.title=element_text(size=rel(1)),
                                 legend.text=element_text(size=rel(1)),
                                 axis.title.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.text.y=element_text(size=rel(1)),
                                 plot.title = element_text(size=rel(1),face = "bold")),
                         figureSC+
                           ggtitle("C")+
                           theme(legend.position = "none",
                                 axis.title.x=element_text(size=rel(1)),
                                 axis.text.x=element_text(size=rel(1)),
                                 axis.text.y=element_text(size=rel(1)),
                                 plot.title = element_text(size=rel(1),face = "bold")),ncol = 1),
             mylegend,heights=c(10, 1))
dev.off()

##### Figures litter composition  ----------

my_sum_comp<- macroscale%>%
  group_by(treatment,elevation,habitat) %>%
  summarise( 
    n=n(),
    mean=mean(decomposition,na.rm=T),
    sd=sd(decomposition,na.rm = T),
    mean_alnus=mean(alnus*100,na.rm=T),
    sd_alnus=sd(alnus*100,na.rm=T),
    mean_piper=mean(piper*100,na.rm=T),
    sd_piper=sd(piper*100,na.rm=T),
    mean_croton=mean(croton*100,na.rm=T),
    sd_croton=sd(croton*100,na.rm=T),
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))%>%
  mutate( se_alnus=sd_alnus/sqrt(n))  %>%
  mutate( ic_alnus=se_alnus * qt((1-0.05)/2 + .5, n-1))%>%
  mutate( se_croton=sd_croton/sqrt(n))  %>%
  mutate( ic_croton=se_croton * qt((1-0.05)/2 + .5, n-1))%>%
  mutate( se_piper=sd_piper/sqrt(n))  %>%
  mutate( ic_piper=se_piper * qt((1-0.05)/2 + .5, n-1))%>%
  mutate(elevation.2=
           case_when(elevation =="High" ~ 3,
                     elevation =="Low" ~ 1,
                     elevation =="Middle" ~ 2))

figure4A <- ggplot(my_sum_comp,
                   aes(x=as.factor(treatment),
                       y=mean,
                       color = as.factor(elevation.2),
                       shape=habitat)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  xlab("Litter composition") +
  scale_x_discrete(limit = c("nomixA","nomixB","nomixC", "mixA","mixB","mixC", "mixD"),
                   labels = c("A","P","C","AP","AC","PC", "APC"),
                   expand = expansion(add = c(0.6)))+
  ylab("Litter mass loss (%)") +
  scale_y_continuous(limit=c(0,100),breaks = seq(0, 100, by = 10)) +
  scale_shape_manual(name = "Habitat  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Elevation",
                     labels = c("Low","Midd","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))

figure4B <- ggplot(my_sum_comp,
                   aes(x=as.factor(treatment),
                       y=mean_alnus,
                       color = as.factor(elevation.2),
                       shape=habitat)) +
  geom_errorbar(aes(ymin=mean_alnus-se_alnus, ymax=mean_alnus+se_alnus),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  xlab("Litter diversity") +
  scale_x_discrete(limit = c("nomixA","mixA","mixB", "mixD"),
                   labels = c("A","AP","AC", "APC"),
                   expand = expansion(add = c(0.6)))+
  ylab(expression(paste("Litter mass loss of"," ",italic ("A. acuminata")," ","(%)")))+
  scale_y_continuous(limit=c(0,100),breaks = seq(0, 100, by = 10)) +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))

figure4C <- ggplot(my_sum_comp,
                   aes(x=as.factor(treatment),
                       y=mean_piper,
                       color = as.factor(elevation.2),
                       shape=habitat)) +
  geom_errorbar(aes(ymin=mean_piper-se_piper, ymax=mean_piper+se_piper),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  xlab("Litter composition") +
  scale_x_discrete(limit = c("nomixB","mixA","mixC", "mixD"),
                   labels = c("P","AP","PC", "APC"),
                   expand = expansion(add = c(0.6)))+
  ylab(expression(paste("Litter mass loss of"," ",italic ("P. imperiale")," ","(%)")))+
  scale_y_continuous(limit=c(0,100),breaks = seq(0, 100, by = 10)) +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))

figure4D <- ggplot(my_sum_comp,
                   aes(x=as.factor(treatment),
                       y=mean_croton,
                       color = as.factor(elevation.2),
                       shape=shabitat)) +
  geom_errorbar(aes(ymin=mean_croton-se_croton, ymax=mean_croton+se_croton),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  xlab("Litter composition") +
  scale_x_discrete(limit = c("nomixC","mixB","mixC", "mixD"),
                   labels = c("C","AC","PC", "APC"),
                   expand = expansion(add = c(0.6)))+
  ylab(expression(paste("Litter mass loss of"," ",italic ("C. magdalenensis")," ","(%)")))+
  scale_y_continuous(limit=c(0,100),breaks = seq(0, 100, by = 10)) +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))


tiff("Figure3final.tiff",
     height = 9,
     width = 7,
     unit="in",
     res = 600)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend<-g_legend(figure4A+theme(legend.position="bottom"))


figure4 <- grid.arrange(arrangeGrob(figure4A+
                                      ggtitle("a") +
                                      theme(legend.position = "none",
                                            legend.title=element_text(size=rel(1)),
                                            legend.text=element_text(size=rel(1)),
                                            axis.title.x=element_blank(),
                                            axis.text.x=element_text(size=rel(1),family="Arial"),
                                            axis.ticks.x=element_blank(),
                                            axis.text.y=element_text(size=rel(1)),
                                            plot.title = element_text(size=rel(1),face = "bold")),
                                    figure4B+
                                      ggtitle("b")+
                                      theme(legend.position = "none",
                                            legend.title=element_text(size=rel(1)),
                                            legend.text=element_text(size=rel(1)),
                                            axis.title.x=element_blank(),
                                            axis.text.x=element_text(size=rel(1),family="Arial"),
                                            axis.ticks.x=element_blank(),
                                            axis.text.y=element_text(size=rel(1)),
                                            plot.title = element_text(size=rel(1),face = "bold")),
                                    figure4C+
                                      ggtitle("c")+
                                      theme(legend.position = "none",
                                            axis.title.x=element_text(size=rel(1)),
                                            axis.text.x=element_text(size=rel(1)),
                                            axis.text.y=element_text(size=rel(1)),
                                            plot.title = element_text(size=rel(1),face = "bold")),
                                    figure4D+
                                      ggtitle("d")+
                                      theme(legend.position = "none",
                                            axis.title.x=element_text(size=rel(1)),
                                            axis.text.x=element_text(size=rel(1)),
                                            axis.text.y=element_text(size=rel(1)),
                                            plot.title = element_text(size=rel(1),face = "bold")),ncol = 2),
                        mylegend,heights=c(10, 1))

dev.off()
###litter composicion and environmentla conditions

my_sum_env1<- macroscale%>%
  group_by(treatment,elevation,habitat) %>%
  summarise( 
    n=n(),
    mean_T =mean(temperature  ,na.rm=T),
    sd_T =sd(temperature  ,na.rm = T),
    mean_PH=mean(pH ,na.rm=T),
    sd_PH=sd(pH ,na.rm=T),
    mean_C=mean(conductivity,na.rm=T),
    sd_C=sd(conductivity,na.rm=T)) %>%
  mutate( se_T=sd_T/sqrt(n))  %>%
  mutate( se_PH=sd_PH/sqrt(n))  %>%
  mutate( se_C=sd_C/sqrt(n))  %>%
  mutate(elevation.2=
           case_when(elevation =="High" ~ 3,
                     elevation =="Low" ~ 1,
                     elevation =="Middle" ~ 2))

figureSA1 <- ggplot(my_sum_env1,
                    aes(x=as.factor(treatment),
                        y=mean_T,
                        color = as.factor(elevation.2),
                        shape=habitat)) +
  geom_errorbar(aes(ymin=mean_T-se_T, ymax=mean_T+se_T),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  scale_x_discrete(limit = c("nomixA","nomixB","nomixC", "mixA","mixB","mixC", "mixD"),
                   labels = c("A","P","C","AP","AC","PC", "APC"),
                   expand = expansion(add = c(0.6)))+
  ylab(expression(paste("Water temperature"," ","("^o,"C",")")))+
  scale_y_continuous(limits=c(10,17),breaks = seq(0, 17, by = 1)) +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))

figureSB1 <- ggplot(my_sum_env1,
                    aes(x=as.factor(treatment),
                        y=mean_PH,
                        color = as.factor(elevation.2),
                        shape=habitat)) +
  geom_errorbar(aes(ymin=mean_PH-se_PH, ymax=mean_PH+se_PH),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  scale_x_discrete(limit = c("nomixA","nomixB","nomixC", "mixA","mixB","mixC", "mixD"),
                   labels = c("A","P","C","AP","AC","PC", "APC"),
                   expand = expansion(add = c(0.6)))+
  ylab(expression(paste("Water pH")))+
  scale_y_continuous(limits=c(5,7),breaks = seq(5, 7, by = 0.5)) +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))

figureSC1 <- ggplot(my_sum_env1,
                    aes(x=as.factor(treatment),
                        y=mean_C,
                        color = as.factor(elevation.2),
                        shape=habitat)) +
  geom_errorbar(aes(ymin=mean_C-se_C, ymax=mean_C+se_C),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  xlab("Litter composition") +
  scale_x_discrete(limit = c("nomixA","nomixB","nomixC", "mixA","mixB","mixC", "mixD"),
                   labels = c("A","P","C","AP","AC","PC", "APC"),
                   expand = expansion(add = c(0.6)))+
  ylab(expression(paste("Water total dissolved solid "," ", "(","mgl"^"-1",")"))) +
  scale_y_continuous(limits=c(0,900),breaks = seq(0, 900, by = 100)) +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))



tiff("FigureS2.tiff",
     height = 9,
     width = 7,
     unit="in",
     res = 1200)

grid.arrange(arrangeGrob(figureSA1+
                           ggtitle("a") +
                           theme(legend.position = "none",
                                 legend.title=element_text(size=rel(1)),
                                 legend.text=element_text(size=rel(1)),
                                 axis.title.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.text.y=element_text(size=rel(1)),
                                 plot.title = element_text(size=rel(1))),
                         figureSB1+
                           ggtitle("b")+
                           theme(legend.position = "none",
                                 legend.title=element_text(size=rel(1)),
                                 legend.text=element_text(size=rel(1)),
                                 axis.title.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.text.y=element_text(size=rel(1)),
                                 plot.title = element_text(size=rel(1))),
                         figureSC1+
                           ggtitle("c")+
                           theme(legend.position = "none",
                                 axis.title.x=element_text(size=rel(1),family="Arial"),
                                 axis.text.x=element_text(size=rel(1),family="Arial"),
                                 axis.text.y=element_text(size=rel(1)),
                                 plot.title = element_text(size=rel(1),family="Arial")),ncol = 1),
             mylegend,heights=c(10, 1))
dev.off()

##litter composition and net biodiversity effect, complmentarity effect and selection effect 


my_sum_bio1<- result.diversity2%>%
  group_by(treatment,elevation,habitat) %>%
  summarise( 
    n=n(),
    mean_CE =mean(CE ,na.rm=T),
    sd_CE =sd(CE ,na.rm = T),
    mean_SE=mean(SE,na.rm=T),
    sd_SE=sd(SE,na.rm=T),
    mean_NBE=mean(NBE,na.rm=T),
    sd_NBE=sd(NBE,na.rm=T)
  ) %>%
  mutate( se_CE=sd_CE/sqrt(n))  %>%
  mutate( se_SE=sd_SE/sqrt(n))  %>%
  mutate( ic_SE=se_SE * qt((1-0.05)/2 + .5, n-1))%>%
  mutate( se_NBE=sd_NBE/sqrt(n))  %>%
  mutate( ic_NBE=se_NBE * qt((1-0.05)/2 + .5, n-1))%>%
  mutate(elevation.2=
           case_when(elevation =="High" ~ 3,
                     elevation =="Low" ~ 1,
                     elevation =="Middle" ~ 2))

figure5A <- ggplot(my_sum_bio1,
                   aes(x=as.factor(treatment),
                       y=mean_NBE,
                       color = as.factor(elevation.2),
                       shape=habitat)) +
  geom_errorbar(aes(ymin=mean_NBE-se_NBE, ymax=mean_NBE+se_NBE),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  scale_x_discrete(limit = c( "mixA","mixB","mixC", "mixD"),
                   labels = c("AP","AC","PC", "APC"),
                   expand = expansion(add = c(0.6)))+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("Net biodiversity effect") +
  scale_y_continuous(limit=c(-1,1),breaks = seq(-1, 1, by = 1))  +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))

figure5B <- ggplot(my_sum_bio1,
                   aes(x=as.factor(treatment),
                       y=mean_CE,
                       color = as.factor(elevation.2),
                       shape=habitat)) +
  geom_errorbar(aes(ymin=mean_CE-se_CE, ymax=mean_CE+se_CE),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  scale_x_discrete(limit = c( "mixA","mixB","mixC", "mixD"),
                   labels = c("AP","AC","PC", "APC"),
                   expand = expansion(add = c(0.6)))+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("Complementarity effect") +
  scale_y_continuous(limit=c(-1,1.1),breaks = seq(-1, 1.1, by = 1))  +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))

figure5C <- ggplot(my_sum_bio1,
                   aes(x=as.factor(treatment),
                       y=mean_SE,
                       color = as.factor(elevation.2),
                       shape=habitat)) +
  geom_errorbar(aes(ymin=mean_SE-se_SE, ymax=mean_SE+se_SE),
                width=0.1,
                lwd = 0.6,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5),
             lwd =2) + ggtitle("") +
  geom_hline(yintercept = 0,linetype="dashed")+
  xlab("Litter composition") +
  ylab("Selection effect") +
  scale_x_discrete(limit = c( "mixA","mixB","mixC", "mixD"),
                   labels = c("AP","AC","PC", "APC"),
                   expand = expansion(add = c(0.6)))+
  scale_y_continuous(limit=c(-0.5,0.5),breaks = seq(-0.5, 0.5, by = 1))  +
  scale_shape_manual(name = "  ",
                     labels = c("edge", "interior"),
                     values = c(16,17)) + 
  scale_color_manual(name = "Habitat",
                     labels = c("Low","Mid","High"),
                     values = c("black","#479734","#F0A00E"))+
  theme(legend.position = c(0.9,0.9),
        legend.text=element_text(size=12,family = "Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(family="Arial",size=12),
        axis.text.y=element_text(size=12,family = "Arial"),
        axis.text.x=element_text(size=12,family = "Arial"))


tiff("Figure5.tiff",
     height = 9,
     width = 7,
     unit="in",
     res = 1200)

grid.arrange(arrangeGrob(figure5A+
                           ggtitle("a") +
                           theme(legend.position = "none",
                                 legend.title=element_text(size=rel(1)),
                                 legend.text=element_text(size=rel(1)),
                                 axis.title.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.text.y=element_text(size=rel(1)),
                                 plot.title = element_text(size=rel(1))),
                         figure5B+
                           ggtitle("b")+
                           theme(legend.position = "none",
                                 legend.title=element_text(size=rel(1)),
                                 legend.text=element_text(size=rel(1)),
                                 axis.title.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.text.y=element_text(size=rel(1)),
                                 plot.title = element_text(size=rel(1))),
                         figure5C+
                           ggtitle("c")+
                           theme(legend.position = "none",
                                 axis.title.x=element_text(size=rel(1),family="Arial"),
                                 axis.text.x=element_text(size=rel(1),family="Arial"),
                                 axis.text.y=element_text(size=rel(1)),
                                 plot.title = element_text(size=rel(1),family="Arial")),ncol = 1),
             mylegend,heights=c(10, 1))

dev.off()

