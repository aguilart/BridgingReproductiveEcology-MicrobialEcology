#############################################################################################################
#######                  Data assembly and data analysis described in the manuscript:                 #######
######   "Bridging reproductive and microbial ecology: A case study in arbuscular mycorrhizal fungi"#########
#############################################################################################################

#loading packages

library(tidyverse)
library(scales)
library(phytools)
library(smatr)
library(picante)

####################################################################################

#1- Loading files:

  #1.1. AMF spore data. The species names are based on the taxonomy available in Arthur Schussler website
  #(as stated in the manuscript). Data on diameter dimension (i.e. ranges) was manually digitized from
  #original descriptions of the species. Volumes were calculated assuming sphere shape for globose spores,
  #or prolate spheroids (like an "american football") shape for spores with two distinct diameter ranges:
  
  AMF_All_Copy<-read.csv("Aguilar_etal_SuppMatt_AMF_Spore_Database.csv",header = T, stringsAsFactors = F)
  
  #1.2. Soil ascomycete conidia (spore) data. The taxonomy and diameter range values were manually digitized
  #from the book"Compendium of Soil Fungi" by Domsch, Gams and Anderson. Volumes were calculated using 
  #similar logic as for AMF spores.

  ConidiaDataAscos<-read.csv("Aguilar_etal_SuppMat_CompSoilFungData.csv",header = T, stringsAsFactors = F)
  
  #1.3. Bird egg data. Data was obtained from the supplementary material of Stoddard etal 2017 Science:
  #"Avian egg shape: Form, function, and evolution". Egg volume was estimated based on Baker equation as 
  #used in Stoddard etal 2017 Science. Look at reference 11 in their supplementary material on how to 
  #apply the equation based on Ellipticity and Asymmetry values.

  BirdEgg<-read.csv("Aguilar_etal_SuppMatt_Stoddard_etal_EggData.csv",header = T,stringsAsFactors = F)
            BirdEgg$T_value<-1/(BirdEgg$Ellipticity+1)
            BirdEgg$Lambda<-BirdEgg$Asymmetry+1
            step_1<-sapply(list(-0.75,-0.5,-0.25,0,0.25,0.5,0.75),
                           function(x){((BirdEgg$T_value*(1+x)^(1/(1+BirdEgg$Lambda)))*
                                          ((1-x)^(BirdEgg$Lambda/(1+BirdEgg$Lambda))))^2})
            step_1[,1]<-step_1[,1]*4
            step_1[,2]<-step_1[,2]*2
            step_1[,3]<-step_1[,3]*4
            step_1[,4]<-step_1[,4]*2
            step_1[,5]<-step_1[,5]*4
            step_1[,6]<-step_1[,6]*2
            step_1[,7]<-step_1[,7]*4
            BirdEgg$Volume<-(pi*(BirdEgg$AvgLength..cm.^3)*rowSums(step_1))/96
            
  #1.4. Seed size data. Seed mass was retrieved from the database of the Kew Botanical Gardens. Data is given
  #in miligrams
  
  kewData<-read.csv("Aguilar_etal_SuppMat_SeedSize.csv",header = T,stringsAsFactors = F)
  
  #1.5. AMF life history traits. The data correpond to the series of paper of Hart and Reader on greehouse 
  #experiments with 14 AMF species (see manuscript for full reference details).
  
  Hart_Reader2002_Data<-read.csv("Aguilar_etal_SuppMatt_Hart_Reader2002_Data.csv",header = T,stringsAsFactors = F)
                
######################################################################################################             

#Testing the correlation between spore volume calculated either from inner ranges or outer ranges           

                          ###############
                          ###Figure S1###
                          ###############

AMF_All_Copy%>%
  ggplot()+
  aes(x=SporeVolume,y=SporeVolume_2,label=good.names)+
  geom_point()+
  geom_text(size=1.8)+
  labs(x="Spore volume from inner ranges",
       y="Spore volume from outer ranges")+
  scale_x_continuous(breaks=c(10^7.5,10^8),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_continuous(breaks=c(10^7.5,10^8,10^8.2),
                     labels = trans_format("log10", math_format(10^.x)))+
  geom_smooth(method="lm",se=F)+
  ggtitle(label="Spore volume from outer and inner ranges",
          subtitle = "r2 = 0.71  (p << 0.001)" )+
  theme(title = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size=10))

summary(
  lm(SporeVolume_2~SporeVolume,data=AMF_All_Copy)
)



#####################################################################################################
####Descrition of spore size and variation within the Glomeromycota and comparison to other taxa#####
#####################################################################################################




                                      ###############
                                      ###Figure 1a###
                                      ###############

#This figure correspond to Fig 1a in the manucript
                  rbind(AMF_All_Copy%>%
                          select(good.names,SporeVolume)%>%
                          rename(offsrpingSize=SporeVolume)%>%
                          filter(good.names!="Glomus_tenue")%>%
                          mutate(Taxa="Glomeromycotina")%>%
                          mutate(offsrpingSize=offsrpingSize/((10000)^3))%>%
                          mutate(Organism="Microorganism"),
                        
                        ConidiaDataAscos[grep("conidia$",ConidiaDataAscos$SporeName),]%>%
                          filter(Phylum=="Ascomycota")%>%
                          select(good.names,SporeVolume)%>%
                          rename(offsrpingSize=SporeVolume)%>%
                          mutate(Taxa="Ascomycota")%>%
                          mutate(offsrpingSize=offsrpingSize/((10000)^3))%>%
                          mutate(Organism="Microorganism"),
                        kewData%>%
                          select(species,value)%>%
                          rename(good.names=species, offsrpingSize=value)%>%
                          mutate(Taxa="Angiosperms")%>%
                          mutate(Organism="Macroorganism"),
                        BirdEgg%>%
                          select(Species,Volume)%>%
                          rename(good.names=Species,offsrpingSize=Volume)%>%
                          mutate(Taxa="Birds")%>%
                          mutate(Organism="Macroorganism"))%>%
                    
                    ggplot()+
                    aes(Taxa,offsrpingSize,fill=Taxa)+
                    geom_jitter(size=0.5, width = 0.3,alpha=0.2)+
                    geom_violin(alpha=0.8, draw_quantiles=c(0.25, 0.5, 0.75))+
                    facet_grid(. ~ Organism, scales = "free")+
                    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))+
                    labs(x="Taxa",y=expression("Offspring size"))+
                    theme(title = element_text(size = 25),
                          axis.title.x=element_blank(),
                          axis.text.x = element_text(size = 15),
                          axis.text.y = element_text(size = 20),
                          strip.text.x = element_text(size = 20),
                          legend.position = "none")
              

                  

###################################################################################################
##############Phylogenetic conservatism throughout the entire Glomeromycota########################
###################################################################################################

#Loading the tree
AMF_Tree<-read.tree("tree.for.carlos.tre")

#Selecting trait (spore volume) values of species in the tree
#The phylogenetic tree includes half of the entire spore database
                  
SporeVolumeAll<-AMF_All_Copy[AMF_All_Copy$good.names%in%AMF_Tree$tip.label,18]
                names(SporeVolumeAll)<-AMF_All_Copy[
                  AMF_All_Copy$good.names%in%AMF_Tree$tip.label,1]


#Testing for phylogenetic conservatism                
phylosig(AMF_Tree,SporeVolumeAll,method = "lambda",test = T)
                  
#We subset spore data for species present in the tree for making figure 2b
#(The phylogenetic tree includes half of the entire spore database).

AMF_SubsetInPhyloy<-AMF_All_Copy[
                        AMF_All_Copy$good.names%in%AMF_Tree$tip.label,];
                      rownames(AMF_SubsetInPhyloy)<-NULL
                      AMF_SubsetInPhyloy[                    
                        grep("Scutellospora",AMF_SubsetInPhyloy$good.names),4]<-"Gigasporaceae"
                      AMF_SubsetInPhyloy[                    
                        grep("Cetraspora",AMF_SubsetInPhyloy$good.names),4]<-"Gigasporaceae"
                      AMF_SubsetInPhyloy[                    
                        grep("Racocetra",AMF_SubsetInPhyloy$good.names),4]<-"Gigasporaceae"
                      AMF_SubsetInPhyloy[                    
                        grep("Acaulospora",AMF_SubsetInPhyloy$good.names),4]<-"Acaulosporaceae"
                      AMF_SubsetInPhyloy[                    
                        grep("Glomus",AMF_SubsetInPhyloy$good.names),4]<-"Glomeraceae"
                      AMF_SubsetInPhyloy[                    
                        grep("Rhizophagus",AMF_SubsetInPhyloy$good.names),4]<-"Glomeraceae"
                      
                      AMF_SubsetInPhyloy$Family<-as.character(AMF_SubsetInPhyloy$Family)
                      AMF_SubsetInPhyloy$Family<-factor(AMF_SubsetInPhyloy$Family,
                                                        levels = c("Archaeosporaceae",
                                                                   "Geosiphonaceae",
                                                                   "Ambisporaceae",
                                                                   "Claroideoglomeraceae",
                                                                   "Glomeraceae",
                                                                   "Pacisporaceae",
                                                                   "Gigasporaceae",
                                                                   "Diversisporaceae",
                                                                   "Acaulosporaceae",
                                                                   "Paraglomeraceae",
                                                                   "",
                                                                   "sequences cluster inDiversisporaceae",
                                                                   "Entrophosporaceae"))



                                      #############
                                      ##Figure 2b##
                                      #############
AMF_SubsetInPhyloy%>%
  filter(Family!="")%>%
  filter(Family!="Entrophosporaceae")%>%
  filter(Family!="sequences cluster inDiversisporaceae")%>%
  ggplot()+
  aes(x=Family,y=SporeVolume,fill=Family)+
  geom_jitter(alpha=0.5)+
  geom_violin(alpha=0.5)+
  labs(x= "Family", y=expression("Spore size"~"("*mu*m^3*")"))+
  scale_y_log10(breaks=c(10^4,10^5,10^6,10^7),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+  
  theme(title = element_text(size = 25),
    axis.text.x = element_text(angle=0,vjust = 0.7,size=20),
        axis.text.y = element_text(size = 20),
        axis.line = element_line(color = "black"),
        #strip.text.x = element_text(size = 20),
        legend.position = "none"#,panel.background = element_blank()
    )+
        coord_flip()
  
#########################################################################
###Correlations spore size with other traits#############################
#########################################################################

#Conversion factors of different units for mycelia and spores to biomass

LengthToBiomass_Cf<-((4*10^-6)^2)*pi*((100)^3)*1.09*.21*((10)^6)
#Based on the assumption that hyphae are perfect cylinders of 4um 
#radius and conversion factors calculated on Bakken and Olsen (1983)
#This conversion factors transform Length units (in meters) to micrograms

SporeVolumeToBiomass_Cf<-1.5/4113430 #Based on Beilby and Kidby (1980)

#Data preparation. This data correspond to the ones published in Hart and Reader 2002
#and Hart and Reader 2004 as mentioned in the manuscript. 

                  Hart_Reader2002_Data<-left_join(Hart_Reader2002_Data,AMF_All_Copy[,c(1,18)])
                  
                  Hart_Reader2002_Data$root_fungalBiomass<-((Hart_Reader2002_Data$root_ergo+0.18)/0.4)*LengthToBiomass_Cf
                  Hart_Reader2002_Data$soil_fungalBiomass<-Hart_Reader2002_Data$hyphlength*LengthToBiomass_Cf
                  
                  Hart_Reader2002_Data$total_fungalBiomass<-Hart_Reader2002_Data$root_fungalBiomass+
                  Hart_Reader2002_Data$soil_fungalBiomass

#Creating a summary of this data, with just the means of each variable for each fungus per host.
Hart_Reader2002_Summary<-aggregate(Hart_Reader2002_Data[5:12],
                             by=list(Hart_Reader2002_Data$host,Hart_Reader2002_Data$good.names),
                             mean,na.rm=T);
                      names(Hart_Reader2002_Summary)[1]<-"host";
                      names(Hart_Reader2002_Summary)[2]<-"good.names"

                    Hart_Reader2002_Summary$SporeMass<-Hart_Reader2002_Summary$SporeVolume*
                    SporeVolumeToBiomass_Cf

                    Hart_Reader2002_Summary$no_spore_FungMass<-Hart_Reader2002_Summary$no_spores/
                    Hart_Reader2002_Summary$total_fungalBiomass
                    
                    Hart_Reader2002_Summary$no_spore_FungLength<-Hart_Reader2002_Summary$no_spores/
                    Hart_Reader2002_Summary$hyphlength
                    
                    Hart_Reader2002_Summary$Total_MassForSpores<-Hart_Reader2002_Summary$no_spores*
                    Hart_Reader2002_Summary$SporeMass

######################################################################################                  
#####  SPORE OUTPUT~SPORE MASS #######################################################
######################################################################################

#First correlation, Spore output and spore mass. We expect a negative correlation indicating that, under
# a given amount of resources, the fungus has to decide in between producing a lot of small spores or few
# large one. Testing this relationship in the log-log scales also indicate if the total amount of resources
# allocated to reproduction is constant, that is regardless which strategy, the fungus allocate equal amount
# of resources in spore production


#First a simple correlation between the two, no corrections

#############
##Figure 3a##
#############
#This one correspond to Figure 3a in the Manuscript
Hart_Reader2002_Summary%>%
          ggplot()+
          aes(x=SporeMass,y=no_spores,
              color=host,size=2)+
          geom_point()+
          scale_y_log10(breaks=c(10^0.5,10^1,10^1.5),
                        labels = scales::trans_format("log10", scales::math_format(10^.x)))+
          scale_x_log10(
            labels = scales::trans_format("log10", scales::math_format(10^.x)))+
          #geom_text(aes(label=good.names))+
          labs(x=expression("Spore mass"~"("*mu*g*")"),
               y=expression(Number~spores~fungus^{-1}))+
          ggtitle(label="Spore output vs. spore size",
                  subtitle = "slope= -0.6 (C.I= -0.81,-0.50), r2 = 0.13  (p = 0.004)" )+
          geom_abline(slope = -0.64,intercept = 1.07)+
          #geom_abline(slope = -0.23,intercept = 1.07,lty=2)+
          geom_abline(slope=-1,intercept = 1.07,lty=2)+
          theme(title = element_text(size = 25),
                axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 20),
                legend.text = element_text(size=10),
                legend.position = "none")
          # scale_y_continuous(limits = c(0,2.5))+
          # scale_x_continuous(limits = c(-1.5,1.5))


#Normal regression
confint(
#summary.lm(
lm(log10(no_spores)~log10(SporeMass),data = 
     #Hart_Reader2002_Data))
     Hart_Reader2002_Summary))

#Standarized major axis series of analysis:
# 1. Testing wehther host influces the slope of the relationship & wehther the relationship
# is inversely proportinal (testing wehther slope = -1)
summary(
sma(no_spores ~ SporeMass*host,
    data = Hart_Reader2002_Summary, log = "xy",slope.test = -1))
# 2. Testing whether the host influences the intercept (aka elevation of the relationship) 
summary(
sma(no_spores ~ SporeMass+host,
    data = Hart_Reader2002_Summary, log = "xy",slope.test = -1))

#Results: 1) Plant host does not alter the slope of the relationship between 
#           spore production and spore size, nor the intercept (FYI, in this package
#           intercept is called elevation)
#         2) The slope is different from -1, meaning the relationship is not 
#           inversely proportional. 
# Based on this, I decided to ignore host as a driver, however I decided not just average
# the spore output for each host, as they are real replicates. 
# Thus, the model in the end looks like:
summary(
sma(no_spores ~ SporeMass,
    data = Hart_Reader2002_Summary, log = "xy",slope.test = -1))
# Which returns the following parameters:
# Coefficients:
#             elevation      slope
# estimate     1.038594 -0.6412302
# lower limit  0.921332 -0.8170478
# upper limit  1.155856 -0.5032462
# 
# H0 : variables uncorrelated
# R-squared : 0.1333479 
# P-value : 0.0041194 



# Second. 

#A) Correcting for total fungal size (measured in hyphal length in soil)

###############
###Figure 3b###      
###############
Hart_Reader2002_Summary%>%
  ggplot()+
  aes(x=SporeMass,
      y=(no_spores/hyphlength),color=host,size=2)+
  geom_point()+
  scale_y_log10(#breaks=c(10^-1.5,10^-1,10^-0.5),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_x_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  #geom_text(aes(label=good.names))+
  labs(x=expression("Spore mass"~"("*mu*g*")"),
       y=expression(Number~spores~(meters.extradical.mycelia)^{-1}))+
  ggtitle(label="Spore output per meter of extraradical mycelia vs. spore size",
          subtitle = "slope= -0.91 (CI=0.44,0.73), r2= 0.2 (p=0.00002)"
  )+
  geom_abline(slope = -.91,intercept = 0.59)+
  #geom_abline(slope = -0.23,intercept = 1.07,lty=2)+
  geom_abline(slope=-1,intercept = 0.58,lty=2)+
  theme(title = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size=10),
        legend.position = "none"
  )

#Since, as usual host identity does not play a role, this is the analysis I used:  
summary(
  sma(no_spore_FungLength ~ SporeMass,
      data = Hart_Reader2002_Summary, log = "xy",slope.test=-1))

#B) Correcting for total fungal size (measured in biomass of fungi)

###############
###Figure 3c###      
###############

#This one correspond to figure 3c in the manuscript
Hart_Reader2002_Summary%>%
ggplot()+
aes(x=SporeMass,
    y=(no_spores/total_fungalBiomass),color=host,size=2)+
geom_point()+
scale_y_log10(breaks=c(10^-1.5,10^-1,10^-0.5),
              labels = scales::trans_format("log10", scales::math_format(10^.x)))+
scale_x_log10(
  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
#geom_text(aes(label=good.names))+
labs(x=expression("Spore mass"~"("*mu*g*")"),
     y=expression(Number~spores~(gram.fungus)^{-1}))+
ggtitle(label="Spore output per gram of fungus vs. spore size",
        subtitle = "slope= -0.59 (CI=-0.75,-0.47), r2= 0.2 (p=0.0002)")+
geom_abline(slope = -0.59,intercept = -0.94)+
#geom_abline(slope = -0.23,intercept = 1.07,lty=2)+
geom_abline(slope=-1,intercept = -0.94,lty=2)+
theme(title = element_text(size = 25),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.text = element_text(size=10)#,
      #legend.position = "none"
)

#normal regression
confint(
#summary.lm(
lm(log10(no_spore_FungMass)~
     log10(SporeMass),
   data = 
     Hart_Reader2002_Summary))
#Hart_Reader2002_Data))

#Standardized major axis
summary(
sma(no_spore_FungMass ~ SporeMass+host,
    data = Hart_Reader2002_Summary, log = "xy",slope.test=-1))
#Results from analysis:
# Similar when total fungal biomass is not taken into account (see above),
# - Host has not influence neither in the slope nor the intercept of the relationship
# - The relationship is different from -1 (spore output per gram not inversely proportional
# - to spore size)
# - Thus, the anlysis used is: 
summary(
sma(no_spore_FungMass ~ SporeMass,
    data = Hart_Reader2002_Summary, log = "xy",slope.test = -1))
# -Resulting in the following parameters:
# Coefficients:
#             elevation      slope
# estimate    -0.9436634 -0.5969789
# lower limit -1.0447225 -0.7526987
# upper limit -0.8426042 -0.4734746
# 
# H0 : variables uncorrelated
# R-squared : 0.2083189 
# P-value : 0.00024703 

######################################################################################                  
#####  TOTAL AMOUNT OF RESOURCES FOR SPORE PRODUCTION ~ SPORE MASS ###################
######################################################################################  

###############
###FIGURE 3c###
###############

#This one correspond to figure 3c in the manuscript
Hart_Reader2002_Summary%>%
ggplot()+
aes(x=SporeMass,
    y=Total_MassForSpores,
    color=host,size=2)+
geom_point()+
scale_y_log10(breaks=c(10^0,10^0.5,10^1,10^1.5,10^2),
              labels = scales::trans_format("log10", scales::math_format(10^.x)))+
scale_x_log10(
  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
#geom_text(aes(label=good.names))+
labs(x=expression("Spore mass"~"("*mu*g*")"),
     y= expression(Mass~allocated~to~spore~production~"("*mu*g*")"))+
ggtitle(label = "Mass allocated to reproduction vs. spore size",
        subtitle ="slope= 0.97 (CI=0.82,1.14), r2=0.62 (p<<0.001)" )+
#geom_abline(slope=0.76, intercept = 1.07)+
geom_abline(slope = 0.97,intercept = 1.08,lty=2)+
theme(title = element_text(size = 25),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.position =  "none")

#Regression analysis
confint(
#summary.lm(
lm(log10(Total_MassForSpores)~
     log10(SporeMass),
   data = 
     Hart_Reader2002_Summary))
#Hart_Reader2002_Data))

#Standardized major axis
summary(
sma(Total_MassForSpores ~ SporeMass,
    data = Hart_Reader2002_Summary, log = "xy"))

#Results:
# Coefficients:
#            elevation     slope
# estimate    1.0882362 0.9710110
# lower limit 0.9856066 0.8267049
# upper limit 1.1908657 1.1405065
# 
# H0 : variables uncorrelated
# R-squared : 0.6220586 
# P-value : 7.2854e-14 


#Results from this analysis indicate that
#1. As expected based on the previous analysis, host has no influence
#   neither in the slope nor the intercept of this relationship
#2. There is a positive correlation between spore size and total amount of resources
#   allocated to spore production. 
#3. The slope is "quite high" still the slope of the regression is shallower (and congruent
#   with the slope obtained in the analysis no_spores~sporeMass). For some weird reason the 
#   the slope of the SMA is not congruent with the analyis of no_spores~sporeMass.
# The results are:
#Coefficients:
#             elevation     slope
# estimate    1.0882362 0.9710110
# lower limit 0.9856066 0.8267049
# upper limit 1.1908657 1.1405065
# 
# H0 : variables uncorrelated
# R-squared : 0.6220586 
# P-value : 7.2854e-14 


######################################################################################                  
#####  SPORE OUTPUT ~ FUNGAL BIOMASS  ################################################
######################################################################################

###############
###FIGURE 4a###
###############

#This one correspond to figure 4a in the manuscript 
Hart_Reader2002_Summary%>%
ggplot()+
aes(y=no_spores,
    x=total_fungalBiomass,
    color=host,size=2)+
geom_point()+
scale_y_log10(breaks=c(10^0.5,10^1,10^1.5),
              labels = scales::trans_format("log10", scales::math_format(10^.x)))+
scale_x_log10(breaks=c(10^1.7,10^1.8,10^1.9,10^2.0,10^2.1),
              labels = scales::trans_format("log10", scales::math_format(10^.x)))+
labs(x=expression("Total fungal mass"~"("*mu*g*")"),
     y=expression(Number~spores~fungus^{-1}))+
ggtitle(label="Spore output vs Fungal biomass",
        subtitle = "r2= 0.13 (p=0.0035)")+
theme(title = element_text(size = 25),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.position = "none")

summary.lm(
lm(log10(no_spores)~
     log10(total_fungalBiomass),
   data = Hart_Reader2002_Summary))

summary(
sma(no_spores ~ total_fungalBiomass,
    data = Hart_Reader2002_Summary, log = "xy"))


#Results:
# Coefficients:
#            elevation    slope
# estimate    -5.355673 3.246103
# lower limit -6.933817 2.548910
# upper limit -3.777528 4.133997
# 
# H0 : variables uncorrelated
# R-squared : 0.1371333 
# P-value : 0.0035864 



######################################################################################                  
#####  TOTAL FUNGAL BIOMASS ~ SPORE MASS  ############################################
######################################################################################

###############
###FIGURE 4b###
###############

#This one correspond to Figure 4b in the manuscript
Hart_Reader2002_Summary%>%
ggplot()+
aes(x=total_fungalBiomass,
    y=SporeMass,
    color=host,size=2)+
geom_point()+
scale_y_log10(breaks=c(10^-1,10^-0.5,10^0,10^0.5,10^1),
              labels = scales::trans_format("log10", scales::math_format(10^.x)))+
scale_x_log10(breaks=c(10^1.7,10^1.8,10^1.9,10^2.0,10^2.1),
              labels = scales::trans_format("log10", scales::math_format(10^.x)))+
labs(x=expression("Total fungal mass"~"("*mu*g*")"),
     y=expression("Spore mass"~"("*mu*g*")"))+
ggtitle(label="Spore mass vs fungal biomass",
        subtitle = "r2=0.03 (p=0.14)")+
theme(title = element_text(size = 25),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.position =  "none")

summary.lm(
lm(log10(SporeMass)~
     log10(total_fungalBiomass),data = Hart_Reader2002_Summary))

summary(
sma(SporeMass ~ total_fungalBiomass,
    data = Hart_Reader2002_Summary, log = "xy"))

#Results
# Coefficients:
#              elevation    slope
# estimate    -10.134039 5.062306
# lower limit -12.734241 3.922583
# upper limit  -7.533837 6.533179
# 
# H0 : variables uncorrelated
# R-squared : 0.03762316 
# P-value : 0.13754 


##############################################################################################  
###############Phylogenetic corrections ######################################################
##############################################################################################


AMF_Tree_Hart_Reader2002Data<-
  drop.tip(AMF_Tree,
           AMF_Tree$tip.label[!AMF_Tree$tip.label%in%Hart_Reader2002_Data$good.names])

#Creating a vector with the variables of interest
Hart_Reader2002_SummaryBySpecies<-aggregate(Hart_Reader2002_Data[5:12],
                                            by=list(Hart_Reader2002_Data$good.names),
                                            mean,na.rm=T);names(Hart_Reader2002_SummaryBySpecies)[1]<-"good.names"


names(Hart_Reader2002_Summary)[1]<-"host";
names(Hart_Reader2002_Summary)[2]<-"good.names"

Hart_Reader2002_SporeSize<-Hart_Reader2002_SummaryBySpecies[,9];# using volume
names(Hart_Reader2002_SporeSize)<-Hart_Reader2002_SummaryBySpecies[,1]

Hart_Reader2002_SporeOutput<-Hart_Reader2002_SummaryBySpecies[,5];# using no_spores
names(Hart_Reader2002_SporeOutput)<-Hart_Reader2002_SummaryBySpecies[,1]

#testing for phylogenetic conservatism for only the Hart_Reader2002 Data
phylosig(AMF_Tree_Hart_Reader2002Data,Hart_Reader2002_SporeSize,method = "lambda",test = T)
phylosig(AMF_Tree_Hart_Reader2002Data,Hart_Reader2002_SporeOutput,method = "lambda",test = T)
phylosig(AMF_Tree_Hart_Reader2002Data,Hart_Reader2002_SporeSize,method = "K",test = T)

#creating a new vector with the phylogenetic independent contrasts, the
#fungus Septoglomus_constrictum and Dentiscutata_heterogama are not in the
#tree and need to be removed from the dataset.

pic_SporeVol<-pic(log10(Hart_Reader2002_SporeSize[-c(6,15)]),AMF_Tree_Hart_Reader2002Data)
pic_SporeOutput<-pic(log10(Hart_Reader2002_SporeOutput[-c(6,15)]),AMF_Tree_Hart_Reader2002Data)

summary.lm(
lm(pic_SporeOutput~pic_SporeVol))

ggplot()+
aes(x=pic_SporeVol,
    y=pic_SporeOutput)+
geom_point()

#Based on this result one can see that the trade-off is the result of trait covariation pattern that
#resulted from limited evolutionary events that have been retained through lineages in the phylogeny.

