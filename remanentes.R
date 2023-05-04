# Data Analysis: Arsenio Vargas-Vazquez, Omar Yaxmehen Bello-Chavolla
# Latest version of Analysis February-2023
# Any question regarding analysis contact Omar Yaxmehen Bello-Chavolla at oyaxbell@yahoo.com.mx

#### Package management ####
pacman::p_load(readxl,haven, tidyverse, ggpubr, lmtest, nortest, gtools, 
               data.table, caret, glmnet, survival, flextable, blandr, 
               BlandAltmanLeh,see,nephro,skedastic,gtsummary,
               rms, bestNormalize, flexsurv, pROC, timeROC, fmsb, factoextra,
               gridExtra,  nhanesA, wesanderson,forestmodel, ggedit,
               FactoMineR, fpc, NbClust, ggimage, glmnet, ggsci, survminer,
               cluster, ggplotify, UpSetR,tseries, simPH,gvlma,
               nortest, viridis, officer, magrittr, mice,dummies,jtools,
               PerformanceAnalytics, corrplot, mice, performance, SAScii, sjPlot)
setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUT?NOMA DE M?XICO/Remanentes de colesterol - Ars")
#setwd("C:/Users/Investigacion/OneDrive - UNIVERSIDAD NACIONAL AUT?NOMA DE M?XICO/Remanentes de colesterol - Ars")

## COMPUTADORA ARSENIO ##
#setwd("~/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTO??NOMADEME??XICO/Remanentes de colesterol - Ars - OMAR YAXMEHEN BELLO CHAVOLLA's files")

## method=="Aalen" -> output Hazard Difference
## method=="Cox"   -> output Hazard Ratio

mediation_ci1 <- function(lambda.s, lambda.g, covar11, covar12,
                          covar22, alpha.s, var_alpha, G=10^4, method){
  require(mvtnorm)
  Omega <- matrix(c(covar11,covar12,covar12,covar22),nrow=2)
  IE <- rep(0,G); DE <- rep(0,G); TE <- rep(0,G); Q <- rep(0,G)
  
  set.seed(137)
  lambda <- rmvnorm(G, mean = c(lambda.s, lambda.g), sigma = Omega)
  alpha <- rnorm(G, mean=alpha.s, sd=sqrt(var_alpha))
  DE <- lambda[,1]
  IE <- lambda[,2] * alpha
  TE <- IE + DE
  
  DE.obs <- lambda.s
  IE.obs <- lambda.g * alpha.s
  TE.obs <- DE.obs+IE.obs
  pval.DE<-2*min(mean((DE-mean(DE))>DE.obs), mean((DE-mean(DE))<DE.obs))
  pval.IE<-2*min(mean((IE-mean(IE))>IE.obs), mean((IE-mean(IE))<IE.obs))
  pval.TE<-2*min(mean((TE-mean(TE))>TE.obs), mean((TE-mean(TE))<TE.obs))
  
  if (method=="Cox") {DE=exp(DE); IE=exp(IE); TE=exp(TE)}
  print("DE:")
  print(ifelse(method=="Aalen", DE.obs, exp(DE.obs)))
  print(quantile(DE, c(0.025, 0.975)))
  print(paste("pval_DE=", pval.DE))
  print("IE:")
  print(ifelse(method=="Aalen", IE.obs, exp(IE.obs)))
  print(quantile(IE, c(0.025, 0.975)))
  print(paste("pval_IE=", pval.IE))
  print("TE:")
  print(ifelse(method=="Aalen", TE.obs, exp(TE.obs)))
  print(quantile(TE, c(0.025, 0.975)))
  print(paste("pval_TE=", pval.TE))
}

mediation_ci2 <- function(lambda.s, lambda.m, lambda.g, Sigma.lambda,
                          alpha.s, alpha.m, Sigma.alpha, delta.s, Sigma.delta,
                          G=10^4, method){
  require(mvtnorm)
  SY <- rep(0,G); SGY <- rep(0,G); SMY <- rep(0,G); TE <- rep(0,G)
  set.seed(137)
  lambda <- rmvnorm(G, mean = c(lambda.s, lambda.m, lambda.g),
                    sigma = Sigma.lambda)
  alpha <- rmvnorm(G, mean = c(alpha.s, alpha.m),
                   sigma = Sigma.alpha)
  delta <- rnorm(G, mean=delta.s, sd=sqrt(Sigma.delta))
  SY <- lambda[,1]
  SGY <- lambda[,3] * alpha[,1]
  SMY <- (lambda[,2] + lambda[,3]*alpha[,2])*delta
  TE <- SY+SGY+SMY
  
  SY.obs <- lambda.s
  SGY.obs <- lambda.g * alpha.s
  SMY.obs <- (lambda.m + lambda.g*alpha.m)*delta.s
  TE.obs <- SY.obs+SGY.obs+SMY.obs
  pval.SY<-2*min(mean((SY-mean(SY))>SY.obs), mean((SY-mean(SY))<SY.obs))
  pval.SGY<-2*min(mean((SGY-mean(SGY))>SGY.obs), mean((SGY-mean(SGY))<SGY.obs))
  pval.SMY<-2*min(mean((SMY-mean(SMY))>SMY.obs), mean((SMY-mean(SMY))<SMY.obs))
  pval.TE<-2*min(mean((TE-mean(TE))>TE.obs), mean((TE-mean(TE))<TE.obs))
  
  if (method=="Cox") {SY=exp(SY); SGY=exp(SGY); SMY=exp(SMY); TE=exp(TE)}
  print("SY:")
  print(ifelse(method=="Aalen", SY.obs, exp(SY.obs)))
  print(quantile(SY, c(0.025, 0.975)))
  print(paste("pval_SY=", pval.SY))
  print("SGY:")
  print(ifelse(method=="Aalen", SGY.obs, exp(SGY.obs)))
  print(quantile(SGY, c(0.025, 0.975)))
  print(paste("pval_SGY=", pval.SGY))
  print("SMY:")
  print(ifelse(method=="Aalen", SMY.obs, exp(SMY.obs)))
  print(quantile(SMY, c(0.025, 0.975)))
  print(paste("pval_SMY=", pval.SMY))
  print("TE:")
  print(ifelse(method=="Aalen", TE.obs, exp(TE.obs)))
  print(quantile(TE, c(0.025, 0.975)))
  print(paste("pval_TE=", pval.TE))
}

#### Building databases ######

#NHANES-III
NHANES3<-fread("Bases/nhanes3.csv", na = c("", "N/A", "NA","na", "#N/A", "88888", 
                                     "8888", "888888", "888","5555","999","9999",
                                     "99999", "999999", "9998"))[,-c(1,2)]
mortalidad3<-fread("Bases/nhanes3_mort.csv"); mortalidad3$SEQN<-mortalidad3$seqn

#NHANES-IV
NHANES<-fread("Bases/nhanes4.csv",na = c("", "N/A", "NA","na", "#N/A", "88888", "8888", "888888", "888","5555","999","9999", "99999", "999999"))[,-c(1)]
mortalidad<-fread("Bases/nhanes_mortalidad.csv"); mortalidad$SEQN<-mortalidad$seqn

### NHANES-III management ###
NHANES0 <- NHANES3 %>% filter(HSAGEIR.x>=20)%>% 
  mutate("BMI"=BMPBMI, "Diabetes"=HAD1,
         "Hypertension"=HAE2, "Asthma"=HAC1E, "Arthritis"=HAC1A,
         "Heart_failure"=HAC1C, "Heart_attack"=HAF10, "Stroke"=HAC1D, 
         "Emphysema"=HAC1G, "Bronchitis"=HAC1F, "Malignancy"=HAC1O,
         "Sex"=factor(HSSEX.x, levels= c(1,2),labels= c("Men", "Women")), 
         "Age"=HSAGEIR.x, "Ethnicity"=DMARETHN.x, "HDL"=HDP, "TG"=TGP, 
         "ApoB"=ABP, "Cholesterol"=TCP, "Glucose"=SGP, "HbA1c"=GHP, 
         "C_peptide"=C1PSI, no_HDL=TCP-HDP,
         "LDL_S"=(TCP/0.948)-(HDP/0.971)-((TGP/8.56)+((TGP*no_HDL)/2140)-
                                            ((TGP*TGP)/16100))-9.4,
         "Rem_S"=TCP-HDP-LDL_S, "Rem_M"=TCP-HDP-LDL_M,
         "VLDL_C"=((TGP/8.59)+((TGP*no_HDL)/2250)-((TGP*TGP)/16500)), 
         "HOMA2IR"=`HOMA2 IR`, "HOMA2B"=`HOMA2 %B`, "HOMA2S"=`HOMA2 %S`,
         "Fuma_cig"=HAR3,"Fuma_pur"=HAR24, "Fuma_pip"=HAR27) 
NHANES0$Smoking<-ifelse((NHANES0$Fuma_cig==1|NHANES0$Fuma_pur==1|NHANES0$Fuma_pip==1),1,0)
NHANES0$Smoking[is.na(NHANES0$Smoking)]<-0
#Cause-specific mortality
NHANES0 <- merge(NHANES0,mortalidad3,by="SEQN")
d1<-dummies::dummy(NHANES0$ucod_leading)
colnames(d1)<-c("Heart_diseases","Malignant_neoplasms",
                "Chronic_lower_respiratory_diseases","Accidents",
                "Cerebrovascular_diseases","Alzheimer_disease",
                "Diabetes_mellitus","Influenza_or_pneumonia",
                "Nephone_diseases","Other_causes","Alive"); NHANES0<-cbind(NHANES0, d1)
nhanes<-NHANES0 %>% dplyr::select(SEQN, Age, Sex, Ethnicity, mortstat,permth_int,
                                  permth_exm,BMI, HDL, TG, 
                                  Cholesterol, Glucose, HbA1c, C_peptide,
                                  no_HDL, LDL_S, LDL_M,Rem_S, Rem_M,
                                  VLDL_C,HOMA2IR, HOMA2B, HOMA2S,"Heart_diseases",
                                  "Malignant_neoplasms",
                                  "Chronic_lower_respiratory_diseases","Accidents",
                                  "Cerebrovascular_diseases","Alzheimer_disease",
                                  "Diabetes_mellitus","Influenza_or_pneumonia",
                                  "Nephone_diseases","Other_causes","Alive",
                                  Diabetes, Hypertension, Asthma, Arthritis,
                                  Heart_failure, 
                                  Heart_attack, Stroke, Emphysema, Bronchitis,
                                  Malignancy,Smoking); nhanes$id <- rep(1,nrow(nhanes))

#Comorbidities
nhanes$Diabetes<-ifelse(nhanes$Diabetes==1, 1, 0); 
nhanes$Asthma<-ifelse(nhanes$Asthma==1, 1, 0)
nhanes$Arthritis<-ifelse(nhanes$Arthritis==1, 1, 0);
nhanes$Heart_failure<-ifelse(nhanes$Heart_failure==1, 1, 0)
nhanes$Heart_attack<-ifelse(nhanes$Heart_attack==1, 1, 0); 
nhanes$Emphysema<-ifelse(nhanes$Emphysema==1, 1, 0)
nhanes$Bronchitis<-ifelse(nhanes$Bronchitis==1, 1, 0); 
nhanes$Malignancy<-ifelse(nhanes$Malignancy==1, 1, 0)
nhanes$Stroke<-ifelse(nhanes$Stroke==1, 1, 0); 
nhanes$Hypertension<-ifelse(nhanes$Hypertension==1, 1, 0)

#Number of comorbidities
nhanes$num_comorb<-nhanes$Asthma+nhanes$Arthritis+
  nhanes$Heart_failure+nhanes$Heart_attack+nhanes$Emphysema+
  nhanes$Bronchitis+nhanes$Malignancy+nhanes$Stroke

#Race/Ethnicity
nhanes$Ethnicity<-factor(nhanes$Ethnicity, labels = c("White", "Black", 
                                                      "Mexican-American", "Other"))

nhanes$ASCVD<-NULL
nhanes$ASCVD[nhanes$Heart_diseases==1 & nhanes$Cerebrovascular_diseases==1]<-1
nhanes$ASCVD[nhanes$Heart_diseases==0 & nhanes$Cerebrovascular_diseases==1]<-1
nhanes$ASCVD[nhanes$Heart_diseases==1 & nhanes$Cerebrovascular_diseases==0]<-1
nhanes$ASCVD[nhanes$Heart_diseases==0 & nhanes$Cerebrovascular_diseases==0]<-0

nhanes$ASCVD_mort <- NULL
nhanes$ASCVD_mort[nhanes$ASCVD==1 & nhanes$mortstat==1]<-1
nhanes$ASCVD_mort[nhanes$ASCVD==0 & nhanes$mortstat==1]<-2
nhanes$ASCVD_mort[nhanes$ASCVD==0 & nhanes$mortstat==0]<-0
nhanes$ASCVD_mort<-factor(nhanes$ASCVD_mort, labels = c("Censored",
                                                          "ASCVD mortality", "Other causes"))
# Diabetes definition
nhanes$dm_hba1c[nhanes$HbA1c>=6.5]<-1; nhanes$dm_hba1c[nhanes$HbA1c<6.5]<-0
nhanes$dm_G[nhanes$Glucose>=126]<-1;nhanes$dm_G[nhanes$Glucose<126]<-0

# General diabetes variable; DMF#

nhanes$DM<-NULL
nhanes$DM[nhanes$Diabetes==1 & nhanes$dm_hba1c==1]<-1
nhanes$DM[nhanes$Diabetes==0 & nhanes$dm_hba1c==1]<-1
nhanes$DM[nhanes$Diabetes==1 & nhanes$dm_hba1c==0]<-1
nhanes$DM[nhanes$Diabetes==0 & nhanes$dm_hba1c==0]<-0

nhanes$DMF<-NULL
nhanes$DMF[nhanes$DM==1 & nhanes$dm_G==1]<-1
nhanes$DMF[nhanes$DM==0 & nhanes$dm_G==1]<-1
nhanes$DMF[nhanes$DM==1 & nhanes$dm_G==0]<-1
nhanes$DMF[nhanes$DM==0 & nhanes$dm_G==0]<-0

NHANES1 <- NHANES %>% filter(RIDAGEYR>=20) %>% 
  mutate("DXDTOFAT_N" = (DXDTOFAT/1000)/((BMXWT/100)^2), "DXXTRFAT_N" = (DXXTRFAT/1000)/((BMXWT/100)^2), "DXXHEFAT_N" = (DXXHEFAT/1000)/((BMXWT/100)^2),
         "DXXLAFAT_N" = (DXXLAFAT/1000)/((BMXWT/100)^2), "DXXRAFAT_N" = (DXXRAFAT/1000)/((BMXWT/100)^2), "DXXLLFAT_N" = (DXXLLFAT/1000)/((BMXWT/100)^2),
         "DXXRLFAT_N" = (DXXRLFAT/1000)/((BMXWT/100)^2), "DXXTRFAT_DXDTOFAT" = (DXXTRFAT)/(DXDTOFAT), "Weight" = BMXWT, "Waist" = BMXWAIST, "Height" = BMXHT,
         "BMI" = BMXBMI, "Calf_circumference" = BMXCALF, "Arm_circumference" = BMXARMC, "Thigh_circumference" = BMXTHICR, "Triceps_skinfold" = BMXTRI,
         "Subscapular_skinfold" = BMXSUB, "Leg_length" = BMXLEG, "ICE" = BMXWAIST/BMXHT, "Arm_length"= BMXARML, "METSIR" = (log(Glucose*2+TG)*BMXBMI)/log(HDL),
         "EXTSUP_FAT" = DXXLAFAT_N + DXXRAFAT_N, "EXTINF_FAT" = DXXLLFAT_N + DXXRLFAT_N, "DXDHELE_N" = (DXDHELE/1000), "DXDLALE_M" = (DXDLALE/1000), "DXDRALE_N" = (DXDRALE/1000),
         "DXDRLLE_N" = (DXDRLLE/1000), "DXDLALE_N" = (DXDLALE/1000), "DXDLLLE_N" = (DXDLLLE/1000), "DXDTOLE_N" = (DXDTOLE/1000), "DXDTRLE_N" = (DXDTRLE/1000)) %>%    
  mutate("METS_VF" = 4.466 + 0.011*(log(METSIR)^3)+ 3.239*(log(ICE)^3)-0.319*(2-RIAGENDR) + 0.594*(log(RIDAGEYR)), "EXTSUP_FAT_N" = (EXTSUP_FAT/1000)/((Height/100)^2),
         "EXTINF_FAT_N" = (EXTINF_FAT/1000)/((Height/100)^2), "EXTSUP_FAT_DXDTOFAT" = (EXTSUP_FAT/1000)/(DXDTOFAT/1000), "EXTINF_FAT_DXDTOFAT" = (EXTINF_FAT/1000)/(DXDTOFAT/1000),
         "HDL"=HDL, "TG"=TG,"Cholesterol"=TC, "Glucose"=Glucose, "HbA1c"=HbA1c, 
         "C_peptide"=Cpeptide, "LDL_M"=LDLCM,no_HDL=TC-HDL,
         "LDL_S"=(TC/0.948)-(HDL/0.971)-((TG/8.56)+((TG*no_HDL)/2140)-
                                            ((TG*TG)/16100))-9.4,
         "Rem_S"=TC-HDL-LDL_S, "Rem_M"=TC-HDL-LDL_M,
         "VLDL_C"=((TG/8.59)+((TG*no_HDL)/2250)-((TG*TG)/16500)), 
         "HOMA2IR"=HOMA2IR, "HOMA2B"=HOMA2B, "HOMA2S"=HOMA2S, Smoking=SMQ040) %>% 
  mutate("Sex"=factor(RIAGENDR, levels= c(1,2),labels= c("Men", "Women")), "Age"=RIDAGEYR, "Diabetes"=DIQ010, "Hypertension"=BPQ020, "Asthma"=MCQ010,
         "Arthritis"=MCQ160A, "Heart_failure"=MCQ160B, "Heart_attack"=MCQ160E, "Stroke"=MCQ160F, "Emphysema"=MCQ160G, "Bronchitis"=MCQ160K, "Malignancy"=MCQ220)

#Recode race/ethnicity
NHANES1$Ethnicity[NHANES1$RIDRETH1==5] <- 4
NHANES1$Ethnicity[NHANES1$RIDRETH1==1] <- 3; NHANES1$Ethnicity[NHANES1$RIDRETH1==2] <- 4
NHANES1$Ethnicity[NHANES1$RIDRETH1==3] <- 1; NHANES1$Ethnicity[NHANES1$RIDRETH1==4] <- 2

#Smoking
NHANES1$Smoking<-ifelse(NHANES1$Smoking==3,0,1)
NHANES1$Smoking[is.na(NHANES1$Smoking)]<-0

#Cause-specific mortality
NHANES1 <- merge(NHANES1,mortalidad,by="SEQN")
NHANES1$ucod_leading[is.na(NHANES1$ucod_leading)] <- 11
d1<-dummies::dummy(NHANES1$ucod_leading)
colnames(d1)<-c("Heart_diseases","Malignant_neoplasms","Chronic_lower_respiratory_diseases","Accidents",
                "Cerebrovascular_diseases","Alzheimer_disease","Diabetes_mellitus","Influenza_or_pneumonia",
                "Nephone_diseases","Other_causes","Alive"); NHANES1<-cbind(NHANES1,d1)

nhanes1<-NHANES1 %>% dplyr::select(SEQN, Age, Sex, Ethnicity, mortstat,permth_int,
                                   permth_exm,BMI, HDL, TG, 
                                   Cholesterol, Glucose, HbA1c, C_peptide,
                                   no_HDL, LDL_S, LDL_M,Rem_S, Rem_M,
                                   VLDL_C,HOMA2IR, HOMA2B, HOMA2S,"Heart_diseases",
                                   "Malignant_neoplasms",
                                   "Chronic_lower_respiratory_diseases","Accidents",
                                   "Cerebrovascular_diseases","Alzheimer_disease",
                                   "Diabetes_mellitus","Influenza_or_pneumonia",
                                   "Nephone_diseases","Other_causes","Alive",
                                   Diabetes, Hypertension, Asthma, Arthritis,
                                   Heart_failure, 
                                   Heart_attack, Stroke, Emphysema, Bronchitis,
                                   Malignancy, Smoking); nhanes1$id<-rep(2, nrow(nhanes1))
#Comorbidities
nhanes1$Diabetes<-ifelse(nhanes1$Diabetes==1, 1, 0); 
nhanes1$Asthma<-ifelse(nhanes1$Asthma==1, 1, 0)
nhanes1$Arthritis<-ifelse(nhanes1$Arthritis==1, 1, 0);
nhanes1$Heart_failure<-ifelse(nhanes1$Heart_failure==1, 1, 0)
nhanes1$Heart_attack<-ifelse(nhanes1$Heart_attack==1, 1, 0); 
nhanes1$Emphysema<-ifelse(nhanes1$Emphysema==1, 1, 0)
nhanes1$Bronchitis<-ifelse(nhanes1$Bronchitis==1, 1, 0); 
nhanes1$Malignancy<-ifelse(nhanes1$Malignancy==1, 1, 0)
nhanes1$Stroke<-ifelse(nhanes1$Stroke==1, 1, 0); 
nhanes1$Hypertension<-ifelse(nhanes1$Hypertension==1, 1, 0)

#Number of comorbidities
nhanes1$num_comorb<-nhanes1$Asthma+nhanes1$Arthritis+
  nhanes1$Heart_failure+nhanes1$Heart_attack+nhanes1$Emphysema+
  nhanes1$Bronchitis+nhanes1$Malignancy+nhanes1$Stroke

#Race/Ethnicity
nhanes1$Ethnicity<-factor(nhanes1$Ethnicity, labels = c("White", "Black", 
                                                      "Mexican-American", "Other"))

nhanes1$ASCVD<-NULL
nhanes1$ASCVD[nhanes1$Heart_diseases==1 & nhanes1$Cerebrovascular_diseases==1]<-1
nhanes1$ASCVD[nhanes1$Heart_diseases==0 & nhanes1$Cerebrovascular_diseases==1]<-1
nhanes1$ASCVD[nhanes1$Heart_diseases==1 & nhanes1$Cerebrovascular_diseases==0]<-1
nhanes1$ASCVD[nhanes1$Heart_diseases==0 & nhanes1$Cerebrovascular_diseases==0]<-0

nhanes1$ASCVD_mort <- NULL
nhanes1$ASCVD_mort[nhanes1$ASCVD==1 & nhanes1$mortstat==1]<-1
nhanes1$ASCVD_mort[nhanes1$ASCVD==0 & nhanes1$mortstat==1]<-2
nhanes1$ASCVD_mort[nhanes1$ASCVD==0 & nhanes1$mortstat==0]<-0
nhanes1$ASCVD_mort<-factor(nhanes1$ASCVD_mort, labels = c("Censored",
                                                    "ASCVD mortality", "Other causes"))
# Diabetes definition
nhanes1$dm_hba1c[nhanes1$HbA1c>=6.5]<-1; nhanes1$dm_hba1c[nhanes1$HbA1c<6.5]<-0
nhanes1$dm_G[nhanes1$Glucose>=126]<-1;nhanes1$dm_G[nhanes1$Glucose<126]<-0

# General diabetes variable; DMF#

nhanes1$DM<-NULL
nhanes1$DM[nhanes1$Diabetes==1 & nhanes1$dm_hba1c==1]<-1
nhanes1$DM[nhanes1$Diabetes==0 & nhanes1$dm_hba1c==1]<-1
nhanes1$DM[nhanes1$Diabetes==1 & nhanes1$dm_hba1c==0]<-1
nhanes1$DM[nhanes1$Diabetes==0 & nhanes1$dm_hba1c==0]<-0

nhanes1$DMF<-NULL
nhanes1$DMF[nhanes1$DM==1 & nhanes1$dm_G==1]<-1
nhanes1$DMF[nhanes1$DM==0 & nhanes1$dm_G==1]<-1
nhanes1$DMF[nhanes1$DM==1 & nhanes1$dm_G==0]<-1
nhanes1$DMF[nhanes1$DM==0 & nhanes1$dm_G==0]<-0

nhanes<-nhanes %>% filter(!is.na(HOMA2IR)) %>% filter(!is.na(LDL_S))
nhanes1<-nhanes1 %>% filter(!is.na(HOMA2IR)) %>% filter(!is.na(LDL_S))

base<-rbind(nhanes, nhanes1)
base<- base %>% filter(!is.na(mortstat))

#### Data imputation #####

##Imputation procedure
imp<-mice(base, m=1, maxit=1, seed = 123)
x1<-complete(imp, "long")[,-c(1:3)]
base<-as.data.frame(x1)


base<-base %>% filter(HOMA2IR>0) %>%
  filter(HOMA2IR<20) %>% filter(TG<800) %>% filter(Rem_S>0) %>% filter(DMF==0)

#### Transformation of variables of interest #####

### Build database for models
base1 <- base %>%
  transmute(ASCVD_mort,permth_int,
            Sex, Ethnicity,
            Age = Age,
            BMI = BMI,
            Comorbidities = num_comorb,
            "LDLC-S" = LDL_S,
            "LDLC-M" = LDL_M,
            "Rem-S" = Rem_S,
            "Rem-M" = Rem_M,
            "HOMA2-IR" = HOMA2IR,
            Hypertension = Hypertension,
            HbA1c = HbA1c,
            Smoking = Smoking)

#Percentila 85 de HOMA2IR#
quantile(base$HOMA2IR, p=0.85)
base1$HOMAIR_cat[base$HOMA2IR<2.5]<-0;
base1$HOMAIR_cat[base$HOMA2IR>=2.5]<-1

#Percentila 50 de remanentes de colesterol 
quantile(base$Rem_S)
base1$RemS_cat[base$Rem_S<30]<-0;
base1$RemS_cat[base$Rem_S>=30]<-1

base1$RemS_cat2[base$Rem_S<quantile(base$Rem_S)[2]]<-0;
base1$RemS_cat2[between(base$Rem_S, quantile(base$Rem_S)[2], quantile(base$Rem_S)[3])]<-1
base1$RemS_cat2[between(base$Rem_S, quantile(base$Rem_S)[3], quantile(base$Rem_S)[4])]<-2
base1$RemS_cat2[base$Rem_S>=quantile(base$Rem_S)[4]]<-3

base1$RemS_cat2<-factor(base1$RemS_cat2)

## modelos con remanentes calculados con Martin 
quantile(base$Rem_M)
base1$RemM_cat[base1$`Rem-M`<30]<-0
base1$RemM_cat[base1$`Rem-M`>=30]<-1

base1$rem_IR<-base1$HOMAIR_cat+base1$RemS_cat*2
base1$rem_IR<-factor(base1$rem_IR, labels = c("NormRem+NoIR", "NormRem+IR", "HighRem+NoIR", "HighRem+IR"))

base1$rem_LDL<-(base1$`LDLC-S`>=130)+base1$RemS_cat*2
base1$rem_LDL<-factor(base1$rem_LDL, labels = c("NormRem+NormLDL", "NormRem+HighLDL", "HighRem+NormLDL", "HighRem+HighLDL"))

base1$LDL_cat<-factor(base1$`LDLC-S`>=130, labels = c("Normal LDL-C", "High LDL-C"))

#### Incidence rates per 1,000 person-years #####

## All-cause mortality, Overall ##
mort<- base%>%group_by(mortstat) %>% 
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$mortstat==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## CVD mortality, Overall ## 
mort<- base%>%group_by(Heart_diseases) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Heart_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## Cerebrovascular mortality, Overall ## 
mort<- base%>%group_by(Cerebrovascular_diseases) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Cerebrovascular_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## ASCVD mortality, Overall ## 
mort<- base%>%group_by(ASCVD_mort) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$ASCVD_mort=="ASCVD mortality"]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## All-cause mortality, With IR ##
mort<- base%>%filter(HOMA2IR<2.5) %>%
  group_by(mortstat) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$mortstat==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## CVD mortality, Without IR ## 
mort<- base%>%filter(HOMA2IR<2.5) %>%
  group_by(Heart_diseases) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Heart_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## Cerebrovascular mortality, Without IR ## 
mort<- base%>%filter(HOMA2IR<2.5) %>%
  group_by(Cerebrovascular_diseases) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Cerebrovascular_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## ASCVD mortality, Without IR ## 
mort<- base%>%filter(HOMA2IR<2.5) %>%
  group_by(ASCVD_mort) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$ASCVD_mort=="ASCVD mortality"]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## All-cause mortality, with IR  ##
mort<- base%>%filter(HOMA2IR>=2.5) %>%
  group_by(mortstat) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$mortstat==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## CVD mortality, with IR ## 
mort<- base%>%filter(HOMA2IR>=2.5) %>%
  group_by(Heart_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Heart_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## Cerebrovascular mortality, with IR ## 
mort<- base%>%filter(HOMA2IR>=2.5) %>%
  group_by(Cerebrovascular_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Cerebrovascular_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## ASCVD mortality, with IR ## 
mort<- base%>%filter(HOMA2IR>=2.5) %>%
  group_by(ASCVD_mort) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$ASCVD_mort=="ASCVD mortality"]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## All-cause mortality, RC<30  ##
mort<- base%>%filter(Rem_S<30) %>%
  group_by(mortstat) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$mortstat==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## CVD mortality, RC<30 ## 
mort<- base%>%filter(Rem_S<30) %>%
  group_by(Heart_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Heart_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## Cerebrovascular mortality, RC<30 ## 
mort<- base%>%filter(Rem_S<30) %>%
  group_by(Cerebrovascular_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Cerebrovascular_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## ASCVD mortality, RC<30## 
mort<- base%>%filter(Rem_S<30) %>%
  group_by(ASCVD_mort) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$ASCVD_mort=="ASCVD mortality"]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## All-cause mortality, RC>=30  ##
mort<- base%>%filter(Rem_S>=30) %>%
  group_by(mortstat) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$mortstat==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## CVD mortality, RC>=30 ## 
mort<- base%>%filter(Rem_S>=30) %>%
  group_by(Heart_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Heart_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## Cerebrovascular mortality, RC>=30 ## 
mort<- base%>%filter(Rem_S>=30) %>%
  group_by(Cerebrovascular_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Cerebrovascular_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## ASCVD mortality, RC>=30## 
mort<- base%>%filter(Rem_S>=30) %>%
  group_by(ASCVD_mort) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$ASCVD_mort=="ASCVD mortality"]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## All-cause mortality, No-IR+RC<30  ##
mort<- base%>%filter(Rem_S<30 & HOMA2IR<2.5) %>%
  group_by(mortstat) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$mortstat==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## CVD mortality, No-IR+RC<30 ## 
mort<- base%>%filter(Rem_S<30 & HOMA2IR<2.5) %>%
  group_by(Heart_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Heart_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## Cerebrovascular mortality, No-IR+RC<30 ## 
mort<- base%>%filter(Rem_S<30  & HOMA2IR<2.5) %>%
  group_by(Cerebrovascular_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Cerebrovascular_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## ASCVD mortality, No-IR+RC<30## 
mort<- base%>%filter(Rem_S<30 & HOMA2IR<2.5) %>%
  group_by(ASCVD_mort) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$ASCVD_mort=="ASCVD mortality"]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## All-cause mortality, No-IR+RC>=30  ##
mort<- base%>%filter(Rem_S>=30 & HOMA2IR<2.5) %>%
  group_by(mortstat) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$mortstat==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## CVD mortality, No-IR+RC>=30 ## 
mort<- base%>%filter(Rem_S>=30 & HOMA2IR<2.5) %>%
  group_by(Heart_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Heart_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## Cerebrovascular mortality, No-IR+RC>=30 ## 
mort<- base%>%filter(Rem_S>=30  & HOMA2IR<2.5) %>%
  group_by(Cerebrovascular_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Cerebrovascular_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## ASCVD mortality, No-IR+RC>=30## 
mort<- base%>%filter(Rem_S>=30 & HOMA2IR<2.5) %>%
  group_by(ASCVD_mort) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$ASCVD_mort=="ASCVD mortality"]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## All-cause mortality, IR+RC<30  ##
mort<- base%>%filter(Rem_S<30 & HOMA2IR>=2.5) %>%
  group_by(mortstat) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$mortstat==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## CVD mortality, IR+RC<30 ## 
mort<- base%>%filter(Rem_S<30 & HOMA2IR>=2.5) %>%
  group_by(Heart_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Heart_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## Cerebrovascular mortality, IR+RC<30 ## 
mort<- base%>%filter(Rem_S<30  & HOMA2IR>=2.5) %>%
  group_by(Cerebrovascular_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Cerebrovascular_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## ASCVD mortality, IR+RC<30## 
mort<- base%>%filter(Rem_S<30 & HOMA2IR>=2.5) %>%
  group_by(ASCVD_mort) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$ASCVD_mort=="ASCVD mortality"]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## All-cause mortality, IR+RC>=30  ##
mort<- base%>%filter(Rem_S>=30 & HOMA2IR>=2.5) %>%
  group_by(mortstat) %>%
  summarise(cases=n(), time=sum(permth_int)/12)
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$mortstat==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## CVD mortality, IR+RC>=30 ## 
mort<- base%>%filter(Rem_S>=30 & HOMA2IR>=2.5) %>%
  group_by(Heart_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Heart_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## Cerebrovascular mortality, IR+RC>=30 ## 
mort<- base%>%filter(Rem_S>=30  & HOMA2IR>=2.5) %>%
  group_by(Cerebrovascular_diseases) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$Cerebrovascular_diseases==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## ASCVD mortality, IR+RC>=30## 
mort<- base%>%filter(Rem_S>=30 & HOMA2IR>=2.5) %>%
  group_by(ASCVD_mort) %>%
  summarise(cases=n(), time=sum(permth_int/12))
#Luego extraes el tiempo de cada caso
time<- sum(mort$time)
cases<- mort$cases[mort$ASCVD_mort=="ASCVD mortality"]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se);
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

#### Table 1 ####

## Age ##
tapply(base1$Age, base1$HOMAIR_cat, quantile)
wilcox.test(base1$Age~base1$HOMAIR_cat)

## Sex ##
table(base1$Sex, base1$HOMAIR_cat)
table(base1$Sex, base1$HOMAIR_cat) %>% prop.table(2)
chisq.test(base1$Sex, base1$HOMAIR_cat)

## BMI ##
tapply(base1$BMI, base1$HOMAIR_cat, quantile)
wilcox.test(base1$BMI~base1$HOMAIR_cat)

##  HbA1c ##
tapply(base1$HbA1c, base1$HOMAIR_cat, quantile)
wilcox.test(base1$HbA1c~base1$HOMAIR_cat)

##  Smoking ##
table(base1$Smoking, base1$HOMAIR_cat)
table(base1$Smoking, base1$HOMAIR_cat) %>% prop.table(2)
chisq.test(base1$Sex, base1$HOMAIR_cat)

##  Hypertension ##
table(base1$Hypertension, base1$HOMAIR_cat)
table(base1$Hypertension, base1$HOMAIR_cat) %>% prop.table(2)
chisq.test(base1$Hypertension, base1$HOMAIR_cat)

##  HbA1c ##
tapply(base1$`LDLC-S`, base1$HOMAIR_cat, quantile)
wilcox.test(base1$`LDLC-S`~base1$HOMAIR_cat)

##  Remnant Cholesterol ##
tapply(base1$`Rem-S`, base1$HOMAIR_cat, quantile)
wilcox.test(base1$`Rem-S`~base1$HOMAIR_cat)

##  RC >30mg/dL ##
table(base1$RemS_cat, base1$HOMAIR_cat)
table(base1$RemS_cat, base1$HOMAIR_cat) %>% prop.table(2)
chisq.test(base1$RemS_cat, base1$HOMAIR_cat)

##  All-cause mortality ##
table(base$mortstat, base1$HOMAIR_cat)
table(base$mortstat, base1$HOMAIR_cat) %>% prop.table(2)
chisq.test(base$mortstat, base1$HOMAIR_cat)

##  CV mortality ##
table(base$Heart_diseases, base1$HOMAIR_cat)
table(base$Heart_diseases, base1$HOMAIR_cat) %>% prop.table(2)
chisq.test(base$Heart_diseases, base1$HOMAIR_cat)

##  Stroke mortality ##
table(base$Cerebrovascular_diseases, base1$HOMAIR_cat)
table(base$Cerebrovascular_diseases, base1$HOMAIR_cat) %>% prop.table(2)
chisq.test(base$Cerebrovascular_diseases, base1$HOMAIR_cat)

##  ASCVD mortality ##
table(base$ASCVD, base1$HOMAIR_cat)
table(base$ASCVD, base1$HOMAIR_cat) %>% prop.table(2)
chisq.test(base$ASCVD, base1$HOMAIR_cat)

#### Build Figure 2 ####

fig_2a<-base1 %>% 
  mutate(LDL_cat=factor(LDL_cat, labels = c("<130 mg/dL", "???130 mg/dL"))) %>%
  ggplot(aes(x=`HOMA2-IR`, y=`Rem-S`, col=LDL_cat))+geom_point(alpha=0.3)+
  geom_smooth(method = "lm")+labs(x="HOMA2-IR", y="Remnant Cholesterol (Rem-S, mg/dL)")+
  theme_classic()+scale_x_log10()+scale_y_log10()+scale_color_jama()+
  labs(col="LDL-c \n levels")

fig_2b<-base1 %>%
  mutate(LDL_cat=factor(LDL_cat, labels = c("LDL-c <130 mg/dL", "LDL-c ???130 mg/dL"))) %>%
  mutate(HOMAIR_cat=factor(HOMAIR_cat, labels = c("No-IR", "IR"))) %>%
  ggplot(aes(y=`Rem-S`, x=HOMAIR_cat, fill=HOMAIR_cat))+geom_boxplot()+
  labs(x=" ", y="Remnant Cholesterol (Rem-S, mg/dL)")+
  theme_classic()+scale_y_log10()+scale_fill_jama()+stat_compare_means()+
  labs(fill="Status")+facet_wrap(~LDL_cat)

fig2<-ggarrange(fig_2a, fig_2b, labels = c("A", "B"))

ggsave(file="Figure2.jpg", fig2, bg="transparent",
       width=30, height=12, units=c("cm"), dpi=300, limitsize = FALSE)

#### Linear regression model - No Diabetes ####

# Model for Cholesterol Remnants using Sampson #
m_S<-lm(`Rem-S`~`HOMA2-IR`+Age+Sex+Ethnicity+BMI+Comorbidities+Hypertension+HbA1c+Smoking+`LDLC-S`, base1)
summary(m_S)

f1<-forest_model(m_S)

# Model for Cholesterol Remnants using Martin #
m_M<-lm(`Rem-M`~`HOMA2-IR`+Age+Sex+Ethnicity+BMI+Comorbidities+Hypertension+HbA1c+Smoking+`LDLC-S`, base1)
summary(m_M)

f2<-forest_model(m_M)

ggarrange(f1, f2, labels = c("A", "B")) %>% ggsave(file="SuppFigure3.jpg", bg="transparent",
       width=40, height=12, units=c("cm"), dpi=300, limitsize = FALSE)

## Comparison of RC values according to IR levels ##

t.test(base1$`Rem-S`~base1$HOMAIR_cat)
tapply(base1$`Rem-S`, base1$HOMAIR_cat, sd)

#### Build Figure 2 ####

fig2A<-forest_model(m_S,breaks = c(-0.1, 0.1, 0.5))
fig2b<-forest_model(m_Sdm,breaks = c(-0.1, 0.1, 0.5))

fig2<-ggarrange(fig2A, fig2b, labels = c("A", "B"))

ggsave(file="Figure2.jpg", fig2, bg="transparent",
       width=50, height=15, units=c("cm"), dpi=300, limitsize = FALSE)

figS1a<-forest_model(m_M, breaks = c(-0.1, 0.1, 0.5))
figS1b<-forest_model(m_Mdm,breaks = c(-0.1, 0.1, 0.5))

figS1<-ggarrange(figS1a, figS1B, labels = c("A", "B"))

ggsave(file="FigureS1.jpg", fig2, bg="transparent",
       width=50, height=15, units=c("cm"), dpi=300, limitsize = FALSE)

#### Fine & Gray Model - No Diabetes #####

## Build Fine & Gray model ##
fgdata_cv <- finegray(Surv(permth_int, ASCVD_mort) ~ ., data=base1, 
                      na.action=na.pass)

## Continous Fine & Gray model without HOMA2-IR ##
m1_s <- coxph(Surv(fgstart, fgstop, fgstatus) ~scale(`Rem-S`),weight=fgwt, data=fgdata_cv)
summary(m1_s)
forest_model(m1_s)
m1_s <- coxph(Surv(fgstart, fgstop, fgstatus) ~scale(`Rem-S`)*scale(Age)+Sex+BMI+Comorbidities+Hypertension+scale(HbA1c)+Smoking+Ethnicity+`LDLC-S`,weight=fgwt, data=fgdata_cv)
summary(m1_s)

m1_m <- coxph(Surv(fgstart, fgstop, fgstatus) ~scale(`Rem-M`)*scale(Age)+Sex+BMI+Comorbidities+Hypertension+scale(HbA1c)+Smoking+Ethnicity+`LDLC-M`,weight=fgwt, data=fgdata_cv)
summary(m1_s)

g1<-plot_summs(m1_s, m1_m, model.names = c("Rem-S model", "Rem-M model"), exp=T)+xlab("Fine-Gray model for ASCVD without HOMA2-IR")


## Continous Fine & Gray model without RC ##
m1_s <- coxph(Surv(fgstart, fgstop, fgstatus) ~scale(`HOMA2-IR`),weight=fgwt, data=fgdata_cv)
summary(m1_s)

m2_s <- coxph(Surv(fgstart, fgstop, fgstatus) ~scale(`HOMA2-IR`)*scale(Age)+Sex+BMI+Comorbidities+Hypertension+scale(HbA1c)+Smoking+Ethnicity+`LDLC-S`,weight=fgwt, data=fgdata_cv)
m2_m <- coxph(Surv(fgstart, fgstop, fgstatus) ~scale(`HOMA2-IR`)*scale(Age)+Sex+BMI+Comorbidities+Hypertension+scale(HbA1c)+Smoking+Ethnicity+`LDLC-M`,weight=fgwt, data=fgdata_cv)

g2<-plot_summs(m2_s, m2_m, model.names = c("LDLC-S model", "LDLC-M model"), exp=T)+xlab("Fine-Gray model for ASCVD without Rem-S/Rem-M")

ggarrange(g1, g2, labels = c("A", "B")) %>%
ggsave(file="SuppFigure4.jpg", bg="transparent",
       width=30, height=12, units=c("cm"), dpi=300, limitsize = FALSE)


## Continous Fine & Gray model with HOMA2-IR ##
m1_s <- coxph(Surv(fgstart, fgstop, fgstatus) ~scale(`Rem-S`)*scale(`HOMA2-IR`),weight=fgwt, data=fgdata_cv)
summary(m1_s)

m3_s <- coxph(Surv(fgstart, fgstop, fgstatus) ~scale(`Rem-S`)+Age*scale(`HOMA2-IR`)+Sex+BMI+Comorbidities+Hypertension+scale(HbA1c)+Smoking+Ethnicity+`LDLC-S`,weight=fgwt, data=fgdata_cv)
summary(m1_s)

m3_m <- coxph(Surv(fgstart, fgstop, fgstatus) ~scale(`Rem-M`)+Age*scale(`HOMA2-IR`)+Sex+BMI+Comorbidities+Hypertension+scale(HbA1c)+Smoking+Ethnicity+`LDLC-M`,weight=fgwt, data=fgdata_cv)
summary(m1_s)

g3<-plot_summs(m3_s, m3_m, model.names = c("Rem-S model", "Rem-M model"), exp=T)+xlab("Fine-Gray model for ASCVD with joint Rem-S/Rem-M and HOMA2-IR")


## Categorical Fine & Gray model ##

m4_s <- coxph(Surv(fgstart, fgstop, fgstatus) ~ rem_IR+Age+Sex+BMI+Comorbidities+Hypertension+Ethnicity+HbA1c+Smoking+Ethnicity+`LDLC-S`,weight=fgwt, data=fgdata_cv)
summary(m1_s1)
m4_m <- coxph(Surv(fgstart, fgstop, fgstatus) ~ rem_IR+Age+Sex+BMI+Comorbidities+Hypertension+Ethnicity+HbA1c+Smoking+Ethnicity+`LDLC-M`,weight=fgwt, data=fgdata_cv)
summary(m1_s1)

g4<-plot_summs(m4_s, m4_m, model.names = c("Rem-S model", "Rem-M model"), exp=T)+xlab("Fine-Gray model for ASCVD with joint Rem-S/Rem-M and HOMA2-IR, categorical")

ggarrange(g3, g4, labels = c("A", "B")) %>%
  ggsave(file="SuppFigure5.jpg", bg="transparent",
         width=30, height=12, units=c("cm"), dpi=300, limitsize = FALSE)

#### Build Figure 3 ####
fgdata_cv <- finegray(Surv(permth_int, ASCVD_mort) ~ ., data=base, 
                      na.action=na.pass)

m1_s <- coxph(Surv(fgstart, fgstop, fgstatus) ~ Rem_S*HOMA2IR,weight=fgwt, data=fgdata_cv)
summary(m1_s)

Sim1 <- coxsimLinear(m1_s, b = "HOMA2IR", Xj = seq(0, 10, by = 1),nsim = 100)
fig3a<-simGG(Sim1)+xlab("HOMA2-IR")

Sim2 <- coxsimLinear(m1_s, b = "Rem_S", Xj = seq(0, 100, by = 1),nsim = 100)
fig3b<-simGG(Sim2)+xlab("Remnant cholesterol (Rem-S, mg/dL)")


m1_m <- coxph(Surv(fgstart, fgstop, fgstatus) ~ Rem_M*HOMA2IR,weight=fgwt, data=fgdata_cv)
summary(m1_s)

Sim1 <- coxsimLinear(m1_m, b = "HOMA2IR", Xj = seq(0, 10, by = 1),nsim = 100)
fig3a<-simGG(Sim1)+xlab("HOMA2-IR")

Sim2 <- coxsimLinear(m1_m, b = "Rem_M", Xj = seq(0, 100, by = 1),nsim = 100)
fig3b<-simGG(Sim2)+xlab("Remnant cholesterol (Rem-M, mg/dL)")

ggarrange(fig3a, fig3b, labels = c("A", "B")) %>%
ggsave(file="SuppFigure6.jpg", bg="transparent",
       width=30, height=10, units=c("cm"), dpi=300, limitsize = FALSE)

ggarrange(fig3a, fig3b, labels = c("A", "B")) %>%
ggsave(file="Figure3.jpg", bg="transparent",
       width=30, height=10, units=c("cm"), dpi=300, limitsize = FALSE)


#### Casual mediation model ####

# Model for Cholesterol Remnants using Sampson #
m_S<-lm(`HOMA2-IR`~`Rem-S`, base1)
m1_s <- coxph(Surv(fgstart, fgstop, fgstatus) ~`HOMA2-IR`+`Rem-S`,weight=fgwt, data=fgdata_cv)
summary(m1_s)
BIC(m1_s)
set.seed(1234)
method="Cox"
if (method=="Aalen"){
  lambdas<-aalen1$gamma
  lambdas.var<-aalen1$robvar.gamma
} else if(method=="Cox"){
  lambdas<-m1_s$coef
  lambdas.var<-m1_s$var
}
mediation_ci1(lambdas[2], lambdas[1], lambdas.var[2,2], lambdas.var[2,1], 
              lambdas.var[1,1],
              m_S$coef[2], vcov(m_S)[2,2], G=10^6, method=method)

#### modelos saturado
m_S<-lm(`HOMA2-IR`~`Rem-S`+Age+Sex+Ethnicity+BMI+Comorbidities+Hypertension+HbA1c+Smoking+`LDLC-S`, base1)
fgdata_cv <- finegray(Surv(permth_int/12, ASCVD_mort) ~ ., data=base1, 
                      na.action=na.pass)
m1_s <- coxph(Surv(fgstart, fgstop, fgstatus) ~`HOMA2-IR`*Age+`Rem-S`+Sex+BMI+Comorbidities+Hypertension+HbA1c+Smoking+`LDLC-S`+Ethnicity,weight=fgwt, data=fgdata_cv)
BIC(m1_s)

set.seed(1234)
method="Cox"
if (method=="Aalen"){
  lambdas<-aalen1$gamma
  lambdas.var<-aalen1$robvar.gamma
} else if(method=="Cox"){
  lambdas<-m1_s$coef
  lambdas.var<-m1_s$var
}
mediation_ci1(lambdas[3], lambdas[1], lambdas.var[3,3], lambdas.var[3,1], 
              lambdas.var[1,1],
              m_S$coef[2], vcov(m_S)[2,2], G=10^6, method=method)

IC_Sup<-(log(1.005668)/log(1.005246))*100
IC_Inf<-(log(1.011201)/log(1.014243))*100
Prop_med<-(log(1.008415)/log(1.009731))*100
Prop_med;IC_Inf;IC_Sup

# Model for Cholesterol Remnants using Martin #
m_S<-lm(`HOMA2-IR`~`Rem-M`+Age+Sex+Ethnicity+BMI+Comorbidities+Hypertension+HbA1c+Smoking+`LDLC-M`, base1)
summary(m_S)
fgdata_cv <- finegray(Surv(permth_int/12, ASCVD_mort) ~ ., data=base1, 
                      na.action=na.pass)
m1_s <- coxph(Surv(fgstart, fgstop, fgstatus) ~`HOMA2-IR`*Age+`Rem-M`+Sex+BMI+Comorbidities+Hypertension+HbA1c+Smoking+`LDLC-M`+Ethnicity,weight=fgwt, data=fgdata_cv)
summary(m1_s)

set.seed(1234)
method="Cox"
if (method=="Aalen"){
  lambdas<-aalen1$gamma
  lambdas.var<-aalen1$robvar.gamma
} else if(method=="Cox"){
  lambdas<-m1_s$coef
  lambdas.var<-m1_s$var
}
mediation_ci1(lambdas[3], lambdas[1], lambdas.var[3,3], lambdas.var[3,1], 
              lambdas.var[1,1],
              m_S$coef[2], vcov(m_S)[2,2], G=10^6, method=method)


IC_infm<-(log(1.016420)/log(1.019307))*100
IC_Supm<-(log(1.008307)/log(1.006086))*100
Prop_medm<-(log(1.012332)/log(1.01267))*100
Prop_medm;IC_infm;IC_Supm

#### Modelos RC mediador de HOMA

m_S<-lm(`Rem-S`~`HOMA2-IR`+Age+Sex+Ethnicity+BMI+Comorbidities+Hypertension+HbA1c+Smoking+`LDLC-M`, base1)
summary(m_S)
fgdata_cv <- finegray(Surv(permth_int/12, ASCVD_mort) ~ ., data=base1, 
                      na.action=na.pass)
m1_s <- coxph(Surv(fgstart, fgstop, fgstatus) ~HOMAIR_cat*Age+`Rem-S`+Sex+BMI+Comorbidities+Hypertension+HbA1c+Smoking+`LDLC-M`+Ethnicity,weight=fgwt, data=fgdata_cv)
summary(m1_s)

set.seed(1234)
method="Cox"
if (method=="Aalen"){
  lambdas<-aalen1$gamma
  lambdas.var<-aalen1$robvar.gamma
} else if(method=="Cox"){
  lambdas<-m1_s$coef
  lambdas.var<-m1_s$var
}
mediation_ci1(lambdas[1], lambdas[3], lambdas.var[3,1], lambdas.var[3,3], 
              lambdas.var[1,1],
              m_S$coef[2], vcov(m_S)[2,2], G=10^6, method=method)


#### Supplementary Figure 1 ####

rphA<-base1  %>%select(
  Age, BMI, `HOMA2-IR`, `Rem-S`, `Rem-M`, `LDLC-S`, `LDLC-M`, HbA1c
) %>% as.matrix %>% rcorr(type = "spearman");  F_names <- c(
  "Age", "BMI", "HOMA2-IR", "Rem-S", "Rem-M", "LDL-C (S)","LDL-C (M)", "HbA1c"
); rownames(rphA$r) <- F_names; colnames(rphA$r) <- F_names

corrplot_full<-function(x,y){print(x); print(y)}
FigSX_A <- as.ggplot(
  ~corrplot_full(corrplot(rphA$r,method="number",type="lower",add=F,p.mat=rphA$p,sig.level=.05/50,tl.pos="tl", cl.pos="r"),
                 corrplot(rphA$r, method="circle", type="upper", add=T, p.mat=rphA$p, sig.level=.05/50, tl.pos = "n"))) %>% 
  annotate_figure(top = text_grob("Correlations in individuals without diabetes", face = "bold", size = 15, vjust = 1.5))

FigSX_A %>%
  ggsave(file = "SuppFig1.jpg", bg = "transparent", width = 25,  height = 20, units=c("cm"), dpi = 300, limitsize = FALSE)

#### Build Supp Figure 2 ####

fig_2a<-base1 %>% 
  mutate(LDL_cat=factor(`LDLC-M`>=130, labels = c("<130 mg/dL", "???130 mg/dL"))) %>%
  ggplot(aes(x=`HOMA2-IR`, y=`Rem-S`, col=LDL_cat))+geom_point(alpha=0.3)+
  geom_smooth(method = "lm")+labs(x="HOMA2-IR", y="Remnant Cholesterol (Rem-M, mg/dL)")+
  theme_classic()+scale_x_log10()+scale_y_log10()+scale_color_jama()+
  labs(col="LDL-c \n levels")

fig_2b<-base1 %>%
  mutate(LDL_cat=factor(`LDLC-M`>=130, labels = c("LDL-c (M) <130 mg/dL", "LDL-c (M) ???130 mg/dL"))) %>%
  mutate(HOMAIR_cat=factor(HOMAIR_cat, labels = c("No-IR", "IR"))) %>%
  ggplot(aes(y=`Rem-S`, x=HOMAIR_cat, fill=HOMAIR_cat))+geom_boxplot()+
  labs(x=" ", y="Remnant Cholesterol (Rem-M, mg/dL)")+
  theme_classic()+scale_y_log10()+scale_fill_jama()+stat_compare_means()+
  labs(fill="Status")+facet_wrap(~LDL_cat)

fig2<-ggarrange(fig_2a, fig_2b, labels = c("A", "B"))

ggsave(file="SuppFigure2.jpg", fig2, bg="transparent",
       width=30, height=12, units=c("cm"), dpi=300, limitsize = FALSE)
