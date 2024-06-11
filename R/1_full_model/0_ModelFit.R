################################################################################
#required packages
list.of.packages <- c(
    "doParallel", "parallel","foreach","pdftools","plotly",
    "pbapply","dplyr", "tidyr", "parallel",
    "scales","effects","psych", "glmmTMB", "lme4", "lmerTest","here","rlist","ggtext","gridExtra","grid","lattice","viridis") 


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)
################################################################################
#define the data repository
dir_in="/home/rbertrand/W/Bioshift/GD_study/data_Brunno_v02022024/" #to change accordingly to the location of the data
dir_out="/home/rbertrand/W/Bioshift/GD_study/boot_analysis/REMLF2" #to change. It's the repository where the results are saved

# Load data
setwd(dir_in)
mydataset <- read.csv2("gen_data_final_fonseca.csv",sep=",",dec=".",h=T) #file path in GitHub: /adaptive-potential/Data

#Data selection
## Latitude data
mydatatogo <- mydataset  %>%
    dplyr::filter(Type == "LAT",
                  shift_vel_sign == "pospos" | shift_vel_sign == "negneg", # Select only shifts in the same direction of velocity
                  SHIFT != 0, # remove non-significant shifts
                  # Nucleotide_diversity > 0 # select only GD values > 0
    ) %>% 
    mutate(
        # Climate velocity
        vel = as.numeric(velocity),
        vel_abs = abs(vel),
        vel_abs_log = log(vel_abs),
        vel_abs_log1p = log1p(vel_abs),
        
        # Shift
        SHIFT = SHIFT, # to deal with zero shift
        SHIFT_abs = abs(SHIFT),
        SHIFT_abs_log = log(SHIFT_abs),
        SHIFT_abs_log1p = log1p(SHIFT_abs),
        
        # Genetic diversity
        GD = Nucleotide_diversity, 
        GD_sqrt = sqrt(GD),
        GD_log = log(GD),
        GD_log1p = log1p(GD),
        
        # Methods
        Lat = abs(Latitude),
        Lat_band = round(Lat,0),
        ID.area = scale(ID.area),
        DUR = scale(DUR),
        LogExtent = log(Extent),
        START = scale(START),
        Param = factor(Param),
        Group = factor(Group),
        spp = factor(spp), 
        ExtentF = cut(Extent,
                      ordered_result = TRUE,
                      breaks=seq(min(Extent), max(Extent), length.out=10),
                      include.lowest=TRUE),
        NtempUnitsF = cut(NtempUnits,
                          ordered_result = TRUE,
                          breaks=seq(min(NtempUnits), max(NtempUnits), length.out=10),
                          include.lowest=TRUE)) %>%
    
    dplyr::select(
        # Genetic diversity
        GD, GD_log, GD_log1p, GD_sqrt, TajimasD,
        # Shift
        SHIFT, SHIFT_abs, SHIFT_abs_log, SHIFT_abs_log1p, 
        # SHIFT_cor, SHIFT_cor_abs, SHIFT_cor_abs_log, SHIFT_cor_raw, SHIFT_abs_log_scale,
        # Velocity
        vel, vel_abs, vel_abs_log, vel_abs_log1p, 
        trend.mean,
        # Methods + Taxonomy
        Article_ID, 
        Hemisphere,
        shift_vel_sign,
        Lat, # latitudinal position where GD was collected
        Lat_band,
        DUR, Nperiodes, LogNtempUnits, NtempUnits, Extent, LogExtent, ContinuousGrain, Quality, PrAb, ExtentF, NtempUnitsF,
        Param, Group, spp, Class, Order, Family, Genus, 
        ECO, Uncertainty_Parameter, Uncertainty_Distribution, Grain_size, Data, Article_ID,
        # whereas GD data overlaps with the range shift's study area
        overlap_SA,
        # Distance from GD data to Bioshifts SA
        dist_cent, # distance to the centroid of the SA
        dist_N, # distance to the North edge of the SA (equivalent to the distance to the LE at the North hemisphere)
        dist_S, # distance to the South edge of the SA (equivalent to the distance to the TE at the North hemisphere)
        dist_edge # distance to the closest edge of the SA
    ) 


# transform continuous variables
cont_vars <- c(1:14, 20:26, 44:47)

mydatatogo[,cont_vars] <- lapply(mydatatogo[,cont_vars], as.numeric)
mydatatogo[,-cont_vars] <- lapply(mydatatogo[,-cont_vars], function(x) factor(x, levels = unique(x)))


## Filter Classes with at least 5 species per Param

n_sps = 10

test <- mydatatogo %>%
    group_by(Class,Param) %>%
    dplyr::summarise(N_spp = length(unique(spp))) %>% # how many species per parameter?
    dplyr::filter(N_spp >= n_sps) # select classes with > n_sps per param

test <- mutate(test, Class_Param = paste(Class, Param))

mydatatogo <- mydatatogo %>%
    mutate(Class_Param = paste(Class, Param)) %>%
    filter(Class_Param %in% test$Class_Param) %>%
    select(-Class_Param)

mydatatogo[,-cont_vars] <- lapply(mydatatogo[,-cont_vars], function(x) factor(x, levels = unique(x)))

## Extra fixes

# Fix param levels
mydatatogo$Param <- relevel(mydatatogo$Param, ref = "O") 

####
# Calculate N obs per species for model weights
N_obs_spp <- mydatatogo %>%
    group_by(spp) %>%
    tally()

N_obs_spp$weight_obs_spp <- (1/N_obs_spp$n)*(1/nrow(N_obs_spp))

mydatatogo <- merge(mydatatogo,N_obs_spp,by="spp")
sum(mydatatogo$weight_obs_spp)

###############################################
#below four models are fitted
##allW: a model that considers the whole dataset (the position is accounting for in the model as a fixed and additive effect only due to the imbalance in data among positions) weighting for the imbalance in he number of obs among species (ie, in the model every species will have the same weight whatever the number of obs)
##CEW, LEW and TEW: three independent models that model centroid, leading and trailing edges range shifts, respectively. Observations used by each model are weighted to control for imbalance in he number of obs among species (ie, in the model every species will have the same weight whatever the number of obs)
##allW2: similar model as allW but it differs in two ways: (i) the obs are weighted in order to control for the imbalance in both the number of obs among species (ie, in the model every species will have the same weight whatever the number of obs) and the number of obs among the positions.
###(ii) as imbalance in data is controled the position is accounted for an additive and interaction effect (with climate velocity and genetic diversity in the model).

##allW fit
gam_formula1 <- as.formula(
    "SHIFT_abs ~ scale(vel_abs) * scale(GD) + 
    Param + LogNtempUnits + LogExtent + ContinuousGrain + PrAb + Quality + 
    (1 | Class)")
x1=data.frame(table(mydatatogo$spp))
x1=subset(x1,Freq>0)
x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
mydatatogo=merge(mydatatogo[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")

gam1 <- glmmTMB(gam_formula1, 
                family = Gamma(link = "log"),
                weights = weight_obs_spp,
                data = mydatatogo)
summary(gam1)
round(MuMIn::r.squaredGLMM(gam1),2)
AIC(gam1)

#bootstrap model computation
#homemade cluster
no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)
s1=10 #the seed to make reproducible random staff
nB=11000 #the number of bootstrap
library(foreach)
library(parallel)
library(doParallel)
ex=foreach(i=1:nB,.inorder=F,.packages=c('glmmTMB','lme4','MuMIn'))%dopar%{
    gc(reset=T)
    print(i)
    set.seed(s1+i)
    n1=sample(1:nrow(mydatatogo),rep=T,size=nrow(mydatatogo))
    mydatatogo2=mydatatogo[n1,]
    x1=data.frame(table(mydatatogo2$spp))
    x1=subset(x1,Freq>0)
    x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
    mydatatogo2=merge(mydatatogo2[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")
    
    gamX1 <- glmmTMB(gam_formula1, 
                     family = Gamma(link = "log"),
                     weights = weight_obs_spp,
                     data = mydatatogo2)
    r2=round(MuMIn::r.squaredGLMM(gamX1),2)
    r2=data.frame(r2m=r2[1,1],r2c=r2[1,2],moy_GD=mean(mydatatogo2$GD), sd_GD=sd(mydatatogo2$GD),moy_VC=mean(mydatatogo2$vel_abs), sd_VC=sd(mydatatogo2$vel_abs),nB=i)
    coeff=data.frame(summary(gamX1)$coeff$cond,var=row.names(summary(gamX1)$coeff$cond),nB=i)
    randEff=as.data.frame(ranef(gamX1,condVar=T))
    randEff$nB=i
    return(list(r2,coeff,randEff))
}
parallel::stopCluster(cl = my.cluster2)

nom="allW"
randEff=rlist::list.rbind(lapply(ex,"[[",3))
coeff=rlist::list.rbind(lapply(ex,"[[",2))
r2=rlist::list.rbind(lapply(ex,"[[",1))
setwd(dir_out)
write.table(r2,paste0("R2_",nom,".csv"),sep=";",dec=".",row=F)
write.table(coeff,paste0("coeff_",nom,".csv"),sep=";",dec=".",row=F)
write.table(randEff,paste0("randEff_",nom,".csv"),sep=";",dec=".",row=F)

##LEW model fit
gam_formulaLE <- as.formula(
    "SHIFT_abs ~ scale(vel_abs) * scale(GD) + 
    LogNtempUnits + LogExtent + ContinuousGrain + PrAb + Quality + 
    (1 | Class)")
mydatatogoS=subset(mydatatogo,Param=="LE")
x1=data.frame(table(mydatatogoS$spp))
x1=subset(x1,Freq>0)
x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
mydatatogoS=merge(mydatatogoS[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")

gamLE <- glmmTMB(gam_formulaLE, 
                 family = Gamma(link = "log"),
                 weights = weight_obs_spp,
                 data = mydatatogoS)

summary(gamLE)
round(MuMIn::r.squaredGLMM(gamLE),2)
AIC(gamLE)


s1=10 #the seed to make reproducible random staff
nB=11000 #the number of bootstrap

#bootstrap model computation for LE model
#homemade cluster
no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

s1=10 #the seed to make reproducible random staff
nB=11000 #the number of bootstrap
ex=foreach(i=1:nB,.inorder=F,.packages=c('glmmTMB','lme4','MuMIn'))%dopar%{
    gc(reset=T)
    print(i)
    set.seed(s1+i)
    n1=sample(1:nrow(mydatatogoS),rep=T,size=nrow(mydatatogoS))
    mydatatogo2=mydatatogoS[n1,]
    x1=data.frame(table(mydatatogo2$spp))
    x1=subset(x1,Freq>0)
    x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
    mydatatogo2=merge(mydatatogo2[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")
    
    gamX1 <- glmmTMB(gam_formulaLE, 
                     family = Gamma(link = "log"),
                     weights = weight_obs_spp,
                     data = mydatatogo2)
    r2=round(MuMIn::r.squaredGLMM(gamX1),2)
    r2=data.frame(r2m=r2[1,1],r2c=r2[1,2],moy_GD=mean(mydatatogo2$GD), sd_GD=sd(mydatatogo2$GD),moy_VC=mean(mydatatogo2$vel_abs), sd_VC=sd(mydatatogo2$vel_abs),nB=i)
    coeff=data.frame(summary(gamX1)$coeff$cond,var=row.names(summary(gamX1)$coeff$cond),nB=i)
    randEff=as.data.frame(ranef(gamX1,condVar=T))
    randEff$nB=i
    return(list(r2,coeff,randEff))
}
parallel::stopCluster(cl = my.cluster2)

nom="LEW"
randEff=rlist::list.rbind(lapply(ex,"[[",3))
coeff=rlist::list.rbind(lapply(ex,"[[",2))
r2=rlist::list.rbind(lapply(ex,"[[",1))
setwd(dir_out)
write.table(r2,paste0("R2_",nom,".csv"),sep=";",dec=".",row=F)
write.table(coeff,paste0("coeff_",nom,".csv"),sep=";",dec=".",row=F)
write.table(randEff,paste0("randEff_",nom,".csv"),sep=";",dec=".",row=F)

##CEW model fit
gam_formulaCE <- as.formula(
    "SHIFT_abs ~ scale(vel_abs) * scale(GD) + 
    LogNtempUnits + LogExtent + ContinuousGrain + PrAb + Quality + 
    (1 | Class)")
mydatatogoS=subset(mydatatogo,Param=="O")
x1=data.frame(table(mydatatogoS$spp))
x1=subset(x1,Freq>0)
x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
mydatatogoS=merge(mydatatogoS[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")

gamCE <- glmmTMB(gam_formulaCE, 
                 family = Gamma(link = "log"),
                 weights = weight_obs_spp,
                 data = mydatatogoS)

summary(gamCE)
round(MuMIn::r.squaredGLMM(gamCE),2)
AIC(gamCE)

#CE model fit
#homemade cluster
no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

s1=10 #the seed to make reproducible random staff
nB=11000 #the number of bootstrap
ex=foreach(i=1:nB,.inorder=F,.packages=c('glmmTMB','lme4','MuMIn'))%dopar%{
    gc(reset=T)
    print(i)
    set.seed(s1+i)
    n1=sample(1:nrow(mydatatogoS),rep=T,size=nrow(mydatatogoS))
    mydatatogo2=mydatatogoS[n1,]
    x1=data.frame(table(mydatatogo2$spp))
    x1=subset(x1,Freq>0)
    x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
    mydatatogo2=merge(mydatatogo2[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")
    
    gamX1 <- glmmTMB(gam_formulaCE, 
                     family = Gamma(link = "log"),
                     weights = weight_obs_spp,
                     data = mydatatogo2)
    r2=round(MuMIn::r.squaredGLMM(gamX1),2)
    r2=data.frame(r2m=r2[1,1],r2c=r2[1,2],moy_GD=mean(mydatatogo2$GD), sd_GD=sd(mydatatogo2$GD),moy_VC=mean(mydatatogo2$vel_abs), sd_VC=sd(mydatatogo2$vel_abs),nB=i)
    coeff=data.frame(summary(gamX1)$coeff$cond,var=row.names(summary(gamX1)$coeff$cond),nB=i)
    randEff=as.data.frame(ranef(gamX1,condVar=T))
    randEff$nB=i
    return(list(r2,coeff,randEff))
}
parallel::stopCluster(cl = my.cluster2)

nom="CEW"
randEff=rlist::list.rbind(lapply(ex,"[[",3))
coeff=rlist::list.rbind(lapply(ex,"[[",2))
r2=rlist::list.rbind(lapply(ex,"[[",1))
setwd(dir_out)
write.table(r2,paste0("R2_",nom,".csv"),sep=";",dec=".",row=F)
write.table(coeff,paste0("coeff_",nom,".csv"),sep=";",dec=".",row=F)
write.table(randEff,paste0("randEff_",nom,".csv"),sep=";",dec=".",row=F)

#TEW model fit
gam_formulaTE <- as.formula(
    "SHIFT_abs ~ scale(vel_abs) * scale(GD) + 
    LogNtempUnits + LogExtent + ContinuousGrain + PrAb + Quality + 
    (1 | Class)")
mydatatogoS=subset(mydatatogo,Param=="TE")
x1=data.frame(table(mydatatogoS$spp))
x1=subset(x1,Freq>0)
x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
mydatatogoS=merge(mydatatogoS[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")

gamTE <- glmmTMB(gam_formulaTE, 
                 family = Gamma(link = "log"),
                 weights = weight_obs_spp,
                 data = mydatatogoS)

summary(gamTE)
round(MuMIn::r.squaredGLMM(gamTE),2)
AIC(gamTE)

#bootstrap model computation for LTEW
#homemade cluster
no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

s1=10 #the seed to make reproducible random staff
nB=11000 #the number of bootstrap
ex=foreach(i=1:nB,.inorder=F,.packages=c('glmmTMB','lme4','MuMIn'))%dopar%{
    gc(reset=T)
    print(i)
    set.seed(s1+i)
    n1=sample(1:nrow(mydatatogoS),rep=T,size=nrow(mydatatogoS))
    mydatatogo2=mydatatogoS[n1,]
    x1=data.frame(table(mydatatogo2$spp))
    x1=subset(x1,Freq>0)
    x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
    mydatatogo2=merge(mydatatogo2[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")
    
    gamX1 <- glmmTMB(gam_formulaTE, 
                     family = Gamma(link = "log"),
                     weights = weight_obs_spp,
                     data = mydatatogo2)
    r2=round(MuMIn::r.squaredGLMM(gamX1),2)
    r2=data.frame(r2m=r2[1,1],r2c=r2[1,2],moy_GD=mean(mydatatogo2$GD), sd_GD=sd(mydatatogo2$GD),moy_VC=mean(mydatatogo2$vel_abs), sd_VC=sd(mydatatogo2$vel_abs),nB=i)
    coeff=data.frame(summary(gamX1)$coeff$cond,var=row.names(summary(gamX1)$coeff$cond),nB=i)
    randEff=as.data.frame(ranef(gamX1,condVar=T))
    randEff$nB=i
    return(list(r2,coeff,randEff))
}
parallel::stopCluster(cl = my.cluster2)

nom="TEW"
randEff=rlist::list.rbind(lapply(ex,"[[",3))
coeff=rlist::list.rbind(lapply(ex,"[[",2))
r2=rlist::list.rbind(lapply(ex,"[[",1))
setwd(dir_out)
write.table(r2,paste0("R2_",nom,".csv"),sep=";",dec=".",row=F)
write.table(coeff,paste0("coeff_",nom,".csv"),sep=";",dec=".",row=F)
write.table(randEff,paste0("randEff_",nom,".csv"),sep=";",dec=".",row=F)

##allW2 model fit
gam_formula2 <- as.formula(
    "SHIFT_abs ~ scale(vel_abs) * scale(GD) * Param +
    LogNtempUnits + LogExtent + ContinuousGrain + PrAb + Quality + 
    (1 | Class)")
mydatatogoCE=subset(mydatatogo,Param=="O")
x1=data.frame(table(mydatatogoCE$spp))
x1=subset(x1,Freq>0)
x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
mydatatogoCE=merge(mydatatogoCE[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")

mydatatogoLE=subset(mydatatogo,Param=="LE")
x1=data.frame(table(mydatatogoLE$spp))
x1=subset(x1,Freq>0)
x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
mydatatogoLE=merge(mydatatogoLE[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")

mydatatogoTE=subset(mydatatogo,Param=="TE")
x1=data.frame(table(mydatatogoTE$spp))
x1=subset(x1,Freq>0)
x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
mydatatogoTE=merge(mydatatogoTE[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")

mydatatogoTE$weight_obs_spp=mydatatogoTE$weight_obs_spp*(1/3)
mydatatogoCE$weight_obs_spp=mydatatogoCE$weight_obs_spp*(1/3)
mydatatogoLE$weight_obs_spp=mydatatogoLE$weight_obs_spp*(1/3)

mydatatogo1=rbind(mydatatogoTE,mydatatogoCE,mydatatogoLE)
sum(mydatatogo1$weight_obs_spp)
# Model with zero GD + fixed-effect methods
gam2 <- glmmTMB(gam_formula2, 
                family = Gamma(link = "log"),
                weights = weight_obs_spp,
                data = mydatatogo)
summary(gam2)
round(MuMIn::r.squaredGLMM(gam2),2)
AIC(gam2)

#bootstrap model computation for allW2
#homemade cluster
no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

s1=10 #the seed to make reproducible random staff
nB=11000 #the number of bootstrap
ex=foreach(i=1:nB,.inorder=F,.packages=c('glmmTMB','lme4','MuMIn'))%dopar%{
    gc(reset=T)
    print(i)
    set.seed(s1+i)
    n1=sample(1:nrow(mydatatogo),rep=T,size=nrow(mydatatogo))
    mydatatogo2=mydatatogo[n1,]
    mydatatogoCE=subset(mydatatogo2,Param=="O")
    x1=data.frame(table(mydatatogoCE$spp))
    x1=subset(x1,Freq>0)
    x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
    mydatatogoCE=merge(mydatatogoCE[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")
    
    mydatatogoLE=subset(mydatatogo2,Param=="LE")
    x1=data.frame(table(mydatatogoLE$spp))
    x1=subset(x1,Freq>0)
    x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
    mydatatogoLE=merge(mydatatogoLE[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")
    
    mydatatogoTE=subset(mydatatogo2,Param=="TE")
    x1=data.frame(table(mydatatogoTE$spp))
    x1=subset(x1,Freq>0)
    x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
    mydatatogoTE=merge(mydatatogoTE[,1:48],x1[,c(1,3)],by.x="spp",by.y="Var1")
    
    mydatatogoTE$weight_obs_spp=mydatatogoTE$weight_obs_spp*(1/3)
    mydatatogoCE$weight_obs_spp=mydatatogoCE$weight_obs_spp*(1/3)
    mydatatogoLE$weight_obs_spp=mydatatogoLE$weight_obs_spp*(1/3)
    
    mydatatogo2=rbind(mydatatogoTE,mydatatogoCE,mydatatogoLE)
    #sum(mydatatogo2$weight_obs_spp)
    
    gamX1 <- glmmTMB(gam_formula2, 
                     family = Gamma(link = "log"),
                     weights = weight_obs_spp,
                     data = mydatatogo2)
    r2=round(MuMIn::r.squaredGLMM(gamX1),2)
    r2=data.frame(r2m=r2[1,1],r2c=r2[1,2],moy_GD=mean(mydatatogo2$GD), sd_GD=sd(mydatatogo2$GD),moy_VC=mean(mydatatogo2$vel_abs), sd_VC=sd(mydatatogo2$vel_abs),nB=i)
    coeff=data.frame(summary(gamX1)$coeff$cond,var=row.names(summary(gamX1)$coeff$cond),nB=i)
    randEff=as.data.frame(ranef(gamX1,condVar=T))
    randEff$nB=i
    return(list(r2,coeff,randEff))
}
parallel::stopCluster(cl = my.cluster2)

nom="allW2"
randEff=rlist::list.rbind(lapply(ex,"[[",3))
coeff=rlist::list.rbind(lapply(ex,"[[",2))
r2=rlist::list.rbind(lapply(ex,"[[",1))
setwd(dir_out)
write.table(r2,paste0("R2_",nom,".csv"),sep=";",dec=".",row=F)
write.table(coeff,paste0("coeff_",nom,".csv"),sep=";",dec=".",row=F)
write.table(randEff,paste0("randEff_",nom,".csv"),sep=";",dec=".",row=F)


#################################
#### Bootstrap model analysis####
#################################
#bootstrap test: H0= coefficient equal to 0; H1= coefficient higher than 0
test1=function(x,mu=0){
    return(1-(length(x[x>mu])/length(x)))
}

#bootstrap test: H0= coefficient equal to 0; H1= coefficient lower than 0
test2=function(x,mu=0){
    return(1-(length(x[x<mu])/length(x)))
}

setwd(dir_out)
mod=c("allW","CEW","TEW","LEW","allW2")

a=1
s1=10 #the seed to make reproducible random staff
nB=5000 #the number of bootstrap
for(i in 1:length(mod)){
    print(a)
    r2=read.csv2(paste0("R2_",mod[i],".csv"),sep=";",dec=".")
    #r2$nB=1:nrow(r2)
    r2a=subset(r2,is.na(r2m)==F)
    print(paste(mod[i],': ',nrow(r2a),sep=""))
    if(nrow(r2a)>=nS){ #random selection of nS bootstrap model (criteria: model convergence)
        set.seed(s1)
        r2a=r2a[sample(1:nrow(r2a),size=nS,replace=F),]
    }
    if(nrow(r2a)<nS){
        print(paste0("sampling size is less than ",nS," bootstraps"))
    }else{
        resR2=data.frame(type="r2m",moy=mean(r2a$r2m),sd=sd(r2a$r2m),median=median(r2a$r2m),q025=quantile(r2a$r2m,probs=0.025),q975=quantile(r2a$r2m,probs=0.975))
        resR2=rbind(resR2,data.frame(type="r2c",moy=mean(r2a$r2c),sd=sd(r2a$r2c),median=median(r2a$r2c),q025=quantile(r2a$r2c,probs=0.025),q975=quantile(r2a$r2c,probs=0.975)))
        #resR2$n_sing=nrow(subset(r2a,sing==T))
        resR2$model=mod[i]
        
        c1=read.csv2(paste("coeff_",mod[i],".csv",sep=""),sep=";",dec=".",h=T)
        #c1$nB=rep(1:nrow(r2),each=nrow(c1)/6000)
        c1=merge(c1,data.frame(nB=r2a$nB),by.x="nB",by.y="nB")
        rr3=c1
        
        res=data.frame(var=names(tapply(rr3$Estimate,rr3$var,mean)),moy=tapply(rr3$Estimate,rr3$var,mean),sd=tapply(rr3$Estimate,rr3$var,sd),median=tapply(rr3$Estimate,rr3$var,median),q025=tapply(rr3$Estimate,rr3$var,quantile,probs=0.025),
                       p975=tapply(rr3$Estimate,rr3$var,quantile,probs=0.975), pv.inf0=tapply(rr3$Estimate,rr3$var,test2,mu=0),pv.sup0=tapply(rr3$Estimate,rr3$var,test1,mu=0))
        res$model=mod[i]
        if(a==1){
            resR2_ok=resR2
            res_ok=res
        }
        else{
            resR2_ok=rbind(resR2_ok,resR2)
            res_ok=rbind(res_ok,res)
        }
        #write.table(rr3,paste("coeff_",mod[i],"_",nS,".csv",sep=""),sep=";",dec=".",row=F)
        #Coefficient value for each of the 5000 selected bootstrap models
        #nB= ID of the bootstrap model
        #Estimate= coefficient value computed from every boostrap model
        #Std..Error= Standard error computed from every boostrap model
        #t.value= t-value computed from every boostrap model
        #var= names of the term
        
        #write.table(r2a,paste("r2_",mod[i],"_",nS,".csv",sep=""),sep=";",dec=".",row=F)
        #R2 values for each of the 5000 selected bootstrap models :
        #r2m= marginal R2 of every bootstrap model (each line is a different model)
        #r2c= conditional R2 of every bootstrap model
        #warnM= inform for warning during the model fit (such as convergence issue)
        #sing= output of the model singularity test (FALSE= no signularity; TRUE= singularity issue) 
        #nB= ID of the bootstrap model
        
        a=a+1
    }
}

write.table(resR2_ok,"summary_R2.csv",sep=";",dec=".",row=F) #used to make table S5
#R2 statistic computed from 5000 boostrap model:
#type= R2 type (r2m= marginal R2; r2c= conditional R2)
#moy= mean R2 value
#sd= R2 standard deviation value
#median= median R2 value
#q025= 2.5% quantile of the bootstrap distribution
#q975= 97.5% quantile of the bootstrap distribution
#n_sing= number of singular model
#model= model name abbreviation

write.table(res_ok,"summary_coeff.csv",sep=";",dec=".",row=F) #used to make table S6
#coefficient statistics
#Var =  Term of the model
#moy = mean effect computed from the 5000 bootstrap model
#sd = standard deviation computed from the 5000 bootstrap model
#med = median effect computed from the 5000 bootstrap model
#q025 = 2.5% quantile of the bootstrap coefficient distribution
#q975 = 97.5% quantile of the bootstrap coefficient distribution
#pv.sup0 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 0; H1= coefficient greater than 0)
#pv.inf0 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 0; H1= coefficient lower than 0)
#model= model name abbreviation

#significant positive effect
subset(res_ok,model=="allW" & pv.sup0<0.05)
subset(res_ok,model=="CEW" & pv.sup0<0.05)
subset(res_ok,model=="LEW" & pv.sup0<0.05)
subset(res_ok,model=="TEW" & pv.sup0<0.05)
subset(res_ok,model=="allW2" & pv.sup0<0.05)

#significant negative effect
subset(res_ok,model=="allW" & pv.inf0<0.05)
subset(res_ok,model=="CEW" & pv.inf0<0.05)
subset(res_ok,model=="LEW" & pv.inf0<0.05)
subset(res_ok,model=="TEW" & pv.inf0<0.05)
subset(res_ok,model=="allW2" & pv.inf0<0.05)

#non-significant negative effect
subset(res_ok,model=="allW" & pv.inf0>=0.05 & pv.sup0>=0.05)
subset(res_ok,model=="CEW" & pv.inf0>=0.05 & pv.sup0>=0.05)
subset(res_ok,model=="LEW" & pv.inf0>=0.05 & pv.sup0>=0.05)
subset(res_ok,model=="TEW" & pv.inf0>=0.05 & pv.sup0>=0.05)
subset(res_ok,model=="allW2" & pv.inf0>=0.05 & pv.sup0>=0.05)

#################################################################################
####testing climate velocity value for which genetic diversity is significant####
#################################################################################
setwd(dir_out)
nS=5000#number of boostrap to use to infer stat

#allW model
c1=read.csv2(paste("coeff_allW.csv",sep=""),sep=";",dec=".",h=T)
dsel=mydatatogo
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)

no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

ex=foreach(i=1:nS,.combine=rbind,.inorder=F)%dopar%{
    gc(reset=T)
    print(i)
    c1a=subset(c1,nB==i)
    eff=c1a$Estimate[c1a$var=="scale(GD)"]+(c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]*((v1-mean(dsel$vel_abs))/sd(dsel$vel_abs)))
    tmp=data.frame(vel_abs=v1,GDeff=eff,nB=i)
    return(tmp)
}
parallel::stopCluster(cl = my.cluster2)

res=data.frame(vel_abs=names(tapply(ex$GDeff,ex$vel_abs,mean)),moy=tapply(ex$GDeff,ex$vel_abs,mean),sd=tapply(ex$GDeff,ex$vel_abs,sd),median=tapply(ex$GDeff,ex$vel_abs,median),q025=tapply(ex$GDeff,ex$vel_abs,quantile,probs=0.025),
               p975=tapply(ex$GDeff,ex$vel_abs,quantile,probs=0.975), pv.inf0=tapply(ex$GDeff,ex$vel_abs,test2,mu=0),pv.sup0=tapply(ex$GDeff,ex$vel_abs,test1,mu=0))
setwd(dir_out)
write.table(res,"summary_GDeff_allW.csv",sep=";",dec=".",row=F) #used to make table S6

#significant positive effect
subset(res,pv.sup0<0.05)
#significant negative effect
subset(res,pv.inf0<0.05)
#non-significant negative effect
subset(res,pv.inf0>=0.05 & pv.sup0>=0.05)

#CEW model
c1=read.csv2(paste("coeff_CEW.csv",sep=""),sep=";",dec=".",h=T)
dsel=subset(mydatatogo,Param=="O")
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)

no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

ex=foreach(i=1:nS,.combine=rbind,.inorder=F)%dopar%{
    gc(reset=T)
    print(i)
    c1a=subset(c1,nB==i)
    eff=c1a$Estimate[c1a$var=="scale(GD)"]+(c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]*((v1-mean(dsel$vel_abs))/sd(dsel$vel_abs)))
    tmp=data.frame(vel_abs=v1,GDeff=eff,nB=i)
    return(tmp)
}
parallel::stopCluster(cl = my.cluster2)

res=data.frame(vel_abs=names(tapply(ex$GDeff,ex$vel_abs,mean)),moy=tapply(ex$GDeff,ex$vel_abs,mean),sd=tapply(ex$GDeff,ex$vel_abs,sd),median=tapply(ex$GDeff,ex$vel_abs,median),q025=tapply(ex$GDeff,ex$vel_abs,quantile,probs=0.025),
               p975=tapply(ex$GDeff,ex$vel_abs,quantile,probs=0.975), pv.inf0=tapply(ex$GDeff,ex$vel_abs,test2,mu=0),pv.sup0=tapply(ex$GDeff,ex$vel_abs,test1,mu=0))
setwd(dir_out)
write.table(res,"summary_GDeff_CEW.csv",sep=";",dec=".",row=F) #used to make table S6

#significant positive effect
subset(res,pv.sup0<0.05)
#significant negative effect
subset(res,pv.inf0<0.05)
#non-significant negative effect
subset(res,pv.inf0>=0.05 & pv.sup0>=0.05)


#LEW model
c1=read.csv2(paste("coeff_LEW.csv",sep=""),sep=";",dec=".",h=T)
dsel=subset(mydatatogo,Param=="LE")
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)

no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

ex=foreach(i=1:nS,.combine=rbind,.inorder=F)%dopar%{
    gc(reset=T)
    print(i)
    c1a=subset(c1,nB==i)
    eff=c1a$Estimate[c1a$var=="scale(GD)"]+(c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]*((v1-mean(dsel$vel_abs))/sd(dsel$vel_abs)))
    tmp=data.frame(vel_abs=v1,GDeff=eff,nB=i)
    return(tmp)
}
parallel::stopCluster(cl = my.cluster2)

res=data.frame(vel_abs=names(tapply(ex$GDeff,ex$vel_abs,mean)),moy=tapply(ex$GDeff,ex$vel_abs,mean),sd=tapply(ex$GDeff,ex$vel_abs,sd),median=tapply(ex$GDeff,ex$vel_abs,median),q025=tapply(ex$GDeff,ex$vel_abs,quantile,probs=0.025),
               p975=tapply(ex$GDeff,ex$vel_abs,quantile,probs=0.975), pv.inf0=tapply(ex$GDeff,ex$vel_abs,test2,mu=0),pv.sup0=tapply(ex$GDeff,ex$vel_abs,test1,mu=0))
setwd(dir_out)
write.table(res,"summary_GDeff_LEW.csv",sep=";",dec=".",row=F) #used to make table S6

#significant positive effect
subset(res,pv.sup0<0.05)
#significant negative effect
subset(res,pv.inf0<0.05)
#non-significant negative effect
subset(res,pv.inf0>=0.05 & pv.sup0>=0.05)

#TEW model
c1=read.csv2(paste("coeff_TEW.csv",sep=""),sep=";",dec=".",h=T)
dsel=subset(mydatatogo,Param=="TE")
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)

no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

ex=foreach(i=1:nS,.combine=rbind,.inorder=F)%dopar%{
    gc(reset=T)
    print(i)
    c1a=subset(c1,nB==i)
    eff=c1a$Estimate[c1a$var=="scale(GD)"]+(c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]*((v1-mean(dsel$vel_abs))/sd(dsel$vel_abs)))
    tmp=data.frame(vel_abs=v1,GDeff=eff,nB=i)
    return(tmp)
}
parallel::stopCluster(cl = my.cluster2)

res=data.frame(vel_abs=names(tapply(ex$GDeff,ex$vel_abs,mean)),moy=tapply(ex$GDeff,ex$vel_abs,mean),sd=tapply(ex$GDeff,ex$vel_abs,sd),median=tapply(ex$GDeff,ex$vel_abs,median),q025=tapply(ex$GDeff,ex$vel_abs,quantile,probs=0.025),
               p975=tapply(ex$GDeff,ex$vel_abs,quantile,probs=0.975), pv.inf0=tapply(ex$GDeff,ex$vel_abs,test2,mu=0),pv.sup0=tapply(ex$GDeff,ex$vel_abs,test1,mu=0))
setwd(dir_out)
write.table(res,"summary_GDeff_TEW.csv",sep=";",dec=".",row=F) #used to make table S6

#significant positive effect
subset(res,pv.sup0<0.05)
#significant negative effect
subset(res,pv.inf0<0.05)
#non-significant negative effect
subset(res,pv.inf0>=0.05 & pv.sup0>=0.05)

#allW2 model
c1=read.csv2(paste("coeff_allW2.csv",sep=""),sep=";",dec=".",h=T)
dsel=mydatatogo
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)

no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

ex=foreach(i=1:nS,.combine=rbind,.inorder=F)%dopar%{
    gc(reset=T)
    print(i)
    c1a=subset(c1,nB==i)
    eff=c1a$Estimate[c1a$var=="scale(GD)"]+(c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]*((v1-mean(dsel$vel_abs))/sd(dsel$vel_abs)))
    tmp=data.frame(vel_abs=v1,GDeff_CE=eff)
    eff=c1a$Estimate[c1a$var=="scale(GD)"]+c1a$Estimate[c1a$var=="scale(GD):ParamLE"]+((c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]+c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD):ParamLE"])*((v1-mean(dsel$vel_abs))/sd(dsel$vel_abs)))
    tmp$GDeff_LE=eff
    eff=c1a$Estimate[c1a$var=="scale(GD)"]+c1a$Estimate[c1a$var=="scale(GD):ParamTE"]+((c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]+c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD):ParamTE"])*((v1-mean(dsel$vel_abs))/sd(dsel$vel_abs)))
    tmp$GDeff_TE=eff
    tmp$nB=i
    return(tmp)
}
parallel::stopCluster(cl = my.cluster2)

#summary statistictics
resCE=data.frame(vel_abs=names(tapply(ex$GDeff_CE,ex$vel_abs,mean)),moy=tapply(ex$GDeff_CE,ex$vel_abs,mean),sd=tapply(ex$GDeff_CE,ex$vel_abs,sd),median=tapply(ex$GDeff_CE,ex$vel_abs,median),q025=tapply(ex$GDeff_CE,ex$vel_abs,quantile,probs=0.025),
                 p975=tapply(ex$GDeff_CE,ex$vel_abs,quantile,probs=0.975), pv.inf0=tapply(ex$GDeff_CE,ex$vel_abs,test2,mu=0),pv.sup0=tapply(ex$GDeff_CE,ex$vel_abs,test1,mu=0))
resLE=data.frame(vel_abs=names(tapply(ex$GDeff_LE,ex$vel_abs,mean)),moy=tapply(ex$GDeff_LE,ex$vel_abs,mean),sd=tapply(ex$GDeff_LE,ex$vel_abs,sd),median=tapply(ex$GDeff_LE,ex$vel_abs,median),q025=tapply(ex$GDeff_LE,ex$vel_abs,quantile,probs=0.025),
                 p975=tapply(ex$GDeff_LE,ex$vel_abs,quantile,probs=0.975), pv.inf0=tapply(ex$GDeff_LE,ex$vel_abs,test2,mu=0),pv.sup0=tapply(ex$GDeff_LE,ex$vel_abs,test1,mu=0))
resTE=data.frame(vel_abs=names(tapply(ex$GDeff_TE,ex$vel_abs,mean)),moy=tapply(ex$GDeff_TE,ex$vel_abs,mean),sd=tapply(ex$GDeff_TE,ex$vel_abs,sd),median=tapply(ex$GDeff_TE,ex$vel_abs,median),q025=tapply(ex$GDeff_TE,ex$vel_abs,quantile,probs=0.025),
                 p975=tapply(ex$GDeff_TE,ex$vel_abs,quantile,probs=0.975), pv.inf0=tapply(ex$GDeff_TE,ex$vel_abs,test2,mu=0),pv.sup0=tapply(ex$GDeff_TE,ex$vel_abs,test1,mu=0))


setwd(dir_out)
write.table(resTE,"summary_GDeff_allW2_TE.csv",sep=";",dec=".",row=F) #used to make table S6
write.table(resCE,"summary_GDeff_allW2_CE.csv",sep=";",dec=".",row=F) #used to make table S6
write.table(resLE,"summary_GDeff_allW2_LE.csv",sep=";",dec=".",row=F) #used to make table S6

#significant positive effect
subset(resCE,pv.sup0<0.05)
#significant negative effect
subset(resCE,pv.inf0<0.05)
#non-significant negative effect
subset(resCE,pv.inf0>=0.05 & pv.sup0>=0.05)


#significant positive effect
subset(resLE,pv.sup0<0.05)
#significant negative effect
subset(resLE,pv.inf0<0.05)
#non-significant negative effect
subset(resLE,pv.inf0>=0.05 & pv.sup0>=0.05)


#significant positive effect
subset(resTE,pv.sup0<0.05)
#significant negative effect
subset(resTE,pv.inf0<0.05)
#non-significant negative effect
subset(resTE,pv.inf0>=0.05 & pv.sup0>=0.05)


############################################################################################
####testing genetic diversity value for which the climate velocity effect is significant####
############################################################################################

setwd(dir_out)
nS=5000#number of boostrap to use to infer stat

#allW model
c1=read.csv2(paste("coeff_allW.csv",sep=""),sep=";",dec=".",h=T)
dsel=mydatatogo
v1=seq(round(min(dsel$GD),2),round(max(dsel$GD),2),by=0.01)

no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

ex=foreach(i=1:nS,.combine=rbind,.inorder=F)%dopar%{
    gc(reset=T)
    print(i)
    c1a=subset(c1,nB==i)
    eff=c1a$Estimate[c1a$var=="scale(vel_abs)"]+(c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]*((v1-mean(dsel$GD))/sd(dsel$GD)))
    tmp=data.frame(GD=v1,VAeff=eff,nB=i)
    return(tmp)
}
parallel::stopCluster(cl = my.cluster2)

res=data.frame(GD=names(tapply(ex$VAeff,ex$GD,mean)),moy=tapply(ex$VAeff,ex$GD,mean),sd=tapply(ex$VAeff,ex$GD,sd),median=tapply(ex$VAeff,ex$GD,median),q025=tapply(ex$VAeff,ex$GD,quantile,probs=0.025),
               p975=tapply(ex$VAeff,ex$GD,quantile,probs=0.975), pv.inf0=tapply(ex$VAeff,ex$GD,test2,mu=0),pv.sup0=tapply(ex$VAeff,ex$GD,test1,mu=0))
setwd(dir_out)
write.table(res,"summary_velabs_allW.csv",sep=";",dec=".",row=F) #used to make table S6

#significant positive effect
subset(res,pv.sup0<0.05)
#significant negative effect
subset(res,pv.inf0<0.05)
#non-significant negative effect
subset(res,pv.inf0>=0.05 & pv.sup0>=0.05)

#CEW model
c1=read.csv2(paste("coeff_CEW.csv",sep=""),sep=";",dec=".",h=T)
dsel=subset(mydatatogo,Param=="O")
v1=seq(round(min(dsel$GD),2),round(max(dsel$GD),2),by=0.01)

no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

ex=foreach(i=1:nS,.combine=rbind,.inorder=F)%dopar%{
    gc(reset=T)
    print(i)
    c1a=subset(c1,nB==i)
    eff=c1a$Estimate[c1a$var=="scale(vel_abs)"]+(c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]*((v1-mean(dsel$GD))/sd(dsel$GD)))
    tmp=data.frame(GD=v1,VAeff=eff,nB=i)
    return(tmp)
}
parallel::stopCluster(cl = my.cluster2)

res=data.frame(GD=names(tapply(ex$VAeff,ex$GD,mean)),moy=tapply(ex$VAeff,ex$GD,mean),sd=tapply(ex$VAeff,ex$GD,sd),median=tapply(ex$VAeff,ex$GD,median),q025=tapply(ex$VAeff,ex$GD,quantile,probs=0.025),
               p975=tapply(ex$VAeff,ex$GD,quantile,probs=0.975), pv.inf0=tapply(ex$VAeff,ex$GD,test2,mu=0),pv.sup0=tapply(ex$VAeff,ex$GD,test1,mu=0))
setwd(dir_out)
write.table(res,"summary_velabs_CEW.csv",sep=";",dec=".",row=F) #used to make table S6

#significant positive effect
subset(res,pv.sup0<0.05)
#significant negative effect
subset(res,pv.inf0<0.05)
#non-significant negative effect
subset(res,pv.inf0>=0.05 & pv.sup0>=0.05)


#LEW model
c1=read.csv2(paste("coeff_LEW.csv",sep=""),sep=";",dec=".",h=T)
dsel=subset(mydatatogo,Param=="LE")
v1=seq(round(min(dsel$GD),2),round(max(dsel$GD),2),by=0.01)

no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

ex=foreach(i=1:nS,.combine=rbind,.inorder=F)%dopar%{
    gc(reset=T)
    print(i)
    c1a=subset(c1,nB==i)
    eff=c1a$Estimate[c1a$var=="scale(vel_abs)"]+(c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]*((v1-mean(dsel$GD))/sd(dsel$GD)))
    tmp=data.frame(GD=v1,VAeff=eff,nB=i)
    return(tmp)
}
parallel::stopCluster(cl = my.cluster2)

res=data.frame(GD=names(tapply(ex$VAeff,ex$GD,mean)),moy=tapply(ex$VAeff,ex$GD,mean),sd=tapply(ex$VAeff,ex$GD,sd),median=tapply(ex$VAeff,ex$GD,median),q025=tapply(ex$VAeff,ex$GD,quantile,probs=0.025),
               p975=tapply(ex$VAeff,ex$GD,quantile,probs=0.975), pv.inf0=tapply(ex$VAeff,ex$GD,test2,mu=0),pv.sup0=tapply(ex$VAeff,ex$GD,test1,mu=0))
setwd(dir_out)
write.table(res,"summary_velabs_LEW.csv",sep=";",dec=".",row=F) #used to make table S6

#significant positive effect
subset(res,pv.sup0<0.05)
#significant negative effect
subset(res,pv.inf0<0.05)
#non-significant negative effect
subset(res,pv.inf0>=0.05 & pv.sup0>=0.05)

#TEW model
c1=read.csv2(paste("coeff_TEW.csv",sep=""),sep=";",dec=".",h=T)
dsel=subset(mydatatogo,Param=="TE")
v1=seq(round(min(dsel$GD),2),round(max(dsel$GD),2),by=0.01)

no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

ex=foreach(i=1:nS,.combine=rbind,.inorder=F)%dopar%{
    gc(reset=T)
    print(i)
    c1a=subset(c1,nB==i)
    eff=c1a$Estimate[c1a$var=="scale(vel_abs)"]+(c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]*((v1-mean(dsel$GD))/sd(dsel$GD)))
    tmp=data.frame(GD=v1,VAeff=eff,nB=i)
    return(tmp)
}
parallel::stopCluster(cl = my.cluster2)

res=data.frame(GD=names(tapply(ex$VAeff,ex$GD,mean)),moy=tapply(ex$VAeff,ex$GD,mean),sd=tapply(ex$VAeff,ex$GD,sd),median=tapply(ex$VAeff,ex$GD,median),q025=tapply(ex$VAeff,ex$GD,quantile,probs=0.025),
               p975=tapply(ex$VAeff,ex$GD,quantile,probs=0.975), pv.inf0=tapply(ex$VAeff,ex$GD,test2,mu=0),pv.sup0=tapply(ex$VAeff,ex$GD,test1,mu=0))
setwd(dir_out)
write.table(res,"summary_velabs_TEW.csv",sep=";",dec=".",row=F) #used to make table S6

#significant positive effect
subset(res,pv.sup0<0.05)
#significant negative effect
subset(res,pv.inf0<0.05)
#non-significant negative effect
subset(res,pv.inf0>=0.05 & pv.sup0>=0.05)

#allW2 model
c1=read.csv2(paste("coeff_allW2.csv",sep=""),sep=";",dec=".",h=T)
dsel=mydatatogo
v1=seq(round(min(dsel$GD),3),round(max(dsel$GD),3),by=0.001)

no=8 #number of cpu not used for paralell computing
parallel::detectCores() #maximum CPU avilable for paralell computing
nCPU=parallel::detectCores()-no
#nCPU=12
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster2)

ex=foreach(i=1:nS,.combine=rbind,.inorder=F)%dopar%{
    gc(reset=T)
    print(i)
    c1a=subset(c1,nB==i)
    eff=c1a$Estimate[c1a$var=="scale(vel_abs)"]+(c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]*((v1-mean(dsel$GD))/sd(dsel$GD)))
    tmp=data.frame(GD=v1,VAeff_CE=eff)
    eff=c1a$Estimate[c1a$var=="scale(vel_abs)"]+c1a$Estimate[c1a$var=="scale(vel_abs):ParamLE"]+((c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]+c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD):ParamLE"])*((v1-mean(dsel$GD))/sd(dsel$GD)))
    tmp$VAeff_LE=eff
    eff=c1a$Estimate[c1a$var=="scale(vel_abs)"]+c1a$Estimate[c1a$var=="scale(vel_abs):ParamTE"]+((c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD)"]+c1a$Estimate[c1a$var=="scale(vel_abs):scale(GD):ParamTE"])*((v1-mean(dsel$GD))/sd(dsel$GD)))
    tmp$VAeff_TE=eff
    tmp$nB=i
    return(tmp)
}
parallel::stopCluster(cl = my.cluster2)

#summary statistictics
resCE=data.frame(GD=names(tapply(ex$VAeff_CE,ex$GD,mean)),moy=tapply(ex$VAeff_CE,ex$GD,mean),sd=tapply(ex$VAeff_CE,ex$GD,sd),median=tapply(ex$VAeff_CE,ex$GD,median),q025=tapply(ex$VAeff_CE,ex$GD,quantile,probs=0.025),
                 p975=tapply(ex$VAeff_CE,ex$GD,quantile,probs=0.975), pv.inf0=tapply(ex$VAeff_CE,ex$GD,test2,mu=0),pv.sup0=tapply(ex$VAeff_CE,ex$GD,test1,mu=0))
resLE=data.frame(GD=names(tapply(ex$VAeff_LE,ex$GD,mean)),moy=tapply(ex$VAeff_LE,ex$GD,mean),sd=tapply(ex$VAeff_LE,ex$GD,sd),median=tapply(ex$VAeff_LE,ex$GD,median),q025=tapply(ex$VAeff_LE,ex$GD,quantile,probs=0.025),
                 p975=tapply(ex$VAeff_LE,ex$GD,quantile,probs=0.975), pv.inf0=tapply(ex$VAeff_LE,ex$GD,test2,mu=0),pv.sup0=tapply(ex$VAeff_LE,ex$GD,test1,mu=0))
resTE=data.frame(GD=names(tapply(ex$VAeff_TE,ex$GD,mean)),moy=tapply(ex$VAeff_TE,ex$GD,mean),sd=tapply(ex$VAeff_TE,ex$GD,sd),median=tapply(ex$VAeff_TE,ex$GD,median),q025=tapply(ex$VAeff_TE,ex$GD,quantile,probs=0.025),
                 p975=tapply(ex$VAeff_TE,ex$GD,quantile,probs=0.975), pv.inf0=tapply(ex$VAeff_TE,ex$GD,test2,mu=0),pv.sup0=tapply(ex$VAeff_TE,ex$GD,test1,mu=0))


setwd(dir_out)
write.table(resTE,"summary_VAeff_allW2_TE.csv",sep=";",dec=".",row=F) #used to make table S6
write.table(resLE,"summary_VAeff_allW2_LE.csv",sep=";",dec=".",row=F) #used to make table S6
write.table(resCE,"summary_VAeff_allW2_CE.csv",sep=";",dec=".",row=F) #used to make table S6

#significant positive effect
subset(resCE,pv.sup0<0.05)
#significant negative effect
subset(resCE,pv.inf0<0.05)
#non-significant negative effect
subset(resCE,pv.inf0>=0.05 & pv.sup0>=0.05)


#significant positive effect
subset(resLE,pv.sup0<0.05)
#significant negative effect
subset(resLE,pv.inf0<0.05)
#non-significant negative effect
subset(resLE,pv.inf0>=0.05 & pv.sup0>=0.05)


#significant positive effect
subset(resTE,pv.sup0<0.05)
#significant negative effect
subset(resTE,pv.inf0<0.05)
#non-significant negative effect
subset(resTE,pv.inf0>=0.05 & pv.sup0>=0.05)

################################################################################
#########graphics###############################################################
################################################################################

#from allW2 model output and stat
#link function for gamma(log link) model => exp(coef[1] + coef[2] * X)
#what happens at TE?
setwd(dir_out)
res=read.csv2("summary_GDeff_allW2_TE.csv",sep=";",dec=".",h=T) 
c1=read.csv2(paste("summary_coeff.csv",sep=""),sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW2")
int=c1a$moy[c1a$var=="(Intercept)"]+c1a$moy[c1a$var=="ParamTE"]
dsel=subset(mydatatogo,Param=="TE")
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),by=0.0001)

for(i in 1:nrow(res)){
    #pred1=exp(int+res$moy[i]*((v2-mean(mydatatogo$GD))/sd(mydatatogo$GD)))
    pred1=exp(int+res$moy[i]*((v2-mean(mydatatogo$GD))/sd(mydatatogo$GD))+c1a$moy[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$moy[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)) #adding effects of covariables did not change the relationships, but it changes the predicted value of range shift. Furthermore, the predictive value are representative of the mean number of temporal units and the mean spatial extent found in the database. The predictive values are also reprsentative of a combination of method (that define the intercept values of the model): Presence data and balanced sampling. Of course, the predictions don not consider the random effect  
    pred1=data.frame(GD=v2,pred1,vel_abs=res$vel_abs[i])
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}
#png("/media/rom1/stock/PhyloRangeShift/CESAB_proj/meeting/meet2_1122/traitPerspective/fig1/v042023/2_LE_DDmax_velo_curve.png",unit="cm",width=15,height=15*0.66,res=300)#,width=547,height=360
#resX=subset(res,pv.inf0<0.05 | pv.sup0<0.05)
res$ID=as.character(round(res$vel_abs,1))
pred2$ID=as.character(round(pred2$vel_abs,1))
resX=subset(res,pv.inf0<0.05)
if(nrow(resX>0)){
    selN=round(seq(min(resX$vel_abs),max(resX$vel_abs),le=5),1)
    predX=merge(pred2,data.frame(ID=as.character(selN),signifN=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX=pred2
    predX$signifN=NA
}

resX=subset(res,pv.sup0<0.05)
if(nrow(resX>0)){
    selP=round(seq(min(resX$vel_abs),max(resX$vel_abs),le=5),1)
    predXTE=merge(predX,data.frame(ID=as.character(selP),signifP=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predXTE=predX
    predXTE$signifP=NA
}

predXTE1=subset(predXTE,signifN==1 | signifP==1)
gg= ggplot(predXTE, aes(GD, pred1, color = vel_abs, group=ID)) +
    geom_line(data=predXTE1, aes(GD, pred1, color = vel_abs, group=ID)) +
    theme_minimal() +
    #labs(x=bquote("Nucleotide diversity"), y=bquote("Absolute velocity of range shifts (km.yr<sup>-1</sup>)"))+
    labs(x=bquote(""), y=bquote(""))+
    scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"),
                          name = bquote("Absolute velocity of isotherm shifts (km.yr<sup>-1</sup>)"), 
                          guide = guide_colourbar(title.position = "right",barwidth = 1,barheight = 15),
                          limits=c(0, 7),
                          breaks=c(seq(0,7,by=1)))

gg <- gg + theme_classic(base_size=11,base_family="Helvetica")
gg <- gg + theme(plot.title=element_text(hjust=0.5))
#gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + coord_cartesian(clip="off")
gg <- gg + scale_x_continuous(breaks=c(c(seq(0,0.05,by=.01))),labels=c("0",as.character(seq(0.01,0.05,by=.01))),limits=c(-0,0.055),expand=c(0,0))#expand=expand_scale(mult = c(0, 0),add = c(0.2, 0)))
gg <- gg + scale_y_continuous(breaks=c(seq(0,6,by=2)),limits=c(0,6),expand=c(0,0))
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(axis.text.x=element_text(vjust=1))
gg <- gg + theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
gg <- gg + theme(axis.title.x = element_markdown(margin = margin(t = 19,b=0),size=11))
gg <- gg + theme(axis.title.y = element_markdown(margin = margin(r =19),size=11))
gg <- gg + theme(legend.title = element_markdown(angle = -90),legend.title.align = 0.5)
#gg <- gg + theme(axis.title.y = element_text(margin = margin(r =10),size=11))
#gg <- gg + geom_abline(slope=1, intercept=0,linetype="dashed",expand=F)
#gg <- gg + geom_segment(x = -0.2, y = -0.2, xend = 8, yend = 8, linetype="dashed",colour="black")
#gg <- gg + annotate(geom="text", x=-0.75, y=7.5, label="(a)",size=5)
gg <- gg + labs(tag = '(a)')
gg <- gg + theme(plot.tag.position = c(-0.01, 0.99),plot.tag = element_text(hjust = 0))
#gg <- gg + theme(legend.position = "none")
gg1<-gg
gg1
#dev.off()

#what happens at LE?
setwd(dir_out)
res=read.csv2("summary_GDeff_allW2_LE.csv",sep=";",dec=".",h=T) 
c1=read.csv2(paste("summary_coeff.csv",sep=""),sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW2")
int=c1a$moy[c1a$var=="(Intercept)"]+c1a$moy[c1a$var=="ParamLE"]
dsel=subset(mydatatogo,Param=="LE")
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),by=0.0001)

for(i in 1:nrow(res)){
    #pred1=exp(int+res$moy[i]*((v2-mean(mydatatogo$GD))/sd(mydatatogo$GD)))
    pred1=exp(int+res$moy[i]*((v2-mean(mydatatogo$GD))/sd(mydatatogo$GD))+c1a$moy[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$moy[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)) #adding effects of covariables did not change the relationships, but it changes the predicted value of range shift. Furthermore, the predictive value are representative of the mean number of temporal units and the mean spatial extent found in the database. The predictive values are also reprsentative of a combination of method (that define the intercept values of the model): Presence data and balanced sampling. Of course, the predictions don not consider the random effect  
    pred1=data.frame(GD=v2,pred1,vel_abs=res$vel_abs[i])
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}
#png("/media/rom1/stock/PhyloRangeShift/CESAB_proj/meeting/meet2_1122/traitPerspective/fig1/v042023/2_LE_DDmax_velo_curve.png",unit="cm",width=15,height=15*0.66,res=300)#,width=547,height=360
#resX=subset(res,pv.inf0<0.05 | pv.sup0<0.05)
res$ID=as.character(round(res$vel_abs,1))
pred2$ID=as.character(round(pred2$vel_abs,1))
resX=subset(res,pv.inf0<0.05)
if(nrow(resX>0)){
    selN=round(seq(min(resX$vel_abs),max(resX$vel_abs),le=5),1)
    predX=merge(pred2,data.frame(ID=as.character(selN),signifN=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX=pred2
    predX$signifN=NA
}

resX=subset(res,pv.sup0<0.05)
if(nrow(resX>0)){
    selP=round(seq(min(resX$vel_abs),max(resX$vel_abs),le=5),1)
    predXLE=merge(predX,data.frame(ID=as.character(selP),signifP=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predXLE=predX
    predXLE$signifP=NA
}

predXLE1=subset(predXLE,signifN==1 | signifP==1)
gg= ggplot(predXLE, aes(GD, pred1, color = vel_abs, group=ID)) +
    geom_line(data=predXLE1, aes(GD, pred1, color = vel_abs, group=ID)) +
    theme_minimal() +
    #labs(x=bquote("Nucleotide diversity"), y=bquote("Absolute velocity of range shifts (km.yr<sup>-1</sup>)"))+
    labs(x=bquote(""), y=bquote(""))+
    scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"),
                          name = bquote("Absolute velocity of isotherm shifts (km.yr<sup>-1</sup>)"), 
                          guide = guide_colourbar(title.position = "right",barwidth = 1,barheight = 15),
                          limits=c(0, 7),
                          breaks=c(seq(0,7,by=1)))

gg <- gg + theme_classic(base_size=11,base_family="Helvetica")
gg <- gg + theme(plot.title=element_text(hjust=0.5))
#gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + coord_cartesian(clip="off")
gg <- gg + scale_x_continuous(breaks=c(c(seq(0,0.05,by=.01))),labels=c("0",as.character(seq(0.01,0.05,by=.01))),limits=c(-0,0.055),expand=c(0,0))#expand=expand_scale(mult = c(0, 0),add = c(0.2, 0)))
gg <- gg + scale_y_continuous(breaks=c(seq(0,6,by=2)),limits=c(0,6),expand=c(0,0))
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(axis.text.x=element_text(vjust=1))
gg <- gg + theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
gg <- gg + theme(axis.title.x = element_markdown(margin = margin(t = 19,b=0),size=11))
gg <- gg + theme(axis.title.y = element_markdown(margin = margin(r =19),size=11))
gg <- gg + theme(legend.title = element_markdown(angle = -90),legend.title.align = 0.5)
#gg <- gg + theme(axis.title.y = element_text(margin = margin(r =10),size=11))
#gg <- gg + geom_abline(slope=1, intercept=0,linetype="dashed",expand=F)
#gg <- gg + geom_segment(x = -0.2, y = -0.2, xend = 8, yend = 8, linetype="dashed",colour="black")
#gg <- gg + annotate(geom="text", x=-0.75, y=7.5, label="(a)",size=5)
gg <- gg + labs(tag = '(c)')
gg <- gg + theme(plot.tag.position = c(-0.01, 0.99),plot.tag = element_text(hjust = 0))
#gg <- gg + theme(legend.position = "none")
gg3<-gg
gg3
#dev.off()
#what happens at CE?
setwd(dir_out)
res=read.csv2("summary_GDeff_allW2_CE.csv",sep=";",dec=".",h=T) 
c1=read.csv2(paste("summary_coeff.csv",sep=""),sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW2")
int=c1a$moy[c1a$var=="(Intercept)"]
dsel=subset(mydatatogo,Param=="O")
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),by=0.0001)

for(i in 1:nrow(res)){
    #pred1=exp(int+res$moy[i]*((v2-mean(mydatatogo$GD))/sd(mydatatogo$GD)))
    pred1=exp(int+res$moy[i]*((v2-mean(mydatatogo$GD))/sd(mydatatogo$GD))+c1a$moy[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$moy[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)) #adding effects of covariables did not change the relationships, but it changes the predicted value of range shift. Furthermore, the predictive value are representative of the mean number of temporal units and the mean spatial extent found in the database. The predictive values are also reprsentative of a combination of method (that define the intercept values of the model): Presence data and balanced sampling. Of course, the predictions don not consider the random effect  
    pred1=data.frame(GD=v2,pred1,vel_abs=res$vel_abs[i])
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}


#png("/media/rom1/stock/PhyloRangeShift/CESAB_proj/meeting/meet2_1122/traitPerspective/fig1/v042023/2_LE_DDmax_velo_curve.png",unit="cm",width=15,height=15*0.66,res=300)#,width=547,height=360
#resX=subset(res,pv.inf0<0.05 | pv.sup0<0.05)
res$ID=as.character(round(res$vel_abs,1))
pred2$ID=as.character(round(pred2$vel_abs,1))
resX=subset(res,pv.inf0<0.05)
if(nrow(resX>0)){
    selN=round(seq(min(resX$vel_abs),max(resX$vel_abs),le=5),1)
    predX=merge(pred2,data.frame(ID=as.character(selN),signifN=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX=pred2
    predX$signifN=NA
}

resX=subset(res,pv.sup0<0.05)
if(nrow(resX>0)){
    selP=round(seq(min(resX$vel_abs),max(resX$vel_abs),le=5),1)
    predXCE=merge(predX,data.frame(ID=as.character(selP),signifP=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predXCE=predX
    predXCE$signifP=NA
}

predXCE1=subset(predXCE,signifN==1 | signifP==1)
gg= ggplot(predXCE, aes(GD, pred1, color = vel_abs, group=ID)) +
    geom_line(data=predXCE1, aes(GD, pred1, color = vel_abs, group=ID)) +
    theme_minimal() +
    #labs(x=bquote("Nucleotide diversity"), y=bquote("Absolute velocity of range shifts (km.yr<sup>-1</sup>)"))+
    labs(x=bquote(""), y=bquote(""))+
    scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"),
                          name = bquote("Absolute velocity of isotherm shifts (km.yr<sup>-1</sup>)"), 
                          guide = guide_colourbar(title.position = "right",barwidth = 1,barheight = 15),
                          limits=c(0, 7),
                          breaks=c(seq(0,7,by=1)))

gg <- gg + theme_classic(base_size=11,base_family="Helvetica")
gg <- gg + theme(plot.title=element_text(hjust=0.5))
#gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + coord_cartesian(clip="off")
gg <- gg + scale_x_continuous(breaks=c(c(seq(0,0.05,by=.01))),labels=c("0",as.character(seq(0.01,0.05,by=.01))),limits=c(-0,0.055),expand=c(0,0))#expand=expand_scale(mult = c(0, 0),add = c(0.2, 0)))
gg <- gg + scale_y_continuous(breaks=c(seq(0,6,by=2)),limits=c(0,6),expand=c(0,0))
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(axis.text.x=element_text(vjust=1))
gg <- gg + theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
gg <- gg + theme(axis.title.x = element_markdown(margin = margin(t = 19,b=0),size=11))
gg <- gg + theme(axis.title.y = element_markdown(margin = margin(r =19),size=11))
gg <- gg + theme(legend.title = element_markdown(angle = -90),legend.title.align = 0.5)
#gg <- gg + theme(axis.title.y = element_text(margin = margin(r =10),size=11))
#gg <- gg + geom_abline(slope=1, intercept=0,linetype="dashed",expand=F)
#gg <- gg + geom_segment(x = -0.2, y = -0.2, xend = 8, yend = 8, linetype="dashed",colour="black")
#gg <- gg + annotate(geom="text", x=-0.75, y=7.5, label="(a)",size=5)
gg <- gg + labs(tag = '(b)')
gg <- gg + theme(plot.tag.position = c(-0.01, 0.99),plot.tag = element_text(hjust = 0))
#gg <- gg + theme(legend.position = "none")
gg2<-gg
gg2
#dev.off()



png(paste0(dir_out,"/VICeffectOnGD_v0.png"),unit="cm",width=27,height=11,res=300)#,width=547,height=360
#gs <- lapply(1:4, function(ii) 
#grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
lay <- rbind(c(1,2,3))
grid.arrange(gg1,gg2,gg3, 
             ncol=3, nrow=1, widths=c(1,1,1), heights=c(1),layout_matrix = lay)
dev.off()

################################################################################
###!!!! important from Romain: need to think more about how I predict range shift with inverse link functions
### because the range shift is low compared to what Brunno predicted in Rpub figures.

#from allW2 model output and stat
#link function for gamma(log link) model => exp(coef[1] + coef[2] * X)
#what happens at TE?
setwd(dir_out)
res=read.csv2("summary_VAeff_allW2_TE.csv",sep=";",dec=".",h=T) 
c1=read.csv2(paste("summary_coeff.csv",sep=""),sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW2")
int=c1a$moy[c1a$var=="(Intercept)"]+c1a$moy[c1a$var=="ParamTE"]
dsel=subset(mydatatogo,Param=="TE")
v2=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
#v2=seq(round(min(dsel$GD),3),round(max(dsel$GD),3),by=0.0001)

for(i in 1:nrow(res)){
    #pred1=exp(int+res$moy[i]*((v2-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs)))
    pred1=exp(int+res$moy[i]*((v2-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))+c1a$moy[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$moy[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)) #adding effects of covariables did not change the relationships, but it changes the predicted value of range shift. Furthermore, the predictive value are representative of the mean number of temporal units and the mean spatial extent found in the database. The predictive values are also reprsentative of a combination of method (that define the intercept values of the model): Presence data and balanced sampling. Of course, the predictions don not consider the random effect  
    pred1=data.frame(GD=res$GD[i],pred1,vel_abs=v2)
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}
#png("/media/rom1/stock/PhyloRangeShift/CESAB_proj/meeting/meet2_1122/traitPerspective/fig1/v042023/2_LE_DDmax_velo_curve.png",unit="cm",width=15,height=15*0.66,res=300)#,width=547,height=360
#resX=subset(res,pv.inf0<0.05 | pv.sup0<0.05)
res$ID=as.character(round(res$GD,3))
pred2$ID=as.character(round(pred2$GD,3))
resX=subset(res,pv.inf0<0.05)
if(nrow(resX>0)){
    if(nrow(resX)<5){
        selN=round(seq(min(resX$GD),max(resX$GD),le=nrow(resX)),3)
    }else{
        selN=round(seq(min(resX$GD),max(resX$GD),le=5),3)
    }
    predX=merge(pred2,data.frame(ID=as.character(selN),signifN=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX=pred2
    predX$signifN=NA
}

resX=subset(res,pv.sup0<0.05)
if(nrow(resX>0)){
    if(nrow(resX)<5){
        selP=round(seq(min(resX$GD),max(resX$GD),le=nrow(resX)),3)
    }else{
        selP=round(seq(min(resX$GD),max(resX$GD),le=5),3)
    }
    predXTE=merge(predX,data.frame(ID=as.character(selP),signifP=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predXTE=predX
    predXTE$signifP=NA
}

predXTE1=subset(predXTE,signifN==1 | signifP==1)
gg= ggplot(predXTE, aes(vel_abs, pred1, color = vel_abs, group=ID)) +
    geom_line(data=predXTE1, aes(vel_abs, pred1, color = GD, group=ID)) +
    theme_minimal() +
    #labs(x=bquote("Nucleotide diversity"), y=bquote("Absolute velocity of range shifts (km.yr<sup>-1</sup>)"))+
    labs(x=bquote(""), y=bquote(""))+
    scale_color_viridis(name = bquote("Genetic diversity"), 
                        guide = guide_colourbar(title.position = "right",barwidth = 1,barheight = 15),
                        limits=c(0, 0.055),
                        breaks=c(seq(0,0.05,by=0.01)))

gg <- gg + theme_classic(base_size=11,base_family="Helvetica")
gg <- gg + theme(plot.title=element_text(hjust=0.5))
#gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + coord_cartesian(clip="off")
gg <- gg + scale_x_continuous(breaks=c(c(seq(0,7,by=1))),labels=c("0",as.character(seq(1,7,by=1))),limits=c(0,7),expand=c(0,0))#expand=expand_scale(mult = c(0, 0),add = c(0.2, 0)))
gg <- gg + scale_y_continuous(breaks=c(seq(0,10,by=2)),limits=c(0,10),expand=c(0,0))
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(axis.text.x=element_text(vjust=1))
gg <- gg + theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
gg <- gg + theme(axis.title.x = element_markdown(margin = margin(t = 19,b=0),size=11))
gg <- gg + theme(axis.title.y = element_markdown(margin = margin(r =19),size=11))
gg <- gg + theme(legend.title = element_markdown(angle = -90),legend.title.align = 0.5)
#gg <- gg + theme(axis.title.y = element_text(margin = margin(r =10),size=11))
#gg <- gg + geom_abline(slope=1, intercept=0,linetype="dashed",expand=F)
#gg <- gg + geom_segment(x = -0.2, y = -0.2, xend = 8, yend = 8, linetype="dashed",colour="black")
#gg <- gg + annotate(geom="text", x=-0.75, y=7.5, label="(a)",size=5)
gg <- gg + labs(tag = '(a)')
gg <- gg + theme(plot.tag.position = c(-0.01, 0.99),plot.tag = element_text(hjust = 0))
#gg <- gg + theme(legend.position = "none")
gg1<-gg
gg1
#dev.off()

#what happens at LE?
setwd(dir_out)
res=read.csv2("summary_VAeff_allW2_LE.csv",sep=";",dec=".",h=T) 
c1=read.csv2(paste("summary_coeff.csv",sep=""),sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW2")
int=c1a$moy[c1a$var=="(Intercept)"]+c1a$moy[c1a$var=="ParamLE"]
dsel=subset(mydatatogo,Param=="LE")
v2=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
#v2=seq(round(min(dsel$GD),3),round(max(dsel$GD),3),by=0.0001)

for(i in 1:nrow(res)){
    #pred1=exp(int+res$moy[i]*((v2-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs)))
    pred1=exp(int+res$moy[i]*((v2-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))+c1a$moy[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$moy[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)) #adding effects of covariables did not change the relationships, but it changes the predicted value of range shift. Furthermore, the predictive value are representative of the mean number of temporal units and the mean spatial extent found in the database. The predictive values are also reprsentative of a combination of method (that define the intercept values of the model): Presence data and balanced sampling. Of course, the predictions don not consider the random effect  
    pred1=data.frame(GD=res$GD[i],pred1,vel_abs=v2)
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}
#png("/media/rom1/stock/PhyloRangeShift/CESAB_proj/meeting/meet2_1122/traitPerspective/fig1/v042023/2_LE_DDmax_velo_curve.png",unit="cm",width=15,height=15*0.66,res=300)#,width=547,height=360
#resX=subset(res,pv.inf0<0.05 | pv.sup0<0.05)
res$ID=as.character(round(res$GD,3))
pred2$ID=as.character(round(pred2$GD,3))
resX=subset(res,pv.inf0<0.05)
if(nrow(resX>0)){
    if(nrow(resX)<5){
        selN=round(seq(min(resX$GD),max(resX$GD),le=nrow(resX)),3)
    }else{
        selN=round(seq(min(resX$GD),max(resX$GD),le=5),3)
    }
    predX=merge(pred2,data.frame(ID=as.character(selN),signifN=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX=pred2
    predX$signifN=NA
}

resX=subset(res,pv.sup0<0.05)
if(nrow(resX>0)){
    if(nrow(resX)<5){
        selP=round(seq(min(resX$GD),max(resX$GD),le=nrow(resX)),3)
    }else{
        selP=round(seq(min(resX$GD),max(resX$GD),le=5),3)
    }
    predXLE=merge(predX,data.frame(ID=as.character(selP),signifP=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predXLE=predX
    predXLE$signifP=NA
}

predXLE1=subset(predXLE,signifN==1 | signifP==1)
gg= ggplot(predXLE, aes(vel_abs, pred1, color = vel_abs, group=ID)) +
    geom_line(data=predXLE1, aes(vel_abs, pred1, color = GD, group=ID)) +
    theme_minimal() +
    #labs(x=bquote("Nucleotide diversity"), y=bquote("Absolute velocity of range shifts (km.yr<sup>-1</sup>)"))+
    labs(x=bquote(""), y=bquote(""))+
    scale_color_viridis(name = bquote("Genetic diversity"), 
                        guide = guide_colourbar(title.position = "right",barwidth = 1,barheight = 15),
                        limits=c(0, 0.055),
                        breaks=c(seq(0,0.05,by=0.01)))

gg <- gg + theme_classic(base_size=11,base_family="Helvetica")
gg <- gg + theme(plot.title=element_text(hjust=0.5))
#gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + coord_cartesian(clip="off")
gg <- gg + scale_x_continuous(breaks=c(c(seq(0,7,by=1))),labels=c("0",as.character(seq(1,7,by=1))),limits=c(0,7),expand=c(0,0))#expand=expand_scale(mult = c(0, 0),add = c(0.2, 0)))
gg <- gg + scale_y_continuous(breaks=c(seq(0,10,by=2)),limits=c(0,10),expand=c(0,0))
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(axis.text.x=element_text(vjust=1))
gg <- gg + theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
gg <- gg + theme(axis.title.x = element_markdown(margin = margin(t = 19,b=0),size=11))
gg <- gg + theme(axis.title.y = element_markdown(margin = margin(r =19),size=11))
gg <- gg + theme(legend.title = element_markdown(angle = -90),legend.title.align = 0.5)
#gg <- gg + theme(axis.title.y = element_text(margin = margin(r =10),size=11))
#gg <- gg + geom_abline(slope=1, intercept=0,linetype="dashed",expand=F)
#gg <- gg + geom_segment(x = -0.2, y = -0.2, xend = 8, yend = 8, linetype="dashed",colour="black")
#gg <- gg + annotate(geom="text", x=-0.75, y=7.5, label="(a)",size=5)
gg <- gg + labs(tag = '(c)')
gg <- gg + theme(plot.tag.position = c(-0.01, 0.99),plot.tag = element_text(hjust = 0))
#gg <- gg + theme(legend.position = "none")
gg3<-gg
gg3
#dev.off()

#dev.off()
#what happens at CE?
setwd(dir_out)
res=read.csv2("summary_VAeff_allW2_CE.csv",sep=";",dec=".",h=T) 
c1=read.csv2(paste("summary_coeff.csv",sep=""),sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW2")
int=c1a$moy[c1a$var=="(Intercept)"]
dsel=subset(mydatatogo,Param=="O")
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
v2=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
#v2=seq(round(min(dsel$GD),3),round(max(dsel$GD),3),by=0.0001)

for(i in 1:nrow(res)){
    #pred1=exp(int+res$moy[i]*((v2-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs)))
    pred1=exp(int+res$moy[i]*((v2-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))+c1a$moy[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$moy[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)) #adding effects of covariables did not change the relationships, but it changes the predicted value of range shift. Furthermore, the predictive value are representative of the mean number of temporal units and the mean spatial extent found in the database. The predictive values are also reprsentative of a combination of method (that define the intercept values of the model): Presence data and balanced sampling. Of course, the predictions don not consider the random effect  
    pred1=data.frame(GD=res$GD[i],pred1,vel_abs=v2)
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}
#png("/media/rom1/stock/PhyloRangeShift/CESAB_proj/meeting/meet2_1122/traitPerspective/fig1/v042023/2_LE_DDmax_velo_curve.png",unit="cm",width=15,height=15*0.66,res=300)#,width=547,height=360
#resX=subset(res,pv.inf0<0.05 | pv.sup0<0.05)
res$ID=as.character(round(res$GD,3))
pred2$ID=as.character(round(pred2$GD,3))
resX=subset(res,pv.inf0<0.05)
if(nrow(resX>0)){
    if(nrow(resX)<5){
        selN=round(seq(min(resX$GD),max(resX$GD),le=nrow(resX)),3)
    }else{
        selN=round(seq(min(resX$GD),max(resX$GD),le=5),3)
    }
    predX=merge(pred2,data.frame(ID=as.character(selN),signifN=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX=pred2
    predX$signifN=NA
}

resX=subset(res,pv.sup0<0.05)
if(nrow(resX>0)){
    if(nrow(resX)<5){
        selP=round(seq(min(resX$GD),max(resX$GD),le=nrow(resX)),3)
    }else{
        selP=round(seq(min(resX$GD),max(resX$GD),le=5),3)
    }
    predXCE=merge(predX,data.frame(ID=as.character(selP),signifP=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predXCE=predX
    predXCE$signifP=NA
}

predXCE1=subset(predXCE,signifN==1 | signifP==1)
gg= ggplot(predXCE, aes(vel_abs, pred1, color = vel_abs, group=ID)) +
    geom_line(data=predXCE1, aes(vel_abs, pred1, color = GD, group=ID)) +
    theme_minimal() +
    #labs(x=bquote("Nucleotide diversity"), y=bquote("Absolute velocity of range shifts (km.yr<sup>-1</sup>)"))+
    labs(x=bquote(""), y=bquote(""))+
    scale_color_viridis(name = bquote("Genetic diversity"), 
                        guide = guide_colourbar(title.position = "right",barwidth = 1,barheight = 15),
                        limits=c(0, 0.055),
                        breaks=c(seq(0,0.05,by=0.01)))

gg <- gg + theme_classic(base_size=11,base_family="Helvetica")
gg <- gg + theme(plot.title=element_text(hjust=0.5))
#gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + coord_cartesian(clip="off")
gg <- gg + scale_x_continuous(breaks=c(c(seq(0,7,by=1))),labels=c("0",as.character(seq(1,7,by=1))),limits=c(0,7),expand=c(0,0))#expand=expand_scale(mult = c(0, 0),add = c(0.2, 0)))
gg <- gg + scale_y_continuous(breaks=c(seq(0,10,by=2)),limits=c(0,10),expand=c(0,0))
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(axis.text.x=element_text(vjust=1))
gg <- gg + theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
gg <- gg + theme(axis.title.x = element_markdown(margin = margin(t = 19,b=0),size=11))
gg <- gg + theme(axis.title.y = element_markdown(margin = margin(r =19),size=11))
gg <- gg + theme(legend.title = element_markdown(angle = -90),legend.title.align = 0.5)
#gg <- gg + theme(axis.title.y = element_text(margin = margin(r =10),size=11))
#gg <- gg + geom_abline(slope=1, intercept=0,linetype="dashed",expand=F)
#gg <- gg + geom_segment(x = -0.2, y = -0.2, xend = 8, yend = 8, linetype="dashed",colour="black")
#gg <- gg + annotate(geom="text", x=-0.75, y=7.5, label="(a)",size=5)
gg <- gg + labs(tag = '(b)')
gg <- gg + theme(plot.tag.position = c(-0.01, 0.99),plot.tag = element_text(hjust = 0))
#gg <- gg + theme(legend.position = "none")
gg2<-gg
gg2
#dev.off()



png(paste0(dir_out,"/GDeffectOnVIC_v0.png"),unit="cm",width=27,height=11,res=300)#,width=547,height=360
#gs <- lapply(1:4, function(ii) 
#grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
lay <- rbind(c(1,2,3))
grid.arrange(gg1,gg2,gg3, 
             ncol=3, nrow=1, widths=c(1,1,1), heights=c(1),layout_matrix = lay)
dev.off()