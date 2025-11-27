gc();rm(list=ls())

################################################################################
#required packages
list.of.packages <- c(
    "doParallel", "parallel","foreach","pdftools","plotly",
    "pbapply","dplyr", "tidyr", "data.table",
    "scales","effects","psych", 
    "glmmTMB", "lme4", "lmerTest","here","rlist",
    "ggtext","gridExtra","grid","lattice","viridis","performance","MuMIn") 


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

################################################################################
#define the data repository
if(!dir.exists(here("Output/full_model"))){
    dir.create(here("Output/full_model"),recursive = TRUE)
}

dir.in=here("Data") #to change accordingly to the location of the data
dir.out=here("Output/full_model") #to change. It's the repository where the results are saved

# Load data
mydataset <- read.csv2(here(dir.in,"gen_data_final_fonseca2.csv"),
                       sep=",",dec=".",h=T) #file path in GitHub: /adaptive

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
        ECO, Uncertainty_Parameter, Uncertainty_Distribution, Grain_size, Data, Article_ID
    ) 

# transform continuous variables
cont_vars <- c(1:14,18, 20:26)

mydatatogo[,cont_vars] <- lapply(mydatatogo[,cont_vars], as.numeric)
mydatatogo[,-cont_vars] <- lapply(mydatatogo[,-cont_vars], function(x) factor(x, levels = unique(x)))

##eliminating species with uncertain obs
mydatatogo=subset(mydatatogo,spp!="Agrilus planipennis")
mydatatogo=subset(mydatatogo,spp!="Chrysodeixis eriosoma")
mydatatogo=droplevels(mydatatogo)

## Filter Classes with at least 5 species per Param
n_sps = 10

test <- mydatatogo %>%
    group_by(Class,Param) %>%
    summarise(N_spp = length(unique(spp))) %>% # how many species per parameter?
    dplyr::filter(N_spp >= n_sps) # select classes with > n_sps per param

test <- mutate(test, Class_Param = paste(Class, Param))

mydatatogo <- mydatatogo %>%
    mutate(Class_Param = paste(Class, Param)) %>%
    dplyr::filter(Class_Param %in% test$Class_Param) %>%
    dplyr::select(-Class_Param)

mydatatogo[,-cont_vars] <- lapply(mydatatogo[,-cont_vars], function(x) factor(x, levels = unique(x)))

## Extra fixes
# Set the reference param level to the centroid of species obs
mydatatogo$Param <- relevel(mydatatogo$Param, ref = "O") 

# Taxonomy data
taxatable <- mydatatogo %>%
    group_by(Class, Order, Family) %>%
    summarise(Species = length(unique(spp)),
              Shift = n())

write.csv(taxatable,
          here("Output/taxa_table.csv"),
          row.names = FALSE)

################################################################################
##############################Test if GD is adaptive############################
################################################################################

# plug in codon position
Fonseca <- fread(here::here("Data/Fonseca_etal_2023_EvoLetters.txt"))

all.equal(Fonseca$mut_1_bp, Fonseca$mut_2_bp)
all.equal(Fonseca$mut_1_bp, Fonseca$mut_3_bp)
summary(Fonseca[, .(mut_1_bp, mut_2_bp, mut_3_bp)])

Fonseca <- Fonseca %>% 
    select(spp_new,
           mut_1_bp, mut_2_bp, mut_3_bp)

mydatatogo <- merge(mydatatogo,
                    Fonseca,
                    by.x = "spp",
                    by.y = "spp_new")

mydatatogo$piN  <- (mydatatogo$mut_1_bp + mydatatogo$mut_2_bp) / 2
mydatatogo$piS <- mydatatogo$mut_3_bp
mydatatogo$piN_piS <- mydatatogo$piN / mydatatogo$piS

hist(mydatatogo$piN_piS)

all.equal(mydatatogo$mut_1_bp, mydatatogo$mut_2_bp)

cor.test(mydatatogo$GD, mydatatogo$mut_1_bp)
cor.test(mydatatogo$GD, mydatatogo$mut_2_bp)
cor.test(mydatatogo$GD, mydatatogo$mut_3_bp)

table(mydatatogo$Class)
################################################################################
###############Fit a single model with Param*GD*climate velocity################
################################################################################
#Obs are weighted to correct for imbalance in obs position (TE,CE,LE) and species frequency
#A GLMM model are fitted to control for method heterogeneity that assessed species range shifts and phylogeny proximity
#The model are bootstrapped to compute ccurate mean, CI95% and significance of effect size.
gam_velxGDxedge1 <- as.formula(
    "SHIFT_abs ~ scale(vel_abs) * scale(GD) * Param + 
    LogNtempUnits + LogExtent + ContinuousGrain + PrAb + Quality + 
    (Param|Class)")
s1=69 #the initial random seed
nB=12000 #the number of bootstraps
nCPU=60

#setting a local cluster
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)
clusterExport(cl = my.cluster2, varlist = c("s1","nB","glmmTMB","mydatatogo","gam_velxGDxedge1","ranef","check_singularity"))

#register cluster
doParallel::registerDoParallel(cl = my.cluster2)

ex=pblapply(1:nB, function(i) {
    
    gc(reset=T)
    print(i)
    set.seed(s1+i)
    
    n1=sample(1:nrow(mydatatogo),rep=T,size=nrow(mydatatogo))
    mydatatogo2=mydatatogo[n1,]
    mydatatogoCE=subset(mydatatogo2,Param=="O")
    x1=data.frame(table(mydatatogoCE$spp))
    x1=subset(x1,Freq>0)
    x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
    mydatatogoCE=merge(mydatatogoCE,x1[,c(1,3)],by.x="spp",by.y="Var1")
    
    mydatatogoLE=subset(mydatatogo2,Param=="LE")
    x1=data.frame(table(mydatatogoLE$spp))
    x1=subset(x1,Freq>0)
    x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
    mydatatogoLE=merge(mydatatogoLE,x1[,c(1,3)],by.x="spp",by.y="Var1")
    
    mydatatogoTE=subset(mydatatogo2,Param=="TE")
    x1=data.frame(table(mydatatogoTE$spp))
    x1=subset(x1,Freq>0)
    x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
    mydatatogoTE=merge(mydatatogoTE,x1[,c(1,3)],by.x="spp",by.y="Var1")
    
    mydatatogoTE$weight_obs_spp=mydatatogoTE$weight_obs_spp*(1/3)
    mydatatogoCE$weight_obs_spp=mydatatogoCE$weight_obs_spp*(1/3)
    mydatatogoLE$weight_obs_spp=mydatatogoLE$weight_obs_spp*(1/3)
    
    mydatatogo2=rbind(mydatatogoTE,mydatatogoCE,mydatatogoLE)
    
    gamX1 <- glmmTMB(gam_velxGDxedge1, 
                     family = Gamma(link = "log"),
                     weights = weight_obs_spp,
                     REML=F,
                     data = mydatatogo2)
    
    r2=round(MuMIn::r.squaredGLMM(gamX1),2)
    r2=data.frame(r2m=r2[1,1],
                  r2c=r2[1,2],
                  moy_GD=mean(mydatatogo2$GD), 
                  sd_GD=sd(mydatatogo2$GD),
                  moy_VC=mean(mydatatogo2$vel_abs), 
                  sd_VC=sd(mydatatogo2$vel_abs),
                  singu=check_singularity(gamX1)[[1]],
                  nB=i)
    
    coeff=data.frame(summary(gamX1)$coeff$cond,
                     var=row.names(summary(gamX1)$coeff$cond),
                     nB=i)
    
    randEff=as.data.frame(ranef(gamX1,condVar=T))
    randEff$nB=i
    
    return(list(r2,coeff,randEff))
    
}, cl = my.cluster2)

parallel::stopCluster(cl = my.cluster2)

nom="allW"
randEff=rlist::list.rbind(lapply(ex,"[[",3))
coeff=rlist::list.rbind(lapply(ex,"[[",2))
r2=rlist::list.rbind(lapply(ex,"[[",1))

write.table(r2,
            here(dir.out,paste0("R2_",nom,".csv")),
            sep=";",dec=".",row=F)
write.table(coeff,
            here(dir.out,paste0("coeff_",nom,".csv")),
            sep=";",dec=".",row=F)
write.table(randEff,
            here(dir.out,paste0("randEff_",nom,".csv")),
            sep=";",dec=".",row=F)

#analysis of raw bootstrap output
test1=function(x,mu=0){
    return(1-(length(x[x>mu])/length(x)))
}

#bootstrap test: H0= coefficient equal to 0; H1= coefficient lower than 0
test2=function(x,mu=0){
    return(1-(length(x[x<mu])/length(x)))
}


mod='allW'
a=1
s1=10 #the seed to make reproducible random staff
nB=10000
setwd(dir.out)

for(i in 1:length(mod)){
    print(a)
    r2=read.csv2(paste0("R2_",mod[i],".csv"),sep=";",dec=".")
    #r2$nB=1:nrow(r2)
    r2a=subset(r2,is.na(r2m)==F)
    print(paste(mod[i],': ',nrow(r2a),sep=""))
    if(nrow(r2a)>=nB){ #random selection of nB bootstrap model (criteria: model convergence)
        set.seed(s1)
        r2a=r2a[sample(1:nrow(r2a),size=nB,replace=F),]
    }
    if(nrow(r2a)<nB){
        print(paste0("sampling size is less than ",nB," bootstraps"))
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
        
        a=a+1
    }
}

write.table(r2a,
            here(dir.out,paste0("R2_",nom,"_sel10000.csv")),
            sep=";",dec=".",row=F) #10000 selected bootstrapped models

write.table(resR2_ok,
            here(dir.out,paste0("summary_R2.csv")),
            sep=";",dec=".",row=F) 

#R2 statistic computed from 5000 boostrap model:
#type= R2 type (r2m= marginal R2; r2c= conditional R2)
#moy= mean R2 value
#sd= R2 standard deviation value
#median= median R2 value
#q025= 2.5% quantile of the bootstrap distribution
#q975= 97.5% quantile of the bootstrap distribution
#n_sing= number of singular model
#model= model name abbreviation

write.table(res_ok,
            here(dir.out,"summary_coeff.csv"),
            sep=";",dec=".",row=F) 
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


# #significant positive effect
subset(res_ok,model=="allW" & pv.sup0<0.05)

# #significant negative effect
subset(res_ok,model=="allW" & pv.inf0<0.05)

# #non-significant negative effect
subset(res_ok,model=="allW" & pv.inf0>=0.05 & pv.sup0>=0.05)

### Test climate velocity value for which GD is significant
#### allW2 model
nB=r2a$nB
c1=read.csv2(here(dir.out,"coeff_allW.csv"),
             sep=";",dec=".",h=T)
dsel=mydatatogo
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)

my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

clusterExport(cl = my.cluster2, varlist = c("c1","dsel","v1","test1","test2"))

#register cluster
doParallel::registerDoParallel(cl = my.cluster2)

ex=pblapply(nB,function(i){
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
}, cl = my.cluster2)

parallel::stopCluster(cl = my.cluster2)

ex <- do.call("rbind",ex)

#summary statistics
resCE=data.frame(vel_abs=names(tapply(ex$GDeff_CE,ex$vel_abs,mean)),
                 moy=tapply(ex$GDeff_CE,ex$vel_abs,mean),
                 sd=tapply(ex$GDeff_CE,ex$vel_abs,sd),
                 median=tapply(ex$GDeff_CE,ex$vel_abs,median),
                 q025=tapply(ex$GDeff_CE,ex$vel_abs,quantile,probs=0.025),
                 p975=tapply(ex$GDeff_CE,ex$vel_abs,quantile,probs=0.975), 
                 pv.inf0=tapply(ex$GDeff_CE,ex$vel_abs,test2,mu=0),
                 pv.sup0=tapply(ex$GDeff_CE,ex$vel_abs,test1,mu=0))

resLE=data.frame(vel_abs=names(tapply(ex$GDeff_LE,ex$vel_abs,mean)),
                 moy=tapply(ex$GDeff_LE,ex$vel_abs,mean),
                 sd=tapply(ex$GDeff_LE,ex$vel_abs,sd),
                 median=tapply(ex$GDeff_LE,ex$vel_abs,median),
                 q025=tapply(ex$GDeff_LE,ex$vel_abs,quantile,probs=0.025),
                 p975=tapply(ex$GDeff_LE,ex$vel_abs,quantile,probs=0.975),
                 pv.inf0=tapply(ex$GDeff_LE,ex$vel_abs,test2,mu=0),
                 pv.sup0=tapply(ex$GDeff_LE,ex$vel_abs,test1,mu=0))

resTE=data.frame(vel_abs=names(tapply(ex$GDeff_TE,ex$vel_abs,mean)),
                 moy=tapply(ex$GDeff_TE,ex$vel_abs,mean),
                 sd=tapply(ex$GDeff_TE,ex$vel_abs,sd),
                 median=tapply(ex$GDeff_TE,ex$vel_abs,median),
                 q025=tapply(ex$GDeff_TE,ex$vel_abs,quantile,probs=0.025),
                 p975=tapply(ex$GDeff_TE,ex$vel_abs,quantile,probs=0.975),
                 pv.inf0=tapply(ex$GDeff_TE,ex$vel_abs,test2,mu=0),
                 pv.sup0=tapply(ex$GDeff_TE,ex$vel_abs,test1,mu=0))


write.table(resTE,
            here(dir.out,"summary_GDeff_allW_TE.csv"),
            sep=";",dec=".",row=F) 
write.table(resCE,
            here(dir.out,"summary_GDeff_allW_CE.csv"),
            sep=";",dec=".",row=F) 
write.table(resLE,
            here(dir.out,"summary_GDeff_allW_LE.csv"),
            sep=";",dec=".",row=F) 

# #significant positive effect
subset(resCE,pv.sup0<0.05)
# #significant negative effect
subset(resCE,pv.inf0<0.05)
# #non-significant negative effect
subset(resCE,pv.inf0>=0.05 & pv.sup0>=0.05)


# #significant positive effect
subset(resLE,pv.sup0<0.05)
# #significant negative effect
subset(resLE,pv.inf0<0.05)
# #non-significant negative effect
subset(resLE,pv.inf0>=0.05 & pv.sup0>=0.05)


# #significant positive effect
subset(resTE,pv.sup0<0.05)
# #significant negative effect
subset(resTE,pv.inf0<0.05)
# #non-significant negative effect
subset(resTE,pv.inf0>=0.05 & pv.sup0>=0.05)

### Test GD value for which the climate velocity effect is significant
c1=read.csv2(here(dir.out,"coeff_allW.csv"),
             sep=";",dec=".",h=T)
dsel=mydatatogo
v1=seq(round(min(dsel$GD),3),round(max(dsel$GD),3),by=0.001)

my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)

clusterExport(cl = my.cluster2, varlist = c("c1","dsel","v1","test1","test2"))

ex=pblapply(nB,function(i){
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
}, cl = my.cluster2)

parallel::stopCluster(cl = my.cluster2)

ex <- do.call("rbind",ex)

#summary statistictics
resCE=data.frame(GD=names(tapply(ex$VAeff_CE,ex$GD,mean)),
                 moy=tapply(ex$VAeff_CE,ex$GD,mean),
                 sd=tapply(ex$VAeff_CE,ex$GD,sd),
                 median=tapply(ex$VAeff_CE,ex$GD,median),
                 q025=tapply(ex$VAeff_CE,ex$GD,quantile,probs=0.025),
                 p975=tapply(ex$VAeff_CE,ex$GD,quantile,probs=0.975),
                 pv.inf0=tapply(ex$VAeff_CE,ex$GD,test2,mu=0),
                 pv.sup0=tapply(ex$VAeff_CE,ex$GD,test1,mu=0))

resLE=data.frame(GD=names(tapply(ex$VAeff_LE,ex$GD,mean)),
                 moy=tapply(ex$VAeff_LE,ex$GD,mean),
                 sd=tapply(ex$VAeff_LE,ex$GD,sd),
                 median=tapply(ex$VAeff_LE,ex$GD,median),
                 q025=tapply(ex$VAeff_LE,ex$GD,quantile,probs=0.025),
                 p975=tapply(ex$VAeff_LE,ex$GD,quantile,probs=0.975),
                 pv.inf0=tapply(ex$VAeff_LE,ex$GD,test2,mu=0),
                 pv.sup0=tapply(ex$VAeff_LE,ex$GD,test1,mu=0))

resTE=data.frame(GD=names(tapply(ex$VAeff_TE,ex$GD,mean)),
                 moy=tapply(ex$VAeff_TE,ex$GD,mean),
                 sd=tapply(ex$VAeff_TE,ex$GD,sd),
                 median=tapply(ex$VAeff_TE,ex$GD,median),
                 q025=tapply(ex$VAeff_TE,ex$GD,quantile,probs=0.025),
                 p975=tapply(ex$VAeff_TE,ex$GD,quantile,probs=0.975),
                 pv.inf0=tapply(ex$VAeff_TE,ex$GD,test2,mu=0),
                 pv.sup0=tapply(ex$VAeff_TE,ex$GD,test1,mu=0))

write.table(resTE,
            here(dir.out,"summary_VAeff_allW_TE.csv"),
            sep=";",dec=".",row=F) 
write.table(resLE,
            here(dir.out,"summary_VAeff_allW_LE.csv"),
            sep=";",dec=".",row=F) 
write.table(resCE,
            here(dir.out,"summary_VAeff_allW_CE.csv"),
            sep=";",dec=".",row=F) 

# #significant positive effect
subset(resCE,pv.sup0<0.05)
# #significant negative effect
subset(resCE,pv.inf0<0.05)
# #non-significant negative effect
subset(resCE,pv.inf0>=0.05 & pv.sup0>=0.05)


# #significant positive effect
subset(resLE,pv.sup0<0.05)
# #significant negative effect
subset(resLE,pv.inf0<0.05)
# #non-significant negative effect
subset(resLE,pv.inf0>=0.05 & pv.sup0>=0.05)


# #significant positive effect
subset(resTE,pv.sup0<0.05)
# #significant negative effect
subset(resTE,pv.inf0<0.05)
# #non-significant negative effect
subset(resTE,pv.inf0>=0.05 & pv.sup0>=0.05)

################################################################################
##########exploring latitidinal effect on the fitted relationship###############
################################################################################
plot(Lat~GD,data=mydatatogo)
plot(vel_abs~GD,data=mydatatogo)
plot(Lat~vel_abs,data=mydatatogo)
plot(log(SHIFT_abs)~vel_abs,data=mydatatogo)
#looking at the the effect of latitude on the model residuals
gam_velxGDxedge1 <- as.formula(
    "SHIFT_abs ~ scale(vel_abs) * scale(GD) * Param + 
    LogNtempUnits + LogExtent + ContinuousGrain + PrAb + Quality + 
    (1 | Class)")

#computing weights
mydatatogo2=mydatatogo
mydatatogoCE=subset(mydatatogo2,Param=="O")
x1=data.frame(table(mydatatogoCE$spp))
x1=subset(x1,Freq>0)
x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
mydatatogoCE=merge(mydatatogoCE,x1[,c(1,3)],by.x="spp",by.y="Var1")

mydatatogoLE=subset(mydatatogo2,Param=="LE")
x1=data.frame(table(mydatatogoLE$spp))
x1=subset(x1,Freq>0)
x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
mydatatogoLE=merge(mydatatogoLE,x1[,c(1,3)],by.x="spp",by.y="Var1")

mydatatogoTE=subset(mydatatogo2,Param=="TE")
x1=data.frame(table(mydatatogoTE$spp))
x1=subset(x1,Freq>0)
x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
mydatatogoTE=merge(mydatatogoTE,x1[,c(1,3)],by.x="spp",by.y="Var1")

mydatatogoTE$weight_obs_spp=mydatatogoTE$weight_obs_spp*(1/3)
mydatatogoCE$weight_obs_spp=mydatatogoCE$weight_obs_spp*(1/3)
mydatatogoLE$weight_obs_spp=mydatatogoLE$weight_obs_spp*(1/3)
mydatatogo2=rbind(mydatatogoTE,mydatatogoCE,mydatatogoLE)

#fitting the model
##is latitude explaining residuals of the model? YES, not strongly but yes
gamX1 <- glmmTMB(gam_velxGDxedge1, 
                 family = Gamma(link = "log"),
                 weights = weight_obs_spp,
                 REML=F,
                 data = mydatatogo2)
summary(gamX1)
MuMIn::r.squaredGLMM(gamX1)[1,] #37.9%


p1=predict(gamX1,mydatatogo2,re.form=~0)
plot(p1~log(mydatatogo2$SHIFT_abs))
resid=log(mydatatogo2$SHIFT_abs)-p1
hist(resid)
m1=lmer(resid~Lat+(1|Class),data=mydatatogo2,weights=weight_obs_spp)
MuMIn::r.squaredGLMM(m1)[1,] #0.6%
summary(m1) #significant but low effect on residuals: residuals increase with residuals. In the present case it means that error tends to decrease with latitude

par(mai=c(1,2,1,1))
plot(resid~Lat,data=mydatatogo2,
     xlab = "Latitude",
     ylab = "Residuals",
     main = "Model: Genetic diversity + Climate change velocity")
abline(a=-1.125e+00,b=6.441e-03,col=2)
text(2, -8, 
     "Coeff = 6.44e-03\nP-value = < 0.001\nR2 = 0.6",
     cex = .8, adj = c(0, 1))

resid0=resid

coeff=data.frame(summary(gamX1)$coeff$cond,
                 var=row.names(summary(gamX1)$coeff$cond))
coeff=coeff[,c(1,ncol(coeff))]
names(coeff)[1]="missLat"

##is climate velocity and genetic diversity explaining residuals of a model accounting for methods and latitude? YES
#meaning that climate velocity and genetic diversity drive a significant part of the determinism of the species range shifts, that is not confounded with the latitudinal gradient
form1 <- as.formula(
    "SHIFT_abs ~ scale(Lat) + 
    LogNtempUnits + LogExtent + ContinuousGrain + PrAb + Quality + 
    (1 | Class)")
gamX2 <- glmmTMB(form1, 
                 family = Gamma(link = "log"),
                 weights = weight_obs_spp,
                 REML=F,
                 data = mydatatogo2)
summary(gamX2) #low but positive effect of latitude 
MuMIn::r.squaredGLMM(gamX2)[1,] #22.1%

p1=predict(gamX2,mydatatogo2,re.form=~0)
plot(p1~log(mydatatogo2$SHIFT_abs))
resid=log(mydatatogo2$SHIFT_abs)-p1
hist(resid)
m1=lmer(resid~scale(vel_abs) * scale(GD) * Param +(1|Class),data=mydatatogo2,weights=weight_obs_spp,REML=F)
summary(m1) #GD, vel_abs and Param are significants
MuMIn::r.squaredGLMM(m1)[1,]  #6.5% des résidus
plot(resid~vel_abs,data=mydatatogo2)
plot(resid~GD,data=mydatatogo2)

coeff2=data.frame(resid=summary(m1)$coeff[,1],var=row.names(summary(m1)$coeff))
coeff=merge(coeff,coeff2)
plot(missLat~resid,data=coeff)
lm1=lm(missLat~resid,data=coeff)
summary(lm1)

plot(resid0~resid)
p3=predict(m1,mydatatogo2,re.form=~0)
plot(resid0~p3)
resid2=resid-p3
plot(resid0~resid2)

##is methods and latitude in association with climate velocity and genetic diversity explaining species range shift? YES
#meaning that climate velocity and genetic diversity drive a significant part of the determinism of the species range shifts, that is not confounded with the latitudinal gradient
form2 <- as.formula(
    "SHIFT_abs ~ scale(vel_abs) * scale(GD) * Param + scale(Lat)+
    LogNtempUnits + LogExtent + ContinuousGrain + PrAb + Quality + 
    (1 | Class)")
gamX3 <- glmmTMB(form2, 
                 family = Gamma(link = "log"),
                 weights = weight_obs_spp,
                 REML=F,
                 data = mydatatogo2)
summary(gamX3) #latitude has a positive effect
MuMIn::r.squaredGLMM(gamX3)[1,]  #37.9% => not better than the model without Latitude (gamX1 above)

coeff2=data.frame(summary(gamX3)$coeff$cond,
                 var=row.names(summary(gamX3)$coeff$cond))
coeff2=coeff2[,c(1,ncol(coeff2))]
names(coeff2)[1]="allA"
coeff=merge(coeff,coeff2)
plot(missLat~allA,data=coeff)
lm1=lm(missLat~allA,data=coeff)
summary(lm1) ##coefficeient with and without considering lat are the same which is good and demonstrate that the effect of GD, climate change velocity and position are not dependent on Latitude


form3 <- as.formula(
    "SHIFT_abs ~ scale(vel_abs) * scale(GD) * Param * scale(Lat)+
    LogNtempUnits + LogExtent + ContinuousGrain + PrAb + Quality + 
    (1 | Class)")
gamX4 <- glmmTMB(form3, 
                 family = Gamma(link = "log"),
                 weights = weight_obs_spp,
                 REML=F,
                 data = mydatatogo2,
                 na.action = "na.fail")
summary(gamX4) #latitude has a positive effect
MuMIn::r.squaredGLMM(gamX4)[1,] #38%  => not better than the model without Latitude

coeff2=data.frame(summary(gamX4)$coeff$cond,
                  var=row.names(summary(gamX4)$coeff$cond))
coeff2=coeff2[,c(1,ncol(coeff2))]
names(coeff2)[1]="allI"
coeff=merge(coeff,coeff2)
plot(missLat~allI,data=coeff)
lm1=lm(missLat~allI,data=coeff)
summary(lm1) #coefficeient with and without considering lat are the same which is good and demonstrate that the effect of GD, climate change velocity and position are not dependent on Latitude

AIC(gamX1,gamX2,gamX3,gamX4) 
# df      AIC
# gamX1 20 43.41449
# gamX2 10 23.57997
# gamX3 21 45.41115
# gamX4 32 67.40478
#Lower AIC observed when latitudes is not considered in the model
#what happens if we conducte a model selection based on AIC?

#AIC selection from the full model (not sure we need to keep this analysis)
nCPU=15
#setting a local cluster
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)
clusterExport(cl = my.cluster2, varlist = c("gamX4","pdredge","glmmTMB","mydatatogo2"))

mX=dredge(gamX4,cluster=my.cluster2)
parallel::stopCluster(cl = my.cluster2)
mX1=subset(mX,delta<4)

#correlation among covariates
form4 <- as.formula(
    "Lat ~ scale(vel_abs) * scale(GD) * Param + 
    (1 | Class)")
gamX5 <- glmmTMB(form4, 
                 family = Gamma(link = "log"),
                 weights = weight_obs_spp,
                 REML=F,
                 data = mydatatogo2)
summary(gamX5) #low but positive effect of latitude 
MuMIn::r.squaredGLMM(gamX5)[1,]  #9.7% of shared variance


