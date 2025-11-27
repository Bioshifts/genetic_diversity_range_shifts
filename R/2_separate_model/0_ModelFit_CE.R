################################################################################
#required packages
list.of.packages <- c(
    "doParallel", "parallel","foreach","pdftools","plotly",
    "pbapply","dplyr", "tidyr", "parallel",
    "scales","effects","psych", "glmmTMB", "lme4", "lmerTest","here","rlist","ggtext","gridExtra","grid","lattice","viridis","performance","MuMIn") 


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)
################################################################################
#define the data repository
dir.in=here("Data") #to change accordingly to the location of the data
dir.out=here("Output/single_model") #to change. It's the repository where the results are saved

if(!dir.exists(dir.out)){
    dir.create(dir.out,recursive = TRUE)
}

# Load data

mydataset <- read.csv2(here(dir.in,"gen_data_final_fonseca2.csv"),
                       sep=",",dec=".",h=T) 

# Data selection
## Latitude data
mydatatogo <- mydataset  %>%
    dplyr::filter(Type == "LAT",
                  shift_vel_sign == "pospos" | shift_vel_sign == "negneg", # Select only shifts in the same direction of velocity
                  SHIFT != 0, # remove non-significant shifts
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
        ECO, Uncertainty_Parameter, Uncertainty_Distribution, Grain_size, Data
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
    dplyr::summarise(N_spp = length(unique(spp))) %>% # how many species per parameter?
    dplyr::filter(N_spp >= n_sps) # select classes with > n_sps per param

test <- mutate(test, Class_Param = paste(Class, Param))

mydatatogo <- mydatatogo %>%
    mutate(Class_Param = paste(Class, Param)) %>%
    filter(Class_Param %in% test$Class_Param) %>%
    select(-Class_Param)

mydatatogo[,-cont_vars] <- lapply(mydatatogo[,-cont_vars], function(x) factor(x, levels = unique(x)))

## Extra fixes
# Set the reference param level to the centroid of species obs
mydatatogo$Param <- relevel(mydatatogo$Param, ref = "O") 

################################################################################
###############Trailing edge model analysis#####################################
###############Fit a single model with GD*climate velocity######################
################################################################################
mydatatogo=subset(mydatatogo,Param=="O")
#Obs are weighted to correct for imbalance in species frequency
#A GLMM model are fitted to control for method heterogeneity that assessed species range shifts and phylogeny proximity
#The model are bootstrapped to compute ccurate mean, CI95% and significance of effect size.
gam_velxGDxedge1 <- as.formula(
    "SHIFT_abs ~ scale(vel_abs) * scale(GD) + 
    LogNtempUnits + LogExtent + ContinuousGrain + PrAb + Quality + 
    (1 | Class)")
s1=69 #the initial random seed
nB=12000 #the number of bootstraps
nCPU=detectCores() - 3

#setting a local cluster
my.cluster2 <- parallel::makeCluster(
    nCPU, 
    type = "PSOCK"
)
clusterExport(cl = my.cluster2, 
              varlist = c("s1","nB","glmmTMB","mydatatogo","gam_velxGDxedge1",
                          "ranef","check_singularity"))

#register cluster
doParallel::registerDoParallel(cl = my.cluster2)

ex=pblapply(1:nB, function(i) {
    gc(reset=T)
    print(i)
    set.seed(s1+i)
    
    n1=sample(1:nrow(mydatatogo),rep=T,size=nrow(mydatatogo))
    mydatatogo2=mydatatogo[n1,]
    x1=data.frame(table(mydatatogo2$spp))
    x1=subset(x1,Freq>0)
    x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
    mydatatogo2=merge(mydatatogo2,x1[,c(1,3)],by.x="spp",by.y="Var1")
    
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

nom="CE"
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


mod='CE'
a=1
s1=10 #the seed to make reproducible random staff
nB=10000


for(i in 1:length(mod)){
    print(a)
    r2=read.csv2(here(dir.out,paste0("R2_",mod[i],".csv")))
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
        resR2=data.frame(type="r2m",
                         moy=mean(as.numeric(as.numeric(r2a$r2m))),
                         sd=sd(as.numeric(r2a$r2m)),
                         median=median(as.numeric(r2a$r2m)),
                         q025=quantile(as.numeric(r2a$r2m),probs=0.025),
                         q975=quantile(as.numeric(r2a$r2m),probs=0.975))
        resR2=rbind(resR2,data.frame(type="r2c",
                                     moy=mean(as.numeric(r2a$r2c)),
                                     sd=sd(as.numeric(r2a$r2c)),
                                     median=median(as.numeric(r2a$r2c)),
                                     q025=quantile(as.numeric(r2a$r2c),probs=0.025),
                                     q975=quantile(as.numeric(r2a$r2c),probs=0.975)))
        #resR2$n_sing=nrow(subset(r2a,sing==T))
        resR2$model=mod[i]
        
        c1=read.csv2(here(dir.out,paste("coeff_",mod[i],".csv",sep="")),
                     sep=";",dec=".",h=T)
        #c1$nB=rep(1:nrow(r2),each=nrow(c1)/6000)
        c1=merge(c1,data.frame(nB=r2a$nB),by.x="nB",by.y="nB")
        rr3=c1
        
        res=data.frame(var=names(tapply(rr3$Estimate,rr3$var,mean)),
                       moy=tapply(rr3$Estimate,rr3$var,mean),
                       sd=tapply(rr3$Estimate,rr3$var,sd),
                       median=tapply(rr3$Estimate,rr3$var,median),
                       q025=tapply(rr3$Estimate,rr3$var,quantile,probs=0.025),
                       p975=tapply(rr3$Estimate,rr3$var,quantile,probs=0.975), 
                       pv.inf0=tapply(rr3$Estimate,rr3$var,test2,mu=0),
                       pv.sup0=tapply(rr3$Estimate,rr3$var,test1,mu=0))
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
            here(dir.out,paste0("summary_R2_",nom,".csv")),
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
            here(dir.out,paste0("summary_coeff_",nom,".csv")),
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
subset(res_ok,model=="CE" & pv.sup0<0.05)

# #significant negative effect
subset(res_ok,model=="CE" & pv.inf0<0.05)

# #non-significant negative effect
subset(res_ok,model=="CE" & pv.inf0>=0.05 & pv.sup0>=0.05)

### Test climate velocity value for which GD is significant
#### TE model
nB=r2a$nB
c1=read.csv2(here(dir.out,"coeff_CE.csv"),
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
write.table(resCE,
            here(dir.out,"summary_GDeff_CE.csv"),
            sep=";",dec=".",row=F)


# #significant positive effect
subset(resCE,pv.sup0<0.05)
# #significant negative effect
subset(resCE,pv.inf0<0.05)
# #non-significant negative effect
subset(resCE,pv.inf0>=0.05 & pv.sup0>=0.05)

### Test GD value for which the climate velocity effect is significant
c1=read.csv2(here(dir.out,"coeff_CE.csv"),
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

write.table(resCE,
            here(dir.out,"summary_VAeff_CE.csv"),
            sep=";",dec=".",row=F)

# #significant positive effect
subset(resCE,pv.sup0<0.05)
# #significant negative effect
subset(resCE,pv.inf0<0.05)
# #non-significant negative effect
subset(resCE,pv.inf0>=0.05 & pv.sup0>=0.05)

