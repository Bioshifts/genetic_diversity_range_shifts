################################################################################
#required packages
list.of.packages <- c(
    "doParallel", "parallel","foreach","pdftools","plotly",
    "pbapply","dplyr", "tidyr", "parallel",
    "scales","effects","psych", "glmmTMB", "lme4", "lmerTest","here","rlist","ggtext","gridExtra","grid","lattice","viridis","performance","patchwork","cowplot","ggpubr") 


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)
################################################################################
#define the data repository
dir.in="/home/rbertrand/W/Bioshift/GD_study/data_Brunno_v02022024/" #to change accordingly to the location of the data
dir.out="/home/rbertrand/W/Bioshift/GD_study/boot_analysis/REMLF2" #to change. It's the repository where the results are saved

# Load data
setwd(dir.in)
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
cont_vars <- c(1:14,18, 20:26, 44:47)

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

#Figure v1
###############################################################################
###############exploring the climate change dependence of the genetic diversity effect 
rampPallFunGD <- function(x){
    tmp <- colorRampPalette(colors = c("#00AFBB", "#E7B800", "#FC4E07"))(length(x))
    return(tmp)
}

GD_graph_params <- 
    list(geom_line(linewidth = 1), 
         labs(x = "Genetic diversity", 
              y = bquote('Velocity of range shift '(km.yr^1)),
              color = bquote('Climate change velocity '(km.yr^1))), 
         theme_classic(),  
         scale_x_continuous(n.breaks = 5, 
                            limits=c(0,0.055), expand=c(0,0),
                            labels=c("0",as.character(seq(0.01,0.05,length.out=5)))),  
         scale_y_continuous(breaks = seq(0,6,by=2), 
                            limits=c(0,7), 
                            expand=c(0,0)),  
         theme(legend.title = element_text(angle = -90), 
               legend.title.align = 0.5))

## Computing significant relationship
## allW_TE GDeff
setwd(dir.out)
res=read.csv2("summary_GDeff_allW_TE.csv",
              sep=";",dec=".",h=T) 
c1=read.csv2("summary_coeff.csv",
             sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW")

int=c1a$median[c1a$var=="(Intercept)"]+c1a$median[c1a$var=="ParamTE"]

dsel=subset(mydatatogo,Param=="TE")

v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))

for(i in 1:nrow(res)){
    pred1=exp(int+res$median[i]*((v2-mean(mydatatogo$GD))/sd(mydatatogo$GD))+c1a$median[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$median[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)+c1a$median[c1a$var=="ContinuousGrain"]*2) 
    pred1=data.frame(GD=v2,pred1,
                     vel_abs=res$vel_abs[i],
                     sig = ifelse(res$pv.inf0[i]<.05 | res$pv.sup0[i]<.05,1,0),
                     lower = res$q025[i],
                     upper = res$p975[i])
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}

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
    predX2=merge(predX,data.frame(ID=as.character(selP),signifP=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX2=predX
    predX2$signifP=NA
}

predX1=subset(predX2,signifN==1 | signifP==1)

sigcolor <- predX2 %>% group_by(vel_abs) %>% summarise(x = mean(sig))

gg1 <- ggplot(data=predX1, aes(GD, pred1, color = vel_abs, group=ID)) +
    GD_graph_params +
    scale_color_gradientn(colors = rampPallFunGD(sigcolor$x),
                          guide = guide_colourbar(title.position = "right",
                                                  barwidth = 1, barheight = 8),
                          limits=c(0, 7),
                          breaks=c(seq(0,7,by=1)))  


## allW_CE GDeff
res=read.csv2("summary_GDeff_allW_CE.csv",
              sep=";",dec=".",h=T) 
c1=read.csv2("summary_coeff.csv",
             sep=";",dec=".",h=T)

c1a=subset(c1,model=="allW")
int=c1a$median[c1a$var=="(Intercept)"]
dsel=subset(mydatatogo,Param=="O")
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))

for(i in 1:nrow(res)){
    pred1=exp(int+res$median[i]*((v2-mean(mydatatogo$GD))/sd(mydatatogo$GD))+c1a$median[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$median[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)+c1a$median[c1a$var=="ContinuousGrain"]*2) 
    pred1=data.frame(GD=v2,pred1,
                     vel_abs=res$vel_abs[i],
                     sig = ifelse(res$pv.inf0[i]<.05 | res$pv.sup0[i]<.05,1,0),
                     lower = res$q025[i],
                     upper = res$p975[i])
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}


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
    predX2=merge(predX,data.frame(ID=as.character(selP),signifP=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX2=predX
    predX2$signifP=NA
}

predX1=subset(predX2,signifN==1 | signifP==1)

sigcolor <- predX2 %>% group_by(vel_abs) %>% summarise(x = mean(sig))

gg2 <- ggplot(data=predX1, aes(GD, pred1, color = vel_abs, group=ID)) +
    GD_graph_params +
    scale_color_gradientn(colors = rampPallFunGD(sigcolor$x),
                          guide = guide_colourbar(title.position = "right",
                                                  barwidth = 1, barheight = 8),
                          limits=c(0, 7),
                          breaks=c(seq(0,7,by=1)))  



## allW_LE GDeff
setwd(dir.out)
res=read.csv2("summary_GDeff_allW_LE.csv",
              sep=";",dec=".",h=T) 
c1=read.csv2("summary_coeff.csv",
             sep=";",dec=".",h=T)

c1a=subset(c1,model=="allW")
int=c1a$median[c1a$var=="(Intercept)"]+c1a$median[c1a$var=="ParamLE"]

dsel=subset(mydatatogo,Param=="LE")
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))

for(i in 1:nrow(res)){
    pred1=exp(int+res$median[i]*((v2-mean(mydatatogo$GD))/sd(mydatatogo$GD))+c1a$median[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$median[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)+c1a$median[c1a$var=="ContinuousGrain"]*2) 
    pred1=data.frame(GD=v2,pred1,
                     vel_abs=res$vel_abs[i],
                     sig = ifelse(res$pv.inf0[i]<.05 | res$pv.sup0[i]<.05,1,0),
                     lower = res$q025[i],
                     upper = res$p975[i])
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}

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
    predX2=merge(predX,data.frame(ID=as.character(selP),signifP=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX2=predX
    predX2$signifP=NA
}

predX1=subset(predX2,signifN==1 | signifP==1)

sigcolor <- predX2 %>% group_by(vel_abs) %>% summarise(x = mean(sig))

gg3 <- ggplot(data=predX1, aes(GD, pred1, color = vel_abs, group=ID)) +
    GD_graph_params +
    scale_color_gradientn(colors = rampPallFunGD(sigcolor$x),
                          guide = guide_colourbar(title.position = "right",
                                                  barwidth = 1, barheight = 8),
                          limits=c(0, 7),
                          breaks=c(seq(0,7,by=1)))  +
    guides(colour = guide_colourbar(title.position = "left"))+
    theme(legend.key.height = unit(1.6, "cm"),
          legend.title = element_text(size = 12, angle = 90),
          legend.title.align = 0.5,
          legend.direction = "vertical")


#chart formatting
p1 <- gg1 + 
    theme(legend.position = "none")+ 
    #labs(tag = '(a)') +
    #theme(plot.tag.position = c(0.05, 1))
    ggtitle("Trailing edge (TE)") +labs(tag = '(a)')+
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag.position = c(0.065, 0.975))

p2 <- gg2 + 
    theme(legend.position = "none") +
    #labs(tag = '(b)')+
    #theme(plot.tag.position = c(-0.05, 1))
    ggtitle("Centroid (CE)") +labs(tag = '(b)')+
    theme(plot.title = element_text(hjust = 0.5))

p3 <- gg3 + 
    theme(legend.position = "none") + 
    #labs(tag = '(c)')+
    #theme(plot.tag.position = c(-0.05, 1))
    ggtitle("Leading edge (LE)") + labs(tag = '(c)')+
    theme(plot.title = element_text(hjust = 0.5))

legend <- cowplot::get_legend(gg3)
gg_leg<-ggpubr::as_ggplot(legend) 
gg_leg<-gg_leg + 
    theme(plot.margin = margin(0, 3, 0, 0, "cm"))
gg_leg


png(paste0(dir.out,"/fig3_VICeffectOnGD_v1.png"),unit="cm",width=27,height=11,res=300)#,width=547,height=360
(p1 + p2 + p3 + gg_leg ) + 
    plot_layout( nrow = 1, widths = c(1,1,1,0.4)) #+ # common axes => add axis_titles = "collect"
#plot_annotation(tag_levels = list(c('(a)','(b)','(c)',''))) # figure tags
dev.off()

###############Exploring the modulating effect of gentic diversity
GD_graph_params <- 
    list(geom_line(linewidth = 1), 
         labs(x = bquote('Climate change velocity '(km.yr^1)), 
              y = bquote('Velocity of range shift '(km.yr^1)),
              color = "Genetic diversity"), 
         theme_classic(),  
         scale_x_continuous(breaks = seq(0,7,by=1), 
                            limits=c(0,7), 
                            expand=c(0,0)),  
         scale_y_continuous(breaks = seq(0,7,by=1), 
                            limits=c(0,7), 
                            expand=c(0,0)),  
         theme(legend.title = element_text(angle = -90), 
               legend.title.align = 0.5))

## allW
## allW_TE GDeff
setwd(dir.out)
res=read.csv2("summary_VAeff_allW_TE.csv",
              sep=";",dec=".",h=T) 
c1=read.csv2("summary_coeff.csv",
             sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW")

int=c1a$median[c1a$var=="(Intercept)"]+c1a$median[c1a$var=="ParamTE"]

dsel=subset(mydatatogo,Param=="TE")

v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))

for(i in 1:nrow(res)){
    pred1=exp(int+res$median[i]*((v1-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))+c1a$median[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$median[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)+c1a$median[c1a$var=="ContinuousGrain"]*2) 
    pred1=data.frame(vel_abs=v1,pred1,
                     GD=res$GD[i],
                     sig = ifelse(res$pv.inf0[i]<.05 | res$pv.sup0[i]<.05,1,0),
                     lower = res$q025[i],
                     upper = res$p975[i])
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}

res$ID=as.character(round(res$GD,3))
pred2$ID=as.character(round(pred2$GD,3))
resX=subset(res,pv.inf0<0.05)
if(nrow(resX>0)){
    selN=round(seq(min(resX$GD),max(resX$GD),le=5),3)
    predX=merge(pred2,data.frame(ID=as.character(selN),signifN=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX=pred2
    predX$signifN=NA
}

resX=subset(res,pv.sup0<0.05)
if(nrow(resX>0)){
    selP=round(seq(min(resX$GD),max(resX$GD),le=5),3)
    predX2=merge(predX,data.frame(ID=as.character(selP),signifP=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX2=predX
    predX2$signifP=NA
}

predX1=subset(predX2,signifN==1 | signifP==1)

sigcolor <- predX2 %>% group_by(GD) %>% summarise(x = mean(sig))

gg1 <- ggplot(data=predX1, aes(vel_abs, pred1, color = GD, group=ID)) +
    GD_graph_params +
    scale_color_viridis(name = bquote("Genetic diversity"), 
                        guide = guide_colourbar(title.position = "right",barwidth = 1,barheight = 15),
                        limits=c(0, 0.055),
                        breaks=c(seq(0,0.05,by=0.01)))


## allW_CE GDeff

res=read.csv2("summary_VAeff_allW_CE.csv",
              sep=";",dec=".",h=T) 
c1=read.csv2("summary_coeff.csv",
             sep=";",dec=".",h=T)

c1a=subset(c1,model=="allW")
int=c1a$median[c1a$var=="(Intercept)"]
dsel=subset(mydatatogo,Param=="O")
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))

for(i in 1:nrow(res)){
    pred1=exp(int+res$median[i]*((v1-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))+c1a$median[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$median[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)+c1a$median[c1a$var=="ContinuousGrain"]*2) 
    pred1=data.frame(vel_abs=v1,pred1,
                     GD=res$GD[i],
                     sig = ifelse(res$pv.inf0[i]<.05 | res$pv.sup0[i]<.05,1,0),
                     lower = res$q025[i],
                     upper = res$p975[i])
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}

res$ID=as.character(round(res$GD,3))
pred2$ID=as.character(round(pred2$GD,3))
resX=subset(res,pv.inf0<0.05)
if(nrow(resX>0)){
    selN=round(seq(min(resX$GD),max(resX$GD),le=5),3)
    predX=merge(pred2,data.frame(ID=as.character(selN),signifN=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX=pred2
    predX$signifN=NA
}

resX=subset(res,pv.sup0<0.05)
if(nrow(resX>0)){
    selP=round(seq(min(resX$GD),max(resX$GD),le=5),3)
    predX2=merge(predX,data.frame(ID=as.character(selP),signifP=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX2=predX
    predX2$signifP=NA
}

predX1=subset(predX2,signifN==1 | signifP==1)

sigcolor <- predX2 %>% group_by(GD) %>% summarise(x = mean(sig))

gg2 <- ggplot(data=predX1, aes(vel_abs, pred1, color = GD, group=ID)) +
    GD_graph_params +
    scale_color_viridis(name = bquote("Genetic diversity"), 
                        guide = guide_colourbar(title.position = "right",barwidth = 1,barheight = 15),
                        limits=c(0, 0.055),
                        breaks=c(seq(0,0.05,by=0.01)))


## allW2_LE GDeff
setwd(dir.out)
res=read.csv2("summary_VAeff_allW_LE.csv",
              sep=";",dec=".",h=T) 
c1=read.csv2("summary_coeff.csv",
             sep=";",dec=".",h=T)

c1a=subset(c1,model=="allW")
int=c1a$median[c1a$var=="(Intercept)"]+c1a$median[c1a$var=="ParamLE"]

dsel=subset(mydatatogo,Param=="LE")
v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))

for(i in 1:nrow(res)){
    pred1=exp(int+res$median[i]*((v1-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))+c1a$median[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$median[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)+c1a$median[c1a$var=="ContinuousGrain"]*2) 
    pred1=data.frame(vel_abs=v1,pred1,
                     GD=res$GD[i],
                     sig = ifelse(res$pv.inf0[i]<.05 | res$pv.sup0[i]<.05,1,0),
                     lower = res$q025[i],
                     upper = res$p975[i])
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}

res$ID=as.character(round(res$GD,3))
pred2$ID=as.character(round(pred2$GD,3))
resX=subset(res,pv.inf0<0.05)
if(nrow(resX>0)){
    selN=round(seq(min(resX$GD),max(resX$GD),le=5),3)
    predX=merge(pred2,data.frame(ID=as.character(selN),signifN=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX=pred2
    predX$signifN=NA
}

resX=subset(res,pv.sup0<0.05)
if(nrow(resX>0)){
    selP=round(seq(min(resX$GD),max(resX$GD),le=5),3)
    predX2=merge(predX,data.frame(ID=as.character(selP),signifP=1),by.x="ID",by.y="ID",all.x=T)
}else{
    predX2=predX
    predX2$signifP=NA
}

predX1=subset(predX2,signifN==1 | signifP==1)

sigcolor <- predX2 %>% group_by(GD) %>% summarise(x = mean(sig))

gg3 <- ggplot(data=predX1, aes(vel_abs, pred1, color = GD, group=ID)) +
    GD_graph_params +
    scale_color_viridis(name = bquote("Genetic diversity"), 
                        guide = guide_colourbar(title.position = "right",barwidth = 1,barheight = 15),
                        limits=c(0, 0.055),
                        breaks=c(seq(0,0.05,by=0.01)),
                        labels=c("0",as.character(seq(0.01,0.05,by=0.01)))) +
    guides(colour = guide_colourbar(title.position = "left"))+
    theme(legend.key.height = unit(1.6, "cm"),
          legend.title = element_text(size = 12, angle = 90),
          legend.title.align = 0.5,
          legend.direction = "vertical")


#formatting chart
p1 <- gg1 + 
    theme(legend.position = "none")+ 
    #labs(tag = '(a)') +
    #theme(plot.tag.position = c(0.05, 1))
    ggtitle("Trailing edge (TE)") +labs(tag = '(e)')+
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag.position = c(0.065, 0.975))

p2 <- gg2 + 
    theme(legend.position = "none") +
    #labs(tag = '(b)')+
    #theme(plot.tag.position = c(-0.05, 1))
    ggtitle("Centroid (CE)") +labs(tag = '(f)')+
    theme(plot.title = element_text(hjust = 0.5))

p3 <- gg3 + 
    theme(legend.position = "none") + 
    #labs(tag = '(c)')+
    #theme(plot.tag.position = c(-0.05, 1))
    ggtitle("Leading edge (LE)") + labs(tag = '(g)')+
    theme(plot.title = element_text(hjust = 0.5))

legend <- cowplot::get_legend(gg3)
gg_leg<-ggpubr::as_ggplot(legend) 
gg_leg<-gg_leg + 
    theme(plot.margin = margin(0, 3, 0, 0, "cm"))
gg_leg


png(paste0(dir.out,"/fig3_GDeffectOnVIC_v1.png"),unit="cm",width=27,height=11,res=300)#,width=547,height=360
(p1 + p2 + p3 + gg_leg ) + 
    plot_layout(nrow = 1, widths = c(1,1,1,0.4)) #+ # common axes => add axis_titles = "collect"
#plot_annotation(tag_levels = list(c('(a)','(b)','(c)','(d)'))) # figure tags
dev.off()

####Important the final fig3_v1 has been formatted on ppt
