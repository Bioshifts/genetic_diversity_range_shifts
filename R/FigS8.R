gc();rm(list=ls())

################################################################################
#required packages
list.of.packages <- c(
    "doParallel", "parallel","foreach","pdftools","plotly","here",
    "pbapply","dplyr", "tidyr", "parallel",
    "scales","effects","psych", "glmmTMB", "lme4", "lmerTest","here","rlist",
    "ggtext","gridExtra","grid","lattice","viridis","performance","patchwork","cowplot","ggpubr","ggnewscale") 


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)
################################################################################
#define the data repository
dir.in=here("Output/full_model") #to change accordingly to the location of the data
dir.out=here("Figs") #to change. It's the repository where the results are saved

if(!dir.exists(dir.out)){
    dir.create(dir.out,recursive = TRUE)
}
# Load data
# file path in GitHub: /adaptive-potential/Data
mydataset <- read.csv2(here("Data","gen_data_final_fonseca2.csv"),
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
    dplyr::filter(Class_Param %in% test$Class_Param) %>%
    dplyr::select(-Class_Param)

mydatatogo[,-cont_vars] <- lapply(mydatatogo[,-cont_vars], function(x) factor(x, levels = unique(x)))

## Extra fixes
# Set the reference param level to the centroid of species obs
mydatatogo$Param <- relevel(mydatatogo$Param, ref = "O") 

################################################################################
#####################Fig. S8: adding obs as background (colored according to climate velocity)##########################
################################################################################
###############exploring the climate change dependence of the genetic diversity effect 
rampPallFunGD <- function(x){
    tmp <- colorRampPalette(colors = c("#00AFBB", "#E7B800", "#FC4E07"))(length(x))
    return(tmp)
}

GD_graph_params <- 
    list(geom_line(linewidth = 0.75), 
         labs(x = "Genetic diversity", 
              y = bquote('Velocity of range shift '(km.yr^-1)),
              color = bquote('Climate change velocity '(km.yr^-1))), 
         theme_classic(),  
         scale_x_continuous(n.breaks = 5, 
                            limits=c(0,0.055), expand=c(0,0),
                            labels=c("0",as.character(seq(0.01,0.05,length.out=5)))),  
         scale_y_continuous(breaks = seq(0,6,by=2), 
                            limits=c(0,7), 
                            expand=c(0,0)),  
         theme(legend.title = element_text(angle = -90), 
               legend.title.align = 0.5))

GD_graph_params1 <- 
    list(labs(x = "Genetic diversity", 
              y = bquote('Velocity of range shift '(km.yr^-1))), 
         theme_classic(),  
         theme(legend.title = element_text(angle = -90), 
               legend.title.align = 0.5))

GD_graph_params2 <- 
    list(geom_line(linewidth = 0.75),
         labs(x = "Genetic diversity", 
              y = bquote('Velocity of range shift '(km.yr^-1)),
              color = bquote('Climate change velocity '(km.yr^-1))),  
         theme_classic(),  
         theme(legend.title = element_text(angle = -90), 
               legend.title.align = 0.5))

## Computing significant relationship
## allW_TE GDeff
t1=subset(mydatatogo,Param=="TE") #n=295 obs

res=read.csv2(here(dir.in,"summary_GDeff_allW_TE.csv"),
              sep=";",dec=".",h=T) 
c1=read.csv2(here(dir.in,"summary_coeff.csv"),
             sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW")

int=c1a$median[c1a$var=="(Intercept)"]+c1a$median[c1a$var=="ParamTE"]

dsel=subset(mydatatogo,Param=="TE")
v1=c(0.1,1:5,5.6)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))
nv=res$vel_abs%in%v1
res=subset(res,nv==T)

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


sigcolor <- pred2 %>% group_by(vel_abs) %>% summarise(x = mean(sig))
pred2$ID=as.factor(pred2$vel_abs)
pred2$sig2=as.factor(pred2$sig)
pred2$sig2 <- relevel(pred2$sig2, ref = "1")

gg1 <- ggplot()+
    geom_point(data=t1,aes(x=GD,y=SHIFT_abs,color=vel_abs),alpha=6/10)+
    GD_graph_params1 +
    scale_color_gradientn(colors = rampPallFunGD(sigcolor$x),
                          guide = guide_colourbar(title.position = "right",
                                                  barwidth = 1, barheight = 8),
                          limits=c(0, 7),
                          breaks=c(seq(0,7,by=1)))+  
    scale_y_continuous(breaks = seq(0,50,by=10), 
                       limits=c(0,34), 
                       expand=c(0,0))+
    scale_x_continuous(breaks = seq(0,0.02,by=0.01), 
                       limits=c(0,0.025), expand=c(0,0),
                       labels=c("0",as.character(seq(0.01,0.02,by=0.01))))
gg1

gg1=gg1+new_scale_colour() +
    geom_line(data=pred2, aes(GD, pred1, color = vel_abs, group=ID,linetype = sig2)) +
    GD_graph_params2 +
    scale_color_gradientn(colors = rampPallFunGD(sigcolor$x),
                          guide = guide_colourbar(title.position = "right",
                                                  barwidth = 1, barheight = 8),
                          limits=c(0, 7),
                          breaks=c(seq(0,7,by=1)))  
gg1


## allW_CE GDeff
t1=subset(mydatatogo,Param=="O") #2296 obs

res=read.csv2(here(dir.in,"summary_GDeff_allW_CE.csv"),
              sep=";",dec=".",h=T) 
c1=read.csv2(here(dir.in,"summary_coeff.csv"),
             sep=";",dec=".",h=T)

c1a=subset(c1,model=="allW")
int=c1a$median[c1a$var=="(Intercept)"]
dsel=subset(mydatatogo,Param=="O")
v1=c(0.1,1:6,6.9)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))
nv=res$vel_abs%in%v1
res=subset(res,nv==T)

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

sigcolor <- pred2 %>% group_by(vel_abs) %>% summarise(x = mean(sig))
pred2$ID=as.factor(pred2$vel_abs)
pred2$sig2=as.factor(pred2$sig)
pred2$sig2 <- relevel(pred2$sig2, ref = "1")

gg2 <- ggplot()+
    geom_point(data=t1,aes(x=GD,y=SHIFT_abs,color=vel_abs),alpha=6/10)+
    GD_graph_params1 +
    scale_color_gradientn(colors = rampPallFunGD(sigcolor$x),
                          guide = guide_colourbar(title.position = "right",
                                                  barwidth = 1, barheight = 8),
                          limits=c(0, 7),
                          breaks=c(seq(0,7,by=1)))+  
    scale_y_continuous(breaks = seq(0,50,by=10), 
                       limits=c(0,44), 
                       expand=c(0,0))+
    scale_x_continuous(breaks = seq(0,0.05,by=0.01), 
                       limits=c(0,0.051), expand=c(0,0),
                       labels=c("0",as.character(seq(0.01,0.05,by=0.01))))
gg2

gg2=gg2+new_scale_colour() +
    geom_line(data=pred2, aes(GD, pred1, color = vel_abs, group=ID,linetype = sig2))+
    GD_graph_params2 +
    scale_color_gradientn(colors = rampPallFunGD(sigcolor$x),
                          guide = guide_colourbar(title.position = "right",
                                                  barwidth = 1, barheight = 8),
                          limits=c(0, 7),
                          breaks=c(seq(0,7,by=1)))  
gg2

## allW_LE GDeff
t1=subset(mydatatogo,Param=="LE") #2082 obs

res=read.csv2(here(dir.in,"summary_GDeff_allW_LE.csv"),
              sep=";",dec=".",h=T) 
c1=read.csv2(here(dir.in,"summary_coeff.csv"),
             sep=";",dec=".",h=T)

c1a=subset(c1,model=="allW")
int=c1a$median[c1a$var=="(Intercept)"]+c1a$median[c1a$var=="ParamLE"]

dsel=subset(mydatatogo,Param=="LE")
v1=c(0:5,5.8)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))
nv=res$vel_abs%in%v1
res=subset(res,nv==T)

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

sigcolor <- pred2 %>% group_by(vel_abs) %>% summarise(x = mean(sig))
pred2$ID=as.factor(pred2$vel_abs)
pred2$sig2=as.factor(pred2$sig)
pred2$sig2 <- relevel(pred2$sig2, ref = "1")
gg3 <- ggplot()+
    geom_point(data=t1,aes(x=GD,y=SHIFT_abs,color=vel_abs),alpha=6/10)+
    geom_line(data=pred2, aes(GD, pred1, color = vel_abs, group=ID,linetype = sig2))+
    GD_graph_params2 +
    scale_color_gradientn(colors = rampPallFunGD(sigcolor$x),
                          guide = guide_colourbar(title.position = "right",
                                                  barwidth = 1, barheight = 8),
                          limits=c(0, 7),
                          breaks=c(seq(0,7,by=1)))+
    scale_y_continuous(breaks = seq(0,50,by=10), 
                       limits=c(0,54), 
                       expand=c(0,0))+
    scale_x_continuous(breaks = seq(0,0.05,by=0.01), 
                       limits=c(0,0.055), expand=c(0,0),
                       labels=c("0",as.character(seq(0.01,0.05,by=0.01))))+
    guides(linetype = "none")
gg3

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
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))
gg_leg

p1a=p1
p2a=p2
p3a=p3
gg_lega=gg_leg

png(paste0(dir.out,"/FigS4.png"),unit="cm",width=27,height=11,res=300)#,width=547,height=360
(p1a + p2a + p3a + gg_lega ) + 
    plot_layout( nrow = 1, widths = c(1,1,1,0.4)) #+ # common axes => add axis_titles = "collect"
#plot_annotation(tag_levels = list(c('(a)','(b)','(c)',''))) # figure tags
dev.off()


