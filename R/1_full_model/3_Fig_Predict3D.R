################################################################################
#required packages
list.of.packages <- c(
    "doParallel", "parallel","foreach","pdftools","plotly",
    "pbapply","dplyr", "tidyr", "parallel",
    "scales","effects","psych", "glmmTMB", "lme4", "lmerTest","here","rlist","raster","RColorBrewer",
    "ggtext","gridExtra","grid","lattice","viridis","performance","patchwork","cowplot","ggpubr","ggnewscale","plot.matrix","plotly","DT","ggthemes") 


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


############################Figure v1###########################################
# just a try with volcano dataset - volcano is a numeric matrix that ships with R
fig <- plot_ly(z = ~volcano)
fig <- fig %>% add_surface()
fig  


## allW_TE
setwd(dir.out)
c1=read.csv2("summary_coeff.csv",
             sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW")

int=c1a$median[c1a$var=="(Intercept)"]+c1a$median[c1a$var=="ParamTE"]

dsel=subset(mydatatogo,Param=="TE")

dsel$GD2=round(dsel$GD,4)
dsel$vel_abs2=round(dsel$vel_abs,1)
kd <- with(dsel, MASS::kde2d(GD2,vel_abs2, n = 50))

#3D plots of the 2D kernel (we see where obs. are)
fig <- plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface()
fig

plot(kd$z)


v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))  
for(i in 1:length(kd$y)){
    pred1=exp(int+(c1a$median[c1a$var=="scale(GD)"]+c1a$median[c1a$var=="scale(GD):ParamTE"])*((kd$x-mean(mydatatogo$GD))/sd(mydatatogo$GD))+
                  (c1a$median[c1a$var=="scale(vel_abs)"]+c1a$median[c1a$var=="scale(vel_abs):ParamTE"])*((kd$y[i]-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))+
                  (c1a$median[c1a$var=="scale(vel_abs):scale(GD)"]+c1a$median[c1a$var=="scale(vel_abs):scale(GD):ParamTE"])*((kd$y[i]-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))*((kd$x-mean(mydatatogo$GD))/sd(mydatatogo$GD))+
                  c1a$median[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$median[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)+c1a$median[c1a$var=="ContinuousGrain"]*2) 
    pred1=data.frame(GD=kd$x,pred1,
                     vel_abs=kd$y[i],freq=kd$z[,i])
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}

#3D plots with enveloppe and colored with agradient of the number of obs available
# Response surface (with grid on axes)
x <- plot_ly(pred2, x = ~GD, y = ~vel_abs, z = ~pred1,
             intensity=~freq,
             type='mesh3d',
             colorscale = list(c(0,'white'),
                               c(0.1,'red'),
                               c(0.25, 'orange'),
                               c(0.5,'yellow'),
                               c(1, 'green')),
             colorbar=list(title = "N"))

x_te <- x %>% layout(scene = list(xaxis = list(title = 'Genetic diversity', showgrid = TRUE, gridcolor = "black"),
                                  yaxis = list(title = 'Absolute climate change velocity', showgrid = TRUE, gridcolor = "black"),
                                  zaxis = list(title = 'Predicted absolute species range shift', showgrid = TRUE, gridcolor = "black"),
                                  colorbar = list(title = "N")))
x_te  

kd_te=kd

#######################a 2d plot reresenting predictions of the model
#TE
r1=rasterFromXYZ(pred2[,c("GD","vel_abs","pred1")])
plot(r1,xlim=c(0,0.3),ylim=c(0,6))

r1_df <- as.data.frame(r1, xy = TRUE)
r1_df <- setNames(r1_df, c("GD", "vel_abs", "pred1")) #c("duration (°C)", "rate of climate warming (°C/yr)", "Colonisation/Extinction balance")

#ggplot() + 
#geom_raster(data = r1_df, aes(x =dur, y = ratT, fill = nB)) + 
#scale_fill_distiller(palette = "BrBG") + 
#coord_sf(crs = 4326)
gg=ggplot(data = r1_df, aes(x =GD, y = vel_abs, fill = pred1))  
gg <- gg + geom_tile(color="grey", size=0.1)
#gg<-  gg + theme_minimal()
gg <- gg + scale_fill_viridis(option="E",
                              #name = stringr::str_wrap("Predictions of range shift velocity "),
                              name = bquote("Predictions of range shift velocity "(km.yr^-1)),
                              guide = guide_colourbar(title.position = "left",barwidth = 1,barheight = 15),
                              limits=c(0, 12.1),
                              breaks=c(seq(0,12,by=2)),
                              na.value = "white",
                              direction=-1)
gg <- gg + labs(x="Genetic diversity", y=bquote("Climate change velocity (°C.yr<sup>-1</sup>)"), title="Trailing edge")
gg <- gg + theme_tufte(base_size=11,base_family="Helvetica")
gg <- gg + theme(plot.title=element_text(hjust=0.5))
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + coord_cartesian(expand = FALSE)
gg <- gg + scale_x_continuous(breaks=seq(0,0.025,by=0.005),limits=c(-0.001,0.025),expand=c(0,0),labels=c(0,seq(0.005,0.025,by=0.005)))
gg <- gg + scale_y_continuous(breaks=seq(0,5,by=1),limits=c(0,5.7),expand=c(0,0))
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(axis.text.x=element_text(vjust=1))
gg <- gg + theme(axis.title.x = element_text(margin = margin(t = 8,b=0),size=11))
gg <- gg + theme(axis.title.y = element_markdown(margin = margin(r =10),size=11))
#gg <- gg + theme(axis.title.y = element_text(margin = margin(r =10),size=11))
gg <- gg + theme(legend.title=element_text(size=11,angle=-270))
gg <- gg + theme(legend.title.align = 0.5,
                 legend.direction = "vertical")
gg <- gg + theme(legend.text=element_text(size=10))
gg1=gg
gg1


## allW2_CE
setwd(dir.out)
c1=read.csv2("summary_coeff.csv",
             sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW")

int=c1a$median[c1a$var=="(Intercept)"]

dsel=subset(mydatatogo,Param=="O")

dsel$GD2=round(dsel$GD,4)
dsel$vel_abs2=round(dsel$vel_abs,1)
kd <- with(dsel, MASS::kde2d(GD2,vel_abs2, n = 50))

#3D plots of the 2D kernel (we see where obs. are)
fig <- plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface()
fig

plot(kd$z)


v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))  
for(i in 1:length(kd$y)){
    pred1=exp(int+(c1a$median[c1a$var=="scale(GD)"])*((kd$x-mean(mydatatogo$GD))/sd(mydatatogo$GD))+
                  (c1a$median[c1a$var=="scale(vel_abs)"])*((kd$y[i]-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))+
                  (c1a$median[c1a$var=="scale(vel_abs):scale(GD)"])*((kd$y[i]-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))*((kd$x-mean(mydatatogo$GD))/sd(mydatatogo$GD))+
                  c1a$median[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$median[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)+c1a$median[c1a$var=="ContinuousGrain"]*2) 
    pred1=data.frame(GD=kd$x,pred1,
                     vel_abs=kd$y[i],freq=kd$z[,i])
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}

#3D plots with enveloppe and colored with agradient of the number of obs available
# Response surface (with grid on axes)
x <- plot_ly(pred2, x = ~GD, y = ~vel_abs, z = ~pred1,
             intensity=~freq,
             type='mesh3d',
             colorscale = list(c(0,'white'),
                               c(0.1,'red'),
                               c(0.25, 'orange'),
                               c(0.5,'yellow'),
                               c(1, 'green')),
             colorbar=list(title = "N"))

x_ce <- x %>% layout(scene = list(xaxis = list(title = 'Genetic diversity', showgrid = TRUE, gridcolor = "black"),
                                  yaxis = list(title = 'Absolute climate change velocity', showgrid = TRUE, gridcolor = "black"),
                                  zaxis = list(title = 'Predicted absolute species range shift', showgrid = TRUE, gridcolor = "black"),
                                  colorbar = list(title = "N")))
x_ce  

kd_ce=kd

#######################a 2d plot reresenting predictions of the model
#CE
r1=rasterFromXYZ(pred2[,c("GD","vel_abs","pred1")])
plot(r1,xlim=c(0,0.3),ylim=c(0,6))

r1_df <- as.data.frame(r1, xy = TRUE)
r1_df <- setNames(r1_df, c("GD", "vel_abs", "pred1")) #c("duration (°C)", "rate of climate warming (°C/yr)", "Colonisation/Extinction balance")

#ggplot() + 
#geom_raster(data = r1_df, aes(x =dur, y = ratT, fill = nB)) + 
#scale_fill_distiller(palette = "BrBG") + 
#coord_sf(crs = 4326)
gg=ggplot(data = r1_df, aes(x =GD, y = vel_abs, fill = pred1))  
gg <- gg + geom_tile(color="grey", size=0.1)
#gg<-  gg + theme_minimal()
gg <- gg + scale_fill_viridis(option="E",
                              name = bquote("Predictions of range shift velocity "(km.yr^-1)),
                              guide = guide_colourbar(title.position = "left",barwidth = 1,barheight = 15),
                              limits=c(0, 12.1),
                              breaks=c(seq(0,12,by=2)),
                              na.value = "white",
                              direction=-1)

gg <- gg + labs(x="Genetic diversity", y=bquote("Climate change velocity (°C.yr<sup>-1</sup>)"), title="Center")
gg <- gg + theme_tufte(base_size=11,base_family="Helvetica")
gg <- gg + theme(plot.title=element_text(hjust=0.5))
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + coord_cartesian(expand = FALSE)
gg <- gg + scale_x_continuous(breaks=seq(0,0.05,by=0.01),limits=c(-0.001,0.052),expand=c(0,0),labels=c(0,seq(0.01,0.05,by=0.01)))
gg <- gg + scale_y_continuous(breaks=seq(0,7,by=1),limits=c(0,7),expand=c(0,0))
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(axis.text.x=element_text(vjust=1))
gg <- gg + theme(axis.title.x = element_text(margin = margin(t = 8,b=0),size=11))
gg <- gg + theme(axis.title.y = element_markdown(margin = margin(r =10),size=11))
#gg <- gg + theme(axis.title.y = element_text(margin = margin(r =10),size=11))
gg <- gg + theme(legend.title=element_text(size=11,angle=-270))
gg <- gg + theme(legend.title.align = 0.5,
                 legend.direction = "vertical")
gg <- gg + theme(legend.text=element_text(size=10))
gg2=gg
gg2


## allW2_TE
setwd(dir.out)
c1=read.csv2("summary_coeff.csv",
             sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW")

int=c1a$median[c1a$var=="(Intercept)"]+c1a$median[c1a$var=="ParamLE"]

dsel=subset(mydatatogo,Param=="LE")

dsel$GD2=round(dsel$GD,4)
dsel$vel_abs2=round(dsel$vel_abs,1)
kd <- with(dsel, MASS::kde2d(GD2,vel_abs2, n = 50))

#3D plots of the 2D kernel (we see where obs. are)
fig <- plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface()
fig

plot(kd$z)


v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))  
for(i in 1:length(kd$y)){
    pred1=exp(int+(c1a$median[c1a$var=="scale(GD)"]+c1a$median[c1a$var=="scale(GD):ParamLE"])*((kd$x-mean(mydatatogo$GD))/sd(mydatatogo$GD))+
                  (c1a$median[c1a$var=="scale(vel_abs)"]+c1a$median[c1a$var=="scale(vel_abs):ParamLE"])*((kd$y[i]-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))+
                  (c1a$median[c1a$var=="scale(vel_abs):scale(GD)"]+c1a$median[c1a$var=="scale(vel_abs):scale(GD):ParamLE"])*((kd$y[i]-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))*((kd$x-mean(mydatatogo$GD))/sd(mydatatogo$GD))+
                  c1a$median[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$median[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)+c1a$median[c1a$var=="ContinuousGrain"]*2) 
    pred1=data.frame(GD=kd$x,pred1,
                     vel_abs=kd$y[i],freq=kd$z[,i])
    if(i==1){
        pred2=pred1
    }else{
        pred2=rbind(pred2,pred1)
    }
}

#3D plots with enveloppe and colored with agradient of the number of obs available
# Response surface (with grid on axes)
x <- plot_ly(pred2, x = ~GD, y = ~vel_abs, z = ~pred1,
             intensity=~freq,
             type='mesh3d',
             colorscale = list(c(0,'white'),
                               c(0.1,'red'),
                               c(0.25, 'orange'),
                               c(0.5,'yellow'),
                               c(1, 'green')),
             colorbar=list(title = "N"))

x_le <- x %>% layout(scene = list(xaxis = list(title = 'Genetic diversity', showgrid = TRUE, gridcolor = "black"),
                                  yaxis = list(title = 'Absolute climate change velocity', showgrid = TRUE, gridcolor = "black"),
                                  zaxis = list(title = 'Predicted absolute species range shift', showgrid = TRUE, gridcolor = "black"),
                                  colorbar = list(title = "N")))
x_le  

kd_le=kd

#######################a 2d plot reresenting predictions of the model
#LE
r1=rasterFromXYZ(pred2[,c("GD","vel_abs","pred1")])
plot(r1,xlim=c(0,0.3),ylim=c(0,6))

r1_df <- as.data.frame(r1, xy = TRUE)
r1_df <- setNames(r1_df, c("GD", "vel_abs", "pred1")) #c("duration (°C)", "rate of climate warming (°C/yr)", "Colonisation/Extinction balance")

#ggplot() + 
#geom_raster(data = r1_df, aes(x =dur, y = ratT, fill = nB)) + 
#scale_fill_distiller(palette = "BrBG") + 
#coord_sf(crs = 4326)
gg=ggplot(data = r1_df, aes(x =GD, y = vel_abs, fill = pred1))  
gg <- gg + geom_tile(color="grey", size=0.1)
#gg<-  gg + theme_minimal()
gg <- gg + scale_fill_viridis(option="E",
                              name = bquote("Predictions of range shift velocity "(km.yr^-1)),
                              guide = guide_colourbar(title.position = "left",barwidth = 1,barheight = 15),
                              limits=c(0, 12.1),
                              breaks=c(seq(0,12,by=2)),
                              na.value = "white",
                              direction=-1)

gg <- gg + labs(x="Genetic diversity", y=bquote("Climate change velocity (°C.yr<sup>-1</sup>)"), title="Leading edge")
gg <- gg + theme_tufte(base_size=11,base_family="Helvetica")
gg <- gg + theme(plot.title=element_text(hjust=0.5))
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + coord_cartesian(expand = FALSE)
gg <- gg + scale_x_continuous(breaks=seq(0,0.05,by=0.01),limits=c(-0.001,0.056),expand=c(0,0),labels=c(0,seq(0.01,0.05,by=0.01)))
gg <- gg + scale_y_continuous(breaks=seq(0,5,by=1),limits=c(-0.1,5.9),expand=c(0,0))
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(axis.text.x=element_text(vjust=1))
gg <- gg + theme(axis.title.x = element_text(margin = margin(t = 8,b=0),size=11))
gg <- gg + theme(axis.title.y = element_markdown(margin = margin(r =10),size=11))
#gg <- gg + theme(axis.title.y = element_text(margin = margin(r =10),size=11))
gg <- gg + theme(legend.title=element_text(size=11,angle=-270))
gg <- gg + theme(legend.title.align = 0.5,
                 legend.direction = "vertical")
gg <- gg + theme(legend.text=element_text(size=10))
gg3=gg
gg3

leg <- get_legend(gg)
gg_leg<-ggpubr::as_ggplot(leg) 
gg_leg<-gg_leg + theme(plot.margin = margin(0, 1.5, 0, 1, "cm"))
gg_leg

gg1 <- gg1+ theme(legend.position = "none")+
    ggtitle("Trailing edge (TE)") +labs(tag = '(a)')+
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag.position = c(0.065, 0.975))
gg2 <- gg2+ theme(legend.position = "none")+
    ggtitle("Centroïd (CE)") +labs(tag = '(b)')+
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag.position = c(0.065, 0.975))
gg3 <- gg3+ theme(legend.position = "none")+
    ggtitle("Leading edge (LE)") +labs(tag = '(c)')+
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag.position = c(0.065, 0.975))

png(paste0(dir.out,"/fig_2Dpred.png"),unit="cm",width=27,height=11,res=300)#,width=547,height=360
#gs <- lapply(1:4, function(ii) 
#grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
lay <- rbind(c(1,2,3,4))
gridExtra::grid.arrange(gg1,gg2,gg3,gg_leg, 
                        ncol=4, nrow=1, widths=c(1,1,1,0.2), heights=c(1),layout_matrix = lay)
dev.off()


###################let's try a 2D figures to show where are the obs#############
###################to combine with fig2_2Dpred.png##############################
#TE
a=1
kd=kd_te
for(i in kd_te$y){
    res=data.frame(GD=kd$x,vel_abs=i,n=kd$z[,a])
    if(a==1){
        res2=res
    }else{
        res2=rbind(res2,res)
    }
    a=a+1
}


r1=rasterFromXYZ(res2[,c("GD","vel_abs","n")])
plot(r1,xlim=c(0,0.3),ylim=c(0,6))

r1_df <- as.data.frame(r1, xy = TRUE)
r1_df <- setNames(r1_df, c("GD", "vel_abs", "n")) #c("duration (°C)", "rate of climate warming (°C/yr)", "Colonisation/Extinction balance")

r1_df$n2=r1_df$n
r1_df$n2[r1_df$n<0.1]=NA
gg=ggplot(data = r1_df, aes(x =GD, y = vel_abs, fill = n2))  
gg <- gg + geom_tile(color="grey", size=0.1)
#gg<-  gg + theme_minimal()
gg <- gg + scale_fill_viridis(option="B",
                              name = stringr::str_wrap("Number of observations"), 
                              guide = guide_colourbar(title.position = "left",barwidth = 1,barheight = 15),
                              limits=c(0.1, 210),
                              breaks=c(1,seq(50,200,by=50)),
                              na.value = "white",
                              direction=-1)

gg <- gg + labs(x="Genetic diversity", y=bquote("Climate change velocity (°C.yr<sup>-1</sup>)"), title="")
gg <- gg + theme_tufte(base_size=11,base_family="Helvetica")
gg <- gg + theme(plot.title=element_text(hjust=0.5))
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + coord_cartesian(expand = FALSE)
gg <- gg + scale_x_continuous(breaks=seq(0,0.025,by=0.005),limits=c(-0.001,0.025),expand=c(0,0),labels=c(0,seq(0.005,0.025,by=0.005)))
gg <- gg + scale_y_continuous(breaks=seq(0,5,by=1),limits=c(0,5.7),expand=c(0,0))
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(axis.text.x=element_text(vjust=1))
gg <- gg + theme(axis.title.x = element_text(margin = margin(t = 8,b=0),size=11))
gg <- gg + theme(axis.title.y = element_markdown(margin = margin(r =10),size=11))
#gg <- gg + theme(axis.title.y = element_text(margin = margin(r =10),size=11))
gg <- gg + theme(legend.title=element_text(size=11,angle=-270))
gg <- gg + theme(legend.title.align = 0.5,
                 legend.direction = "vertical")
gg <- gg + theme(legend.text=element_text(size=10))
gg1=gg
gg1


#CE
a=1
kd=kd_ce
for(i in kd$y){
    res=data.frame(GD=kd$x,vel_abs=i,n=kd$z[,a])
    if(a==1){
        res2=res
    }else{
        res2=rbind(res2,res)
    }
    a=a+1
}


r1=rasterFromXYZ(res2[,c("GD","vel_abs","n")])
plot(r1,xlim=c(0,0.3),ylim=c(0,6))

r1_df <- as.data.frame(r1, xy = TRUE)
r1_df <- setNames(r1_df, c("GD", "vel_abs", "n")) #c("duration (°C)", "rate of climate warming (°C/yr)", "Colonisation/Extinction balance")

#ggplot() + 
#geom_raster(data = r1_df, aes(x =dur, y = ratT, fill = nB)) + 
#scale_fill_distiller(palette = "BrBG") + 
#coord_sf(crs = 4326)
r1_df$n2=r1_df$n
r1_df$n2[r1_df$n<0.1]=NA
gg=ggplot(data = r1_df, aes(x =GD, y = vel_abs, fill = n2))  
gg <- gg + geom_tile(color="grey", size=0.1)
#gg<-  gg + theme_minimal()
gg <- gg + scale_fill_viridis(option="B",
                              name = stringr::str_wrap("Number of observations"), 
                              guide = guide_colourbar(title.position = "left",barwidth = 1,barheight = 15),
                              limits=c(0, 210),
                              breaks=c(1,seq(50,200,by=50)),
                              na.value = "white",
                              direction=-1)

gg <- gg + labs(x="Genetic diversity", y=bquote("Climate change velocity (°C.yr<sup>-1</sup>)"), title="")
gg <- gg + theme_tufte(base_size=11,base_family="Helvetica")
gg <- gg + theme(plot.title=element_text(hjust=0.5))
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + coord_cartesian(expand = FALSE)
gg <- gg + scale_x_continuous(breaks=seq(0,0.05,by=0.01),limits=c(-0.001,0.052),expand=c(0,0),labels=c(0,seq(0.01,0.05,by=0.01)))
gg <- gg + scale_y_continuous(breaks=seq(0,7,by=1),limits=c(0,7),expand=c(0,0))
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(axis.text.x=element_text(vjust=1))
gg <- gg + theme(axis.title.x = element_text(margin = margin(t = 8,b=0),size=11))
gg <- gg + theme(axis.title.y = element_markdown(margin = margin(r =10),size=11))
#gg <- gg + theme(axis.title.y = element_text(margin = margin(r =10),size=11))
gg <- gg + theme(legend.title=element_text(size=11,angle=-270))
gg <- gg + theme(legend.title.align = 0.5,
                 legend.direction = "vertical")
gg <- gg + theme(legend.text=element_text(size=10))
gg2=gg
gg2

#LE
a=1
kd=kd_le
for(i in kd$y){
    res=data.frame(GD=kd$x,vel_abs=i,n=kd$z[,a])
    if(a==1){
        res2=res
    }else{
        res2=rbind(res2,res)
    }
    a=a+1
}


r1=rasterFromXYZ(res2[,c("GD","vel_abs","n")])
plot(r1,xlim=c(0,0.3),ylim=c(0,6))

r1_df <- as.data.frame(r1, xy = TRUE)
r1_df <- setNames(r1_df, c("GD", "vel_abs", "n")) #c("duration (°C)", "rate of climate warming (°C/yr)", "Colonisation/Extinction balance")

#ggplot() + 
#geom_raster(data = r1_df, aes(x =dur, y = ratT, fill = nB)) + 
#scale_fill_distiller(palette = "BrBG") + 
#coord_sf(crs = 4326)
r1_df$n2=r1_df$n
r1_df$n2[r1_df$n<0.1]=NA
gg=ggplot(data = r1_df, aes(x =GD, y = vel_abs, fill = n2))  
gg <- gg + geom_tile(color="grey", size=0.1)
#gg<-  gg + theme_minimal()
gg <- gg + scale_fill_viridis(option="B",
                              name = stringr::str_wrap("Number of observations"), 
                              guide = guide_colourbar(title.position = "left",barwidth = 1,barheight = 15),
                              limits=c(0, 210),
                              breaks=c(1,seq(50,200,by=50)),
                              na.value = "white",
                              direction=-1)

gg <- gg + labs(x="Genetic diversity", y=bquote("Climate change velocity (°C.yr<sup>-1</sup>)"), title="")
gg <- gg + theme_tufte(base_size=11,base_family="Helvetica")
gg <- gg + theme(plot.title=element_text(hjust=0.5))
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + coord_cartesian(expand = FALSE)
gg <- gg + scale_x_continuous(breaks=seq(0,0.05,by=0.01),limits=c(-0.001,0.056),expand=c(0,0),labels=c(0,seq(0.01,0.05,by=0.01)))
gg <- gg + scale_y_continuous(breaks=seq(0,5,by=1),limits=c(-0.1,5.9),expand=c(0,0))
gg <- gg + theme(axis.text=element_text(size=10))
gg <- gg + theme(axis.text.x=element_text(vjust=1))
gg <- gg + theme(axis.title.x = element_text(margin = margin(t = 8,b=0),size=11))
gg <- gg + theme(axis.title.y = element_markdown(margin = margin(r =10),size=11))
#gg <- gg + theme(axis.title.y = element_text(margin = margin(r =10),size=11))
gg <- gg + theme(legend.title=element_text(size=11,angle=-270))
gg <- gg + theme(legend.title.align = 0.5,
                 legend.direction = "vertical")
gg <- gg + theme(legend.text=element_text(size=10))
gg3=gg
gg3

leg <- get_legend(gg)
gg_leg<-ggpubr::as_ggplot(leg) 
gg_leg<-gg_leg + theme(plot.margin = margin(0, 1.5, 0, 1, "cm"))
gg_leg

gg1 <- gg1+ theme(legend.position = "none")+
    labs(tag = '(c)')+
    theme(plot.tag.position = c(0.065, 0.975))
gg2 <- gg2+ theme(legend.position = "none")+
    labs(tag = '(d)')+
    theme(plot.tag.position = c(0.065, 0.975))
gg3 <- gg3+ theme(legend.position = "none")+
    labs(tag = '(e)')+
    theme(plot.tag.position = c(0.065, 0.975))

png(paste0(dir.out,"/fig_2Dnobs.png"),unit="cm",width=27,height=11,res=300)#,width=547,height=360
#gs <- lapply(1:4, function(ii) 
#grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
lay <- rbind(c(1,2,3,4))
gridExtra::grid.arrange(gg1,gg2,gg3,gg_leg, 
                        ncol=4, nrow=1, widths=c(1,1,1,0.2), heights=c(1),layout_matrix = lay)
dev.off()
