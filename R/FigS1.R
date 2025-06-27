gc();rm(list=ls())

################################################################################
#required packages
list.of.packages <- c(
    "doParallel", "parallel","foreach","pdftools","plotly","pbapply",
    "dplyr", "tidyr","scales","effects","psych", "glmmTMB", "lme4", 
    "lmerTest","here","rlist","ggtext","metR","gridExtra","grid",
    "lattice","viridis","performance","MuMIn","raster","RColorBrewer",
    "rgl","patchwork","cowplot","ggpubr","ggnewscale","plot.matrix",
    "DT","ggthemes","fields") 

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

################################################################################
#define the data repository
dir.in=here("Output/full_model") #to change accordingly to the location of the data
dir.out=here("Figs") #to change. It's the repository where the results are saved

# Load data
mydataset <- read.csv2(here("Data","gen_data_final_fonseca2.csv"),
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


############################Figure S1###########################################

params <- c("ParamTE","ParamCE","ParamLE")
params2 <- c("TE","O","LE")
params3 <- c(":ParamTE","",":ParamLE")
param4 <- c("Trailing edge", "Centroid", "Leading edge")

png(paste0(dir.out,"/FigS1.png"),unit="cm",width=20,height=6.666,res=300)#,width=547,height=360

par(mfrow = c(1,3)) 

for(j in 1:length(unique(params))){
    
    c1=read.csv2(here(dir.in,"summary_coeff.csv"),
                 sep=";",dec=".",h=T)
    c1a=subset(c1,model=="allW")
    
    if(params2[j] == "O"){
        int=c1a$median[c1a$var=="(Intercept)"]
    } else {
        int=c1a$median[c1a$var=="(Intercept)"]+c1a$median[c1a$var==params[j]]
    }
    
    dsel=subset(mydatatogo,Param==params2[j])
    
    dsel$GD2=round(dsel$GD,4)
    dsel$vel_abs2=round(dsel$vel_abs,1)
    kd <- with(dsel, 
               MASS::kde2d(GD2,vel_abs2, 
                           n = 50,
                           lims=c(round(range(mydatatogo$GD),4),
                                  round(range(mydatatogo$vel_abs),1))))
    
    
    v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
    v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))  
    for(i in 1:length(kd$y)){
        pred1=exp(int+(c1a$median[c1a$var=="scale(GD)"]+c1a$median[c1a$var==paste0("scale(GD)",params3[j])])*((kd$x-mean(mydatatogo$GD))/sd(mydatatogo$GD))+
                      (c1a$median[c1a$var=="scale(vel_abs)"]+c1a$median[c1a$var==paste0("scale(vel_abs)",params3[j])])*((kd$y[i]-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))+
                      (c1a$median[c1a$var=="scale(vel_abs):scale(GD)"]+c1a$median[c1a$var==paste0("scale(vel_abs):scale(GD)",params3[j])])*((kd$y[i]-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))*((kd$x-mean(mydatatogo$GD))/sd(mydatatogo$GD))+
                      c1a$median[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$median[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)+c1a$median[c1a$var=="ContinuousGrain"]*2) 
        pred1=data.frame(GD=kd$x,pred1,
                         vel_abs=kd$y[i],freq=kd$z[,i])
        if(i==1){
            pred2=pred1
        }else{
            pred2=rbind(pred2,pred1)
        }
    }
    
    
    pred2=pred2[order(pred2$vel_abs),]
    pred2=pred2[order(pred2$GD),]
    z1=matrix(pred2$freq,ncol=length(unique(pred2$GD)),byrow=T)
    jet.colors2 <- colorRampPalette( c("red","orange","yellow","green") ) 
    pal2 <- jet.colors2(100)
    z=z1
    z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
    col.ind2 <- cut(z.facet.center, 100)
    col.ind2[z.facet.center<0.01]="#FFFFFF"
    s1=1:round(max(z.facet.center),0)
    col.ind3 <- cut(s1, 100)
    
    par(mar=c(0,2,2,4),
        mgp = c(8, 10, 0))
    
    with(pred2,
         persp(x=sort(unique(GD)),
               y=sort(unique(vel_abs)),
               z=matrix(pred1,ncol=length(unique(GD)),byrow=T),
               col=pal2[col.ind2],
               main = param4[j],
               ticktype = "detailed", phi=30, theta=-40,border="grey80",lwd=0.2,
               xlab='\n\nGenetic diversity',
               ylab='\n\nClimate change velocity',
               zlab='\n\nPredicted\nrange shift velocity',
               cex.lab=0.75,
               cex.axis=0.75))
    
    r1=rasterFromXYZ(pred2[,c("GD","vel_abs","pred1")])
    rX2=r1
    rX2[]=sample(1:100,size=ncell(r1),rep=T)
    
    plot(rX2, horizontal=F, legend.only=TRUE, col=pal2,
         smallplot=c(.85, .87, .4, .8), maxpixels=2000000,
         axis.args=list(tck=-0.5,
                        at=as.numeric(col.ind3[s1%in%c(1,seq(20,100,by=20))]),
                        labels=F),
         legend.args=list(text="Density of observations", 
                          side=4,font=2, 
                          line=1.2, 
                          cex=0.65))
    
    if(params2[j] == "TE"){
        plot(rX2, horizontal=F, legend.only=TRUE, col=pal2,
             smallplot=c(.85, .87, .4, .8),maxpixels=2000000,
             axis.args=list(tck=F,
                            lwd=0,
                            line=-0.50,
                            at=as.numeric(col.ind3[s1%in%c(1,seq(20,100,by=20))]),
                            labels=c(1,seq(20,100,by=20)),
                            cex.axis=0.5))
    }
    if(params2[j] == "O"){
        plot(rX2, horizontal=F, legend.only=TRUE,col=pal2,
             smallplot=c(.85, .87, .4, .8),maxpixels=2000000,
             axis.args=list(tck=F,
                            lwd=0,
                            line=-0.50,
                            at=as.numeric(col.ind3[s1%in%c(1,seq(20,60,by=20),max(s1))]),
                            labels=c(1,seq(20,60,by=20),max(s1)),
                            cex.axis=0.5))
    }
    if(params2[j] == "LE"){
        plot(rX2, horizontal=F, legend.only=TRUE,col=pal2,
             smallplot=c(.85, .87, .4, .8),maxpixels=2000000,
             axis.args=list(tck=F,
                            lwd=0,
                            line=-0.50,
                            at=as.numeric(col.ind3[s1%in%c(1,seq(25,175,by=25),max(s1))]),
                            labels=c(1,seq(25,175,by=25),max(s1)),
                            cex.axis=0.5))
    }
    
}



dev.off()


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
r1_df=subset(r1_df,vel_abs<=5.7 & GD<=0.025)
gg=ggplot(data = r1_df, aes(x =GD, y = vel_abs, fill = pred1))  
gg <- gg + geom_tile(color="grey", size=0.1)
#gg<-  gg + theme_minimal()
gg <- gg + scale_fill_viridis(option="E",
                              #name = stringr::str_wrap("Predictions of range shift velocity "),
                              name = bquote("Predictions of range shift velocity "(km.yr^-1)),
                              guide = guide_colourbar(title.position = "left",barwidth = 1,barheight = 15),
                              limits=c(0, 9.033),
                              breaks=c(seq(0,9,by=1)),
                              na.value = "white",
                              direction=-1)+
    #geom_contour(aes(x =GD, y = vel_abs, z = pred1), breaks = seq(0,12,by=2), colour = "grey85", linewidth = 0.2)+
    geom_contour2(aes(z = pred1), breaks = seq(0,12,by=1),colour = "white", linewidth = 0.4)+
    
    metR::geom_text_contour(aes(z = pred1),
                            #colour = "black", size = 4.5, fontface = "bold",
                            breaks = seq(0,12,by=1),
                            stroke = 0.15, stroke.colour = "white", # 'stroke' controls the width of stroke relative to the size of the text
                            skip = 0, # number of contours to skip
                            rotate = FALSE, # horizontal labelling; if TRUE, rotate text following the contour
                            #label.placer = label_placer_fraction(frac = 0.5)
    )
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
c1=read.csv2(here(dir.in,"summary_coeff.csv"),
             sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW")

int=c1a$median[c1a$var=="(Intercept)"]

dsel=subset(mydatatogo,Param=="O")

dsel$GD2=round(dsel$GD,4)
dsel$vel_abs2=round(dsel$vel_abs,1)
#kd <- with(dsel, MASS::kde2d(GD2,vel_abs2, n = 50))
kd <- with(dsel, MASS::kde2d(GD2,vel_abs2, n = 50,lims=c(round(range(mydatatogo$GD),4),round(range(mydatatogo$vel_abs),1))))


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


pred2=pred2[order(pred2$vel_abs),]
pred2=pred2[order(pred2$GD),]
z1=matrix(pred2$freq,ncol=length(unique(pred2$GD)),byrow=T)
jet.colors2 <- colorRampPalette( c("red","orange","yellow","green") ) 
pal2 <- jet.colors2(100)
#col.ind2 <- col.ind
#col.ind2 <- pal2[col.ind]
#col.ind2[z1<0.5]="#FFFFFF"
z=z1
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
col.ind2 <- cut(z.facet.center, 100)
col.ind2[z.facet.center<0.01]="#FFFFFF"
s1=1:round(max(z.facet.center),0)
col.ind3 <- cut(s1, 100)

png(paste0(dir.out,"/figS1_CE_v062025.png"),unit="cm",width=11,height=11,res=300)#,width=547,height=360

par(mar=c(0,2,0,4))
with(pred2,persp(x=sort(unique(GD)),y=sort(unique(vel_abs)),z=matrix(pred1,ncol=length(unique(GD)),byrow=T),
                 col=pal2[col.ind2],ticktype = "detailed", phi=30, theta=-40,border="grey80",lwd=0.2,
                 xlab="",
                 ylab='',
                 zlab="",
                 cex.lab=0.75,
                 cex.axis=0.75))

r1=rasterFromXYZ(pred2[,c("GD","vel_abs","pred1")])
rX2=r1
rX2[]=sample(1:100,size=ncell(r1),rep=T)

plot(rX2, horizontal=F, legend.only=TRUE,col=pal2,smallplot=c(.85, .87, .4, .8),maxpixels=2000000
     ,axis.args=list(tck=-0.5,at=as.numeric(col.ind3[s1%in%c(1,seq(20,60,by=20),max(s1))]),labels=F)
     ,legend.args=list(text="Density of observations", side=4,font=2, line=1.2, cex=0.65))

plot(rX2, horizontal=F, legend.only=TRUE,col=pal2,smallplot=c(.85, .87, .4, .8),maxpixels=2000000
     ,axis.args=list(tck=F,lwd=0,line=-0.50,at=as.numeric(col.ind3[s1%in%c(1,seq(20,60,by=20),max(s1))]),labels=c(1,seq(20,60,by=20),max(s1)),cex.axis=0.5))

dev.off()


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
r1_df=subset(r1_df,vel_abs<=7 & GD<=0.052)
gg=ggplot(data = r1_df, aes(x =GD, y = vel_abs, fill = pred1))  
gg <- gg + geom_tile(color="grey", size=0.1)
#gg<-  gg + theme_minimal()
gg <- gg + scale_fill_viridis(option="E",
                              #name = stringr::str_wrap("Predictions of range shift velocity "),
                              name = bquote("Predictions of range shift velocity "(km.yr^-1)),
                              guide = guide_colourbar(title.position = "left",barwidth = 1,barheight = 15),
                              limits=c(0, 9.033),
                              breaks=c(seq(0,9,by=1)),
                              na.value = "white",
                              direction=-1)+
    #geom_contour(aes(x =GD, y = vel_abs, z = pred1), breaks = seq(0,12,by=2), colour = "grey85", linewidth = 0.2)+
    geom_contour2(aes(z = pred1), breaks = seq(0,12,by=1),colour = "white", linewidth = 0.4)+
    
    metR::geom_text_contour(aes(z = pred1),
                            #colour = "black", size = 4.5, fontface = "bold",
                            breaks = seq(0,12,by=1),
                            stroke = 0.15, stroke.colour = "white", # 'stroke' controls the width of stroke relative to the size of the text
                            skip = 0, # number of contours to skip
                            rotate = FALSE, # horizontal labelling; if TRUE, rotate text following the contour
                            #label.placer = label_placer_fraction(frac = 0.5)
    )

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


## allW2_LE
c1=read.csv2(here(dir.in,"summary_coeff.csv"),
             sep=";",dec=".",h=T)
c1a=subset(c1,model=="allW")

int=c1a$median[c1a$var=="(Intercept)"]+c1a$median[c1a$var=="ParamLE"]

dsel=subset(mydatatogo,Param=="LE")

dsel$GD2=round(dsel$GD,4)
dsel$vel_abs2=round(dsel$vel_abs,1)
#kd <- with(dsel, MASS::kde2d(GD2,vel_abs2, n = 50))
kd <- with(dsel, MASS::kde2d(GD2,vel_abs2, n = 50,lims=c(round(range(mydatatogo$GD),4),round(range(mydatatogo$vel_abs),1))))


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


pred2=pred2[order(pred2$vel_abs),]
pred2=pred2[order(pred2$GD),]
z1=matrix(pred2$freq,ncol=length(unique(pred2$GD)),byrow=T)
jet.colors2 <- colorRampPalette( c("red","orange","yellow","green") ) 
pal2 <- jet.colors2(100)
#col.ind2 <- col.ind
#col.ind2 <- pal2[col.ind]
#col.ind2[z1<0.5]="#FFFFFF"
z=z1
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
col.ind2 <- cut(z.facet.center, 100)
col.ind2[z.facet.center<0.1]="#FFFFFF"
s1=1:round(max(z.facet.center),0)
col.ind3 <- cut(s1, 100)

png(paste0(dir.out,"/figS1_LE_v062025.png"),unit="cm",width=11,height=11,res=300)#,width=547,height=360

par(mar=c(0,2,0,4))
with(pred2,persp(x=sort(unique(GD)),y=sort(unique(vel_abs)),z=matrix(pred1,ncol=length(unique(GD)),byrow=T),
                 col=pal2[col.ind2],ticktype = "detailed", phi=30, theta=-40,border="grey80",lwd=0.2,
                 xlab="",
                 ylab='',
                 zlab="",
                 cex.lab=0.75,
                 cex.axis=0.75))

r1=rasterFromXYZ(pred2[,c("GD","vel_abs","pred1")])
rX2=r1
rX2[]=sample(1:100,size=ncell(r1),rep=T)

plot(rX2, horizontal=F, legend.only=TRUE,col=pal2,smallplot=c(.85, .87, .4, .8),maxpixels=2000000
     ,axis.args=list(tck=-0.5,at=as.numeric(col.ind3[s1%in%c(1,seq(25,175,by=25),max(s1))]),labels=F)
     ,legend.args=list(text="Density of observations", side=4,font=2, line=1.2, cex=0.65))

plot(rX2, horizontal=F, legend.only=TRUE,col=pal2,smallplot=c(.85, .87, .4, .8),maxpixels=2000000
     ,axis.args=list(tck=F,lwd=0,line=-0.50,at=as.numeric(col.ind3[s1%in%c(1,seq(25,175,by=25),max(s1))]),labels=c(1,seq(25,175,by=25),max(s1)),cex.axis=0.5))

dev.off()

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
r1_df=subset(r1_df,vel_abs<=5.9 & GD<=0.056)
gg=ggplot(data = r1_df, aes(x =GD, y = vel_abs, fill = pred1))  
gg <- gg + geom_tile(color="grey", size=0.1)
#gg<-  gg + theme_minimal()
gg <- gg + scale_fill_viridis(option="E",
                              #name = stringr::str_wrap("Predictions of range shift velocity "),
                              name = bquote("Predictions of range shift velocity "(km.yr^-1)),
                              guide = guide_colourbar(title.position = "left",barwidth = 1,barheight = 15),
                              limits=c(0, 9.033),
                              breaks=c(seq(0,9,by=1)),
                              na.value = "white",
                              direction=-1)+
    #geom_contour(aes(x =GD, y = vel_abs, z = pred1), breaks = seq(0,12,by=2), colour = "grey85", linewidth = 0.2)+
    geom_contour2(aes(z = pred1), breaks = seq(0,12,by=1),colour = "white", linewidth = 0.4)+
    
    metR::geom_text_contour(aes(z = pred1),
                            #colour = "black", size = 4.5, fontface = "bold",
                            breaks = seq(0,12,by=1),
                            stroke = 0.15, stroke.colour = "white", # 'stroke' controls the width of stroke relative to the size of the text
                            skip = 0, # number of contours to skip
                            rotate = FALSE, # horizontal labelling; if TRUE, rotate text following the contour
                            #label.placer = label_placer_fraction(frac = 0.5)
    )

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
    ggtitle("Trailing edge") +labs(tag = '(a)')+
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag.position = c(0.065, 0.975))
gg2 <- gg2+ theme(legend.position = "none")+
    ggtitle("Centroid") +labs(tag = '(b)')+
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag.position = c(0.065, 0.975))
gg3 <- gg3+ theme(legend.position = "none")+
    ggtitle("Leading edge") +labs(tag = '(c)')+
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag.position = c(0.065, 0.975))

png(paste0(dir.out,"/figS2_panelABC_v062025.png"),unit="cm",width=27,height=11,res=300)#,width=547,height=360
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
                              name = stringr::str_wrap("Density of observations"), 
                              guide = guide_colourbar(title.position = "left",barwidth = 1,barheight = 15),
                              limits=c(0.1, 210),
                              breaks=c(1,seq(25,200,by=25)),
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
                              name = stringr::str_wrap("Density of observations"), 
                              guide = guide_colourbar(title.position = "left",barwidth = 1,barheight = 15),
                              limits=c(0.1, 210),
                              breaks=c(1,seq(25,200,by=25)),
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
                              name = stringr::str_wrap("Density of observations"), 
                              guide = guide_colourbar(title.position = "left",barwidth = 1,barheight = 15),
                              limits=c(0.1, 210),
                              breaks=c(1,seq(25,200,by=25)),
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
    labs(tag = '(d)')+
    theme(plot.tag.position = c(0.065, 0.975))
gg2 <- gg2+ theme(legend.position = "none")+
    labs(tag = '(e)')+
    theme(plot.tag.position = c(0.065, 0.975))
gg3 <- gg3+ theme(legend.position = "none")+
    labs(tag = '(f)')+
    theme(plot.tag.position = c(0.065, 0.975))

png(paste0(dir.out,"/figS2_panel DEF.png"),unit="cm",width=27,height=11,res=300)#,width=547,height=360
#gs <- lapply(1:4, function(ii) 
#grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
lay <- rbind(c(1,2,3,4))
gridExtra::grid.arrange(gg1,gg2,gg3,gg_leg, 
                        ncol=4, nrow=1, widths=c(1,1,1,0.2), heights=c(1),layout_matrix = lay)
dev.off()
